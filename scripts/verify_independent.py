"""Independent verification: re-read every raw file with openpyxl/csv (no pandas parsing logic shared with
load_studies.py) and compare every value with processed/<dataset>/expression_matrix.csv.

For each dataset prints: raw values, compared, mismatches, raw samples missing from processed, processed samples
not in raw, duplicate (sample, marker) keys in raw, and what was excluded (controls).
"""
import csv
import math
import os
import re
import statistics
from collections import defaultdict

import openpyxl

BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
PROC = os.path.join(BASE, 'processed')


def num(v):
    if v is None:
        return None
    try:
        x = float(str(v).strip())
        return None if math.isnan(x) else x
    except ValueError:
        return None


def rows_xlsx(path, sheet=None):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet] if sheet else wb.worksheets[0]
    return [list(r) for r in ws.iter_rows(values_only=True)]


def rows_csv(path, delim=','):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        return list(csv.reader(fh, delimiter=delim))


def load_processed(ds):
    d = os.path.join(PROC, ds)
    m = rows_csv(os.path.join(d, 'expression_matrix.csv'))
    cols = m[0][1:]
    mat = {}
    for r in m[1:]:
        for c, v in zip(cols, r[1:]):
            mat[(c, r[0])] = num(v)
    ann = list(csv.DictReader(open(os.path.join(d, 'marker_annotation.csv'), encoding='utf-8')))
    meta = list(csv.DictReader(open(os.path.join(d, 'sample_metadata.csv'), encoding='utf-8')))
    return mat, cols, ann, meta


CTRL = re.compile(r'^(SC\d|IPC|NC\d|Neg|Positive Control|CONTROL|Ctrl|Bridge|Pool|PC\d)', re.I)


def compare(ds, raw, sample_map, marker_map, note=''):
    """raw: dict (raw_sample, raw_marker) -> value (list if duplicated). maps raw ids -> processed ids."""
    mat, cols, ann, meta = load_processed(ds)
    dups = {k: v for k, v in raw.items() if isinstance(v, list)}
    n = mism = miss = 0
    examples = []
    raw_samples = set(s for s, _ in raw)
    excluded = sorted(s for s in raw_samples if sample_map(s) is None)
    for (s, mk), v in raw.items():
        ps, pm = sample_map(s), marker_map(mk)
        if ps is None:
            continue
        if isinstance(v, list):
            v = statistics.mean([x for x in v if x is not None]) if any(x is not None for x in v) else None
        n += 1
        if (ps, pm) not in mat:
            miss += 1
            if len(examples) < 3:
                examples.append(('missing', s, mk, ps, pm))
            continue
        p = mat[(ps, pm)]
        if (v is None) != (p is None) or (v is not None and abs(v - p) > 1e-9 * max(1, abs(v))):
            mism += 1
            if len(examples) < 3:
                examples.append(('diff', s, mk, v, p))
    mapped = set(sample_map(s) for s in raw_samples) - {None}
    extra = set(cols) - mapped
    meta_ok = [r['sample_id'] for r in meta] == cols
    n_markers_raw = len(set(marker_map(mk) for _, mk in raw))
    rows_proc = len(set(m for _, m in mat))
    print(f'{ds:34s} values={n:7d} mismatch={mism} missing={miss} | samples raw={len(raw_samples)} excluded={len(excluded)} '
          f'processed={len(cols)} extra_in_processed={len(extra)} | markers raw={n_markers_raw} processed={rows_proc} | '
          f'meta_order_ok={meta_ok} | raw_dup_keys={len(dups)} {note}')
    if excluded:
        print('     excluded samples:', excluded[:12], '...' if len(excluded) > 12 else '')
    for e in examples:
        print('     e.g.', e)
    if dups:
        print('     duplicate raw keys e.g.', list(dups.items())[:2])
    return dict(ds=ds, values=n, mismatch=mism, missing=miss, extra=len(extra), meta_ok=meta_ok, dups=len(dups))


def put(raw, k, v):
    if k in raw:
        prev = raw[k]
        raw[k] = (prev if isinstance(prev, list) else [prev]) + [v]
    else:
        raw[k] = v


def ann_map(ds, key):
    """map raw marker identifier -> processed row label using marker_annotation column `key`"""
    _, _, ann, _ = load_processed(ds)
    return {a[key]: a['marker_row_label'] for a in ann}


# ---------- export-format files
def export(path):
    R = rows_xlsx(path)
    first = [str(r[0]).strip() if r[0] is not None else '' for r in R]
    ia = first.index('Assay')
    iu = [i for i, x in enumerate(first) if x.startswith('Uniprot')][0]
    io = [i for i, r in enumerate(R) if any(str(x).startswith('OID') for x in r[1:6] if x)][0]
    mcols = [j for j in range(1, len(R[ia])) if R[iu][j] and re.match(r'^[A-Z0-9]', str(R[iu][j])) and str(R[iu][j]).strip() not in ('Det Ctrl', 'Inc Ctrl 2', 'Inc Ctrl 1', 'Ext Ctrl')
             and str(R[ia][j]).strip() not in ('Plate ID', 'QC Warning', 'QC Deviation from median')]
    raw = {}
    for i in range(io + 1, len(R)):
        s = first[i]
        if not s or s in ('LOD', 'MaxLOD', 'PlateLOD', 'Missing Data freq.', 'Normalization', 'nan') or s.startswith('Patient id'):
            continue
        vals = [num(R[i][j]) for j in mcols]
        if all(v is None for v in vals):
            continue
        for j, v in zip(mcols, vals):
            put(raw, (s, str(R[io][j]).strip()), v)
    return raw


def v010():
    ds = 'row010_Meng2026_GastricCa'
    raw = export(os.path.join(BASE, ds, '03_Olink_LC-P20240914054_DYW_240927_NPX.xlsx'))
    _, _, _, meta = load_processed(ds)
    lab = {r['lab_sample_id']: r['sample_id'] for r in meta}
    om = ann_map(ds, 'olink_id')
    return compare(ds, raw, lambda s: None if CTRL.match(s) else lab.get(s, '??' + s), lambda m: om.get(m, '??' + m))


def v_export(ds, f, sample_ok=lambda s: True):
    raw = export(os.path.join(BASE, ds, f))
    om = ann_map(ds, 'olink_id')
    return compare(ds, raw, lambda s: s if (not CTRL.match(s) and sample_ok(s)) else None, lambda m: om.get(m, '??' + m))


def v027():
    ds = 'row027_Tikkanen2025_SolidTumors'
    R = rows_xlsx(os.path.join(BASE, ds, '41598_2025_31229_MOESM10_ESM.xlsx'))
    raw = {}
    for r in R[1:]:
        for s, v in zip(R[0][1:], r[1:]):
            put(raw, (s, r[0]), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v032():
    ds = 'row032_Gao2025_NSCLC'
    R = rows_csv(os.path.join(BASE, ds, 'olink.txt'), '\t')
    clin = {'Response', 'relapse', 'PFS', 'dead', 'OS', 'Age', 'Gender', 'Pathological_type', 'Smoking', 'Chemotherapy', 'Stage'}
    raw = {}
    for r in R[1:]:
        for h, v in zip(R[0][1:], r[1:]):
            if h not in clin:
                put(raw, (r[0], h), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v033():
    ds = 'row033_Xiao2025_TNBC'
    R = rows_xlsx(os.path.join(BASE, ds, 'TableS1-S2_media-9.xlsx'), 'Table S2')
    hdr = R[1]
    raw = {}
    for r in R[2:]:
        if r[0] is None:
            continue
        for h, v in zip(hdr[1:], r[1:]):
            if h is not None:
                put(raw, (str(r[0]), str(h)), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v040(ds, f):
    R = rows_csv(os.path.join(BASE, 'row040_Plummer2025_TIMEPOINT', f))
    raw = {}
    for r in R[1:]:
        for h, v in zip(R[0][1:], r[1:]):
            if not re.match(r'Plate ID|QC', h):
                put(raw, (r[0], h), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v044():
    ds = 'row044_Braun2025_RCCvaccine'
    R = rows_xlsx(os.path.join(BASE, ds, '41586_2024_8507_MOESM5_ESM.xlsx'))
    raw = {}
    for r in R[1:]:
        for h, v in zip(R[0][1:4], r[1:4]):
            put(raw, (h, r[0]), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v_long(ds, files, sid, mcol, vcol, mapcol='original_name', sample_filter=lambda r: True):
    raw = {}
    for f in files:
        with open(os.path.join(BASE, ds.split('_CITN')[0] if 'CITN' in ds else ds, f), newline='', encoding='utf-8-sig') as fh:
            for r in csv.DictReader(fh):
                if sample_filter(r):
                    put(raw, (sid(r), r[mcol]), num(r[vcol]))
    om = ann_map(ds, mapcol)
    return raw, om


def v052():
    ds = 'row052_Tsimberidou2024_AMADEUS'
    raw, om = v_long(ds, ['translational-data_AMADEUS_primarycohort_olink.csv'], lambda r: r['subject.id'] + '_' + r['timepoint.id'], 'protein.name', 'npx')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v_parker(ds, files):
    raw, om = v_long(ds, files, lambda r: r['manuscript.id'] + '_' + r['timepoint.id'], 'protein.name', 'npx')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v070():
    ds = 'row070_Nunez2023_irAE'
    R = rows_xlsx(os.path.join(BASE, ds, 'Proteomics_discovery_log2.xlsx'))
    hdr = R[0]
    seen = defaultdict(int)
    keys = []
    for r in R[1:]:
        seen[(r[0], r[1])] += 1
        keys.append((r[0], r[1], seen[(r[0], r[1])]))
    twice = set(p for (p, t), c in seen.items() if c > 1)
    raw = {}
    for (p, t, s), r in zip(keys, R[1:]):
        sid = f'{p}_s{s}_T{t}' if p in twice else f'{p}_T{t}'
        for h, v in zip(hdr[3:], r[3:]):
            put(raw, (sid, h), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m), note=f'(patients listed twice: {sorted(twice)})')


def v072():
    ds = 'row072_Padron2022_PRINCE'
    raw, cnt = {}, defaultdict(int)
    with open(os.path.join(BASE, ds, 'Olink_NatureMed_Olink_select_proteins.csv'), newline='', encoding='utf-8-sig') as fh:
        for r in csv.DictReader(fh):
            cnt[(r['sample.id'], r['paper.name'])] += 1
            put(raw, (r['sample.id'], f"{r['paper.name']}#{cnt[(r['sample.id'], r['paper.name'])]}"), num(r['value']))
    _, _, ann, _ = load_processed(ds)
    om = {f"{a['paper_name']}#{a['measurement']}": a['marker_row_label'] for a in ann}
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v078():
    ds = 'row078_Barroso-Sousa2021_TNBC'
    raw = {}
    for s, col in [('Plasma biomarker assoc with PFS', 'baseline_hazardratio'), ('Plasma biomarker change', 'change_hazardratio')]:
        R = rows_xlsx(os.path.join(BASE, ds, '41523_2021_287_MOESM1_ESM.xlsx'), s)
        for r in R[3:]:
            if r[0]:
                put(raw, (col, r[0]), num(r[1]))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v086(ds, f):
    raw, om = v_long(ds, [f], lambda r: r['SampleID'], 'OlinkID', 'NPX', mapcol='olink_id')
    return compare(ds, raw, lambda s: s if re.match(r'^1[03]_\d+\.', s) else None, lambda m: om.get(m, '??' + m))


def v103():
    ds = 'row103_Yu2025_LGL'
    R = rows_xlsx(os.path.join(BASE, ds, '41392_2025_2164_MOESM1_ESM.xlsx'))
    raw = {}
    for r in R[1:]:
        for h, v in zip(R[0][1:], r[1:]):
            put(raw, (h, r[0]), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


def v096():
    ds = 'row096_Fedak2026_MG'
    R = rows_xlsx(os.path.join(BASE, ds, '41591_2025_4170_MOESM11_ESM.xlsx'))
    vals = defaultdict(list)
    j = 0
    while j < len(R[0]):
        if R[0][j] and str(R[0][j]).strip():
            k = j + 1
            while k < len(R[0]) and R[0][k] and str(R[0][k]).strip():
                for r in R[1:]:
                    if num(r[k]) is not None:
                        vals[(str(R[0][k]).replace('\xa0', ' ').strip(), str(R[0][j]).strip())].append(num(r[k]))
                k += 1
            j = k
        else:
            j += 1
    S = rows_xlsx(os.path.join(BASE, ds, '41591_2025_4170_MOESM6_ESM.xlsx'), '4N')
    for jj in range(1, len(S[0])):
        for r in S[1:]:
            if num(r[jj]) is not None:
                vals[('D1 vs SCR', str(S[0][jj]).strip())].append(num(r[jj]))
    raw = {k: statistics.median(v) for k, v in vals.items()}
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m), note='(medians recomputed)')


def v062():
    ds = 'row062_Coffey2024_Myeloma'
    R = rows_xlsx(os.path.join(BASE, ds, 'jitc-2023-008110supp002_data_supplement.xlsx'), 'Supplemental table 3')
    hdr = R[2]
    raw = {}
    for r in R[3:]:
        if r[0] is None:
            continue
        sid = 'P%02d_%s' % (int(r[0]), str(r[1]).replace(' ', ''))
        for h, v in zip(hdr[2:], r[2:]):
            if h is not None:
                put(raw, (sid, str(h)), num(v))
    om = ann_map(ds, 'original_name')
    return compare(ds, raw, lambda s: s, lambda m: om.get(m, '??' + m))


if __name__ == '__main__':
    res = [v010(), v062(),
           v_export('row012_Johnson2026_ccRCC', '41467_2026_73336_MOESM5_ESM.xlsx'),
           v_export('row024_Shan2026_OvarianCAR-T', '01_Olink_LC-P20241107004_DYW_241216_NPX.xlsx'),
           v027(), v032(), v033(),
           v040('row040_Plummer2025_TIMEPOINT', 'TIMEPOINT_OLINK_passedQC_C1D1.csv'),
           v040('row040_Plummer2025_OUTREACH', 'OUTREACH_C1D1_inflammation_pride.csv'),
           v044(), v052(), v070(), v072(),
           v_export('row073_Liu2022_TNBC', '41467_2022_30569_MOESM5_ESM.xlsx', sample_ok=lambda s: bool(re.match(r'^\d{3}-\d$', s))),
           v078(),
           v086('row086_Glass2024_CTCL_CITN10', 'olink_citn_10.csv'),
           v086('row086_Glass2024_CTCL_CITN13', 'olink_citn_13.csv'),
           v_parker('row093_OReilly2026_REVOLUTION', ['translational-data_REV_cohortA_olink_public.csv', 'translational-data_REV_cohortB_olink_public.csv']),
           v096(), v103(),
           v_parker('row105_Galsky2025_PORTER', ['translational-data_PORTER_olink_public.csv'])]
    bad = [r for r in res if r['mismatch'] or r['missing'] or r['extra'] or not r['meta_ok']]
    print('\nTOTAL values compared:', sum(r['values'] for r in res), '| datasets with problems:', [r['ds'] for r in bad])
