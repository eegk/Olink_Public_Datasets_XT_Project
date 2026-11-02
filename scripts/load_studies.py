"""Parse each downloaded Olink dataset into a markers x samples matrix, sample metadata and marker info.

Each loader returns a dict:
  mat     DataFrame, rows = marker key, columns = sample id
  meta    DataFrame, index = sample id (same order as mat columns)
  markers DataFrame, index = marker key: original_name, uniprot, panel, olink_id (+ extras)
  kind    'per-sample' or 'summary'
  notes   free text describing scale, caveats
"""
import os
import re

import numpy as np
import pandas as pd

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))


def P(*a):
    return os.path.join(BASE, *a)


ACC = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$')
FOOT = {'LOD', 'MaxLOD', 'PlateLOD', 'Missing Data freq.', 'Normalization', 'Plate LOD', 'Max LOD'}


def is_ctrl(s):
    return bool(re.match(r'^(SC\d|IPC|NC\d|Neg|Positive Control|CONTROL|Ctrl|Bridge|Pool|PC\d)', str(s), re.I))


def read_export(path):
    """Olink NPX Manager / Signature wide export (Panel/Assay/Uniprot ID/OlinkID header rows, one row per sample)."""
    d = pd.read_excel(path, header=None)
    c0 = d.iloc[:, 0].astype(str).str.strip()
    ia = c0[c0 == 'Assay'].index[0]
    iu = c0[c0.str.startswith('Uniprot')].index[0]
    ip = c0[c0 == 'Panel'].index[0]
    io = [i for i in d.index if any(str(x).startswith('OID') for x in d.iloc[i, 1:6])][0]
    def is_acc(v):
        return isinstance(v, str) and all(ACC.match(p) for p in re.split(r'[,_;\s]+', v.strip()) if p)
    cols = [j for j in d.columns[1:] if is_acc(d.iat[iu, j])]
    rows = []
    for i in d.index:
        v = str(d.iat[i, 0]).strip()
        if i <= io or v in FOOT or v == 'nan' or not re.match(r'^[A-Za-z0-9]', v):
            continue
        if pd.to_numeric(d.iloc[i, cols], errors='coerce').notna().sum() > 0:
            rows.append(i)
    mk = pd.DataFrame({
        'original_name': [str(d.iat[ia, j]).strip() for j in cols],
        'uniprot': [str(d.iat[iu, j]).strip() for j in cols],
        'panel': [str(d.iat[ip, j]).strip() for j in cols],
        'olink_id': [str(d.iat[io, j]).strip() for j in cols]})
    vals = d.loc[rows, cols].apply(pd.to_numeric, errors='coerce')
    vals.index = [str(d.iat[i, 0]).strip() for i in rows]
    mat = vals.T
    mat.index = range(len(mk))
    qc = None
    for j in d.columns:
        if str(d.iat[ia, j]).strip() == 'QC Warning':
            qc = pd.Series([d.iat[i, j] for i in rows], index=vals.index)
    if qc is not None:
        qc = qc[~qc.index.duplicated()]
    return mk, mat, qc


def finalize(mk, mat):
    """Unique marker keys: assay name, suffixed with panel (then OlinkID) only when duplicated."""
    key = mk['original_name'].copy()
    pan = mk['panel'].astype(str).str.replace(r'Olink (Target 96 )?', '', regex=True).str.replace(r'\(.*\)', '', regex=True).str.strip()
    dup = key.duplicated(keep=False)
    key[dup] = key[dup] + '|' + pan[dup]
    dup = key.duplicated(keep=False)
    key[dup] = key[dup] + '|' + mk.loc[dup, 'olink_id'].astype(str)
    mk = mk.copy()
    mk.index = key.values
    mat = mat.copy()
    mat.index = key.values
    return mk, mat


def simple_markers(names, panel, uniprot=None):
    return pd.DataFrame({'original_name': list(names), 'uniprot': uniprot if uniprot is not None else np.nan,
                         'panel': panel, 'olink_id': np.nan}, index=list(names))


def snake(cols):
    return [re.sub(r'[^0-9a-zA-Z]+', '_', str(c)).strip('_').lower() for c in cols]


# ---------------------------------------------------------------- studies
def s010():
    mk, mat, qc = read_export(P('row010_Meng2026_GastricCa', '03_Olink_LC-P20240914054_DYW_240927_NPX.xlsx'))
    mat = mat[[c for c in mat.columns if not is_ctrl(c)]]
    info = pd.read_excel(P('row010_Meng2026_GastricCa', 'Sample_Info.xlsx'))
    info['SampleID'] = info['SampleID'].replace({'DYX': 'DYX_1'})
    m = info.set_index('Sample').reindex(mat.columns)
    conds = [c for c in m.columns if c.startswith('COND')]

    def resp(r):
        for c in conds:
            v = str(r[c])
            if v.startswith('NR_'):
                return 'NR'
            if v.startswith('R_'):
                return 'R'
        return np.nan

    meta = pd.DataFrame({'sample_id': m['SampleID'].values, 'lab_sample_id': mat.columns})
    meta['patient_id'] = meta['sample_id'].str.split('_').str[0]
    meta['timepoint_index'] = pd.to_numeric(meta['sample_id'].str.split('_').str[1], errors='coerce')
    meta['timepoint_label'] = meta['timepoint_index'].map({1: 'baseline', 2: 'T1 (after 2 cycles)', 3: 'T2 (after 4 cycles)',
                                                           4: 'T3 (after 6 cycles)'})
    r = pd.Series(m.apply(resp, axis=1).values, index=meta['patient_id'])
    meta['response'] = meta['patient_id'].map(r.dropna().groupby(level=0).first())
    if qc is not None:
        meta['qc_warning'] = qc.reindex(mat.columns).values
    for c in conds:
        meta['author_' + c.lower()] = m[c].values
    meta.index = meta['sample_id']
    mat.columns = meta.index
    mk, mat = finalize(mk, mat)
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes="NPX (log2), Target 96 Immuno-Oncology. Columns renamed to the authors' IDs (patient_timepoint). "
                      "Response R/NR taken from the authors' COND grouping columns (kept as author_cond1..6). "
                      "timepoint_index 1-4 = baseline, then every two treatment cycles (paper: baseline, T1, T2, T3); "
                      "timepoint_label added accordingly. The deposit holds 97 samples from 35 patient codes, whereas the paper "
                      "reports 31 patients (31/31/20/4 samples at baseline/T1/T2/T3); the extra patients cannot be identified "
                      "from the deposit. Sample controls SC1/SC2 removed.")


def s012():
    mk, mat, qc = read_export(P('row012_Johnson2026_ccRCC', '41467_2026_73336_MOESM5_ESM.xlsx'))
    mat = mat[[c for c in mat.columns if not is_ctrl(c)]]
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    meta['patient_id'] = meta.index.str.extract(r'^(\d+-\d+)', expand=False).values
    meta['visit_code'] = meta.index.str.extract(r'^\d+-\d+(.*)$', expand=False).values
    meta['treatment_block'] = pd.to_numeric(meta['visit_code'].str.extract(r'B(\d+)', expand=False), errors='coerce')
    meta['week_in_block'] = pd.to_numeric(meta['visit_code'].str.extract(r'W(\d+)', expand=False), errors='coerce')
    mk, mat = finalize(mk, mat)
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes='NPX (log2), three Target 96 panels (Immuno-Oncology, Oncology II, Cardiometabolic); an assay measured on more '
                      'than one panel is kept as separate rows (suffix |panel). visit_code e.g. B1W1 = treatment block 1, week 1 '
                      '(block 1 = pembrolizumab alone; block 2 = pembrolizumab + high-dose IL-2, per the paper). The paper groups '
                      'patients as EXT/INT/NON responders (86 samples from 26 patients); that grouping is not in the deposited file '
                      '(87 samples, 27 patient codes here).')


def s024():
    mk, mat, qc = read_export(P('row024_Shan2026_OvarianCAR-T', '01_Olink_LC-P20241107004_DYW_241216_NPX.xlsx'))
    mat = mat[[c for c in mat.columns if not is_ctrl(c)]]
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    if qc is not None:
        meta['qc_warning'] = qc.reindex(meta.index).values
    mk, mat = finalize(mk, mat)
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes='NPX (log2), Target 96 Oncology II (not IO). The figshare record describes 9 patients in three groups '
                      '(PD n=2; SD <3 months n=3; SD >=3 months n=4) but gives no sample-to-patient/group/timepoint key, so the '
                      'metadata holds only the 24 lab sample IDs.')


def s027():
    d = pd.read_excel(P('row027_Tikkanen2025_SolidTumors', '41598_2025_31229_MOESM10_ESM.xlsx')).set_index('Assay')
    mat = d.apply(pd.to_numeric, errors='coerce')
    m = pd.read_excel(P('row027_Tikkanen2025_SolidTumors', '41598_2025_31229_MOESM2_ESM.xlsx')).set_index('SampleID')
    m.columns = ['ecog', 'sex', 'tumor_type', 'pd_l1_therapy_line', 'orr']
    meta = m.reindex(mat.columns)
    meta.insert(0, 'sample_id', mat.columns)
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology'), kind='per-sample',
                notes='NPX (log2), baseline plasma. orr as coded in the source (1 = objective response).')


def s032():
    d = pd.read_csv(P('row032_Gao2025_NSCLC', 'olink.txt'), sep='\t').set_index('SampleID')
    clin = ['Response', 'relapse', 'PFS', 'dead', 'OS', 'Age', 'Gender', 'Pathological_type', 'Smoking', 'Chemotherapy', 'Stage']
    clin = [c for c in clin if c in d.columns]
    mat = d.drop(columns=clin).T.apply(pd.to_numeric, errors='coerce')
    meta = d[clin].copy()
    meta.columns = snake(meta.columns)
    meta.insert(0, 'sample_id', meta.index)
    meta['id_series'] = np.where(meta.index.str.startswith('CT'), 'CT-series', '24R-series')
    meta['cohort_inferred'] = np.where(meta.index.str.startswith('CT'), 'Cohort 1 (study cohort)', 'Cohort 2 (validation, pre-treatment)')
    meta = meta.reindex(mat.columns)
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology'), kind='per-sample',
                notes='Values look centred/scaled (negative values present) rather than raw NPX. Clinical/outcome columns '
                      '(response, relapse, PFS, death, OS, age, gender, histology, smoking, chemotherapy, stage) were embedded in '
                      'olink.txt and moved to the metadata; coding as in source. cohort_inferred: the 30 "24R" samples match the '
                      "paper's validation cohort (30 pre-treatment samples, 2023-24) and the 29 \"CT\" samples the study cohort "
                      '(paper: 34 patients, 58 pre/post samples - only 29 are in the deposit, timepoint not stated). This mapping is '
                      'inferred from the counts and ID formats, not stated in the file. Response/demographics are only filled for '
                      'the CT series; PFS/OS are present for all 59.')


def s033():
    f = P('row033_Xiao2025_TNBC', 'TableS1-S2_media-9.xlsx')
    d = pd.read_excel(f, sheet_name='Table S2', header=1).set_index('Sample ID')
    mat = d.T.apply(pd.to_numeric, errors='coerce')
    s1 = pd.read_excel(f, sheet_name='Table S1', header=1).set_index('Patient ID')
    s1.columns = snake(s1.columns)
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    parts = [str(x).split('_') for x in meta.index]
    meta['patient_id'] = [p[0] for p in parts]
    meta['group'] = [p[1] if len(p) > 1 else np.nan for p in parts]
    meta['timepoint'] = ['_'.join(p[2:]) for p in parts]
    meta = meta.join(s1.drop(columns=['group'], errors='ignore'), on='patient_id')
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology'), kind='per-sample',
                notes='Olink IO values from Table S2; clinical data from Table S1. Groups per the paper: A = neoadjuvant (101 '
                      'camrelizumab + chemotherapy and 26 chemotherapy-only - see the treatment column), B = adjuvant '
                      'atezolizumab + chemotherapy, C = advanced disease on camrelizumab-based therapy. Response is coded '
                      'Response/nonResponse (group A) or RECIST category (group C); "/" = not assessed. Timepoints Pre/Post '
                      '(group C also Post_1-3).')


def _s040(file, metafile, cohort):
    d = pd.read_csv(P('row040_Plummer2025_TIMEPOINT', file)).set_index('Patient')
    qc = d['QCWarning'] if 'QCWarning' in d.columns else None
    d = d[[c for c in d.columns if not re.match(r'Plate ID|QC', c)]]
    mat = d.T.apply(pd.to_numeric, errors='coerce')
    m = pd.read_csv(P('row040_Plummer2025_TIMEPOINT', metafile)).set_index('Patient ID')
    m.columns = snake(m.columns)
    meta = m.reindex(mat.columns)
    meta.insert(0, 'sample_id', mat.columns)
    meta.insert(1, 'cohort', cohort)
    if qc is not None:
        meta['qc_warning'] = qc.reindex(mat.columns).values
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Inflammation'), kind='per-sample',
                notes=f'{cohort} cohort, C1D1 (baseline) plasma, NPX. Olink Target 96 Inflammation panel (not IO).')


def s040a():
    return _s040('TIMEPOINT_OLINK_passedQC_C1D1.csv', 'TIMEPOINT_Olink_metadata.csv', 'TIMEPOINT')


def s040b():
    return _s040('OUTREACH_C1D1_inflammation_pride.csv', 'OUTREACH_Olink_metadata.csv', 'OUTREACH')


def s044():
    d = pd.read_excel(P('row044_Braun2025_RCCvaccine', '41586_2024_8507_MOESM5_ESM.xlsx')).set_index('Cytokine')
    mat = d[['Week_00_median', 'Week_03_median', 'Week_20_median']].apply(pd.to_numeric, errors='coerce')
    mk = simple_markers(mat.index, 'Target 96 Immuno-Oncology')
    mk['pvalue_w00_w03'] = d['Pvalue_W00_W03'].values
    mk['pvalue_w00_w20'] = d['Pvalue_W00_W20'].values
    meta = pd.DataFrame({'sample_id': mat.columns, 'timepoint': ['Week 0', 'Week 3', 'Week 20'],
                         'statistic': 'median NPX across patients'}).set_index('sample_id', drop=False)
    return dict(mat=mat, meta=meta, markers=mk, kind='summary',
                notes='SUMMARY ONLY: median NPX across patients per timepoint (no per-sample data released). Columns are timepoints.')


def s052():
    d = pd.read_csv(P('row052_Tsimberidou2024_AMADEUS', 'translational-data_AMADEUS_primarycohort_olink.csv'))
    d['sample_id'] = d['subject.id'] + '_' + d['timepoint.id']
    mat = d.pivot_table(index='protein.name', columns='sample_id', values='npx', aggfunc='mean')
    sm = d.drop_duplicates('sample_id').set_index('sample_id')[['subject.id', 'manuscript.id', 'arm', 'treatment', 'tumor.type', 'timepoint.id']]
    subj = pd.read_csv(P('row052_Tsimberidou2024_AMADEUS', 'clinical-data_AMADEUS_primarycohort_subject.csv')).set_index('subject.id')
    meta = sm.join(subj.drop(columns=['manuscript.id', 'arm', 'treatment', 'tumor.type'], errors='ignore'), on='subject.id')
    meta.insert(0, 'sample_id', meta.index)
    meta = meta.reindex(mat.columns)
    meta.columns = snake(meta.columns)
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology'), kind='per-sample',
                notes='NPX; long file pivoted to markers x (subject_timepoint). Source protein names are UniProt-style short names '
                      '(e.g. AGRG1 = ADGRG1). Outcomes (BOR, PFS, OS, responder flags, CD8 converter) joined from the subject table.')


def s070():
    d = pd.read_excel(P('row070_Nunez2023_irAE', 'Proteomics_discovery_log2.xlsx'))
    # 3 patients (IMIT104, IMIT76, IMIT89) appear twice: identical T1/T2 rows but different T3/T4 rows -> keep both series
    d['series'] = d.groupby(['Patient_ID', 'Timepoint']).cumcount() + 1
    twice = d.loc[d['series'] > 1, 'Patient_ID'].unique()
    d['sample_id'] = [f'{p}_s{s}_T{t}' if p in twice else f'{p}_T{t}' for p, s, t in zip(d['Patient_ID'], d['series'], d['Timepoint'])]
    md = d[['sample_id', 'Patient_ID', 'Timepoint', 'Condition', 'series']].set_index('sample_id')
    md.columns = ['patient_id', 'timepoint', 'condition', 'series']
    md['patient_listed_twice'] = md['patient_id'].isin(twice)
    md['timepoint_label'] = md['timepoint'].map({1: 'T1 pre-therapy', 2: 'T2 1-2 weeks after start',
                                                 3: 'T3 ~1 month before first irAE (matched time for no-irAE)',
                                                 4: 'T4 at first irAE (matched time for no-irAE)'})
    d = d.drop(columns=['series'])
    md['response'] = md['condition'].map({'RT': 'R', 'RNT': 'R', 'NRNT': 'NR'})
    md['irae_toxicity'] = md['condition'].map({'RT': 'Yes', 'RNT': 'No', 'NRNT': 'No'})
    md.insert(0, 'sample_id', md.index)
    x = d.drop(columns=['Patient_ID', 'Timepoint', 'Condition']).set_index('sample_id')
    mat = x.T.apply(pd.to_numeric, errors='coerce')
    return dict(mat=mat, meta=md.reindex(mat.columns), markers=simple_markers(mat.index, 'Target 96 Inflammation'), kind='per-sample',
                notes='Discovery cohort, 172 SERUM samples, log2 Olink values (Inflammation panel, not IO). Condition codes, '
                      'confirmed in the paper: RT = "Tox R" (responder with irAE), RNT = "No Tox R", NRNT = "No Tox NR". '
                      'Timepoints per the paper: T1 pre-therapy, T2 1-2 weeks after start, T3 one month before first irAE, '
                      'T4 at first irAE (matched times for patients without irAE). '
                      'Patients IMIT104, IMIT76 and IMIT89 are listed twice in the source (identical T1/T2 values, different T3/T4); '
                      'both series are kept (sample ids *_s1_* / *_s2_*, patient_listed_twice = True). '
                      'The validation-cohort ratio file (4 proteins, T2/T1) is not included.')


PRINCE_MAP = {
    'Angiopoietin-1 Receptor (TIE2)': 'TIE2', 'Arginase 1 (ARG1)': 'ARG1', 'CD205': 'LY75', 'CXCL10': 'CXCL10', 'CXCL9': 'CXCL9',
    'Carbonic anhydrase (CAH9)': 'CAIX', 'Caspase 8': 'CASP-8', 'Chemokine (C-C motif) ligand 19 (CCL19)': 'CCL19',
    'Chemokine (C-C motif) ligand 23 (CCL23)': 'CCL23', 'Chemokine (C-C motif) ligand 7 (CCL7)': 'MCP-3', 'Galectin 9 (Gal-9)': 'Gal-9',
    'Hematopoietic cell-specific Lyn Substrate 1 (HCLS1)': 'HCLS1', 'Hepatocyte growth factor (HGF)': 'HGF', 'IFN-γ': 'IFN-gamma',
    'Interleukin 10 (IL-10)': 'IL10', 'Interleukin 13 (IL-13)': 'IL13', 'Interleukin 2 (IL-2)': 'IL2', 'Interleukin 4 (IL-4)': 'IL4',
    'Interleukin 6 (IL-6)': 'IL6', 'Interleukin 8 (IL-8)': 'IL8', 'LAMP3': 'LAMP3', 'Matrix metallopeptidase 12 (MMP12)': 'MMP12',
    'Matrix metallopeptidase 7 (MMP7)': 'MMP7', 'Nitric oxide synthase 3 (NOS3)': 'NOS3', 'Soluble CD137': 'TNFRSF9',
    'Soluble CD27': 'CD27', 'Soluble CD83': 'CD83', 'Soluble ICOSL': 'ICOSLG', 'Soluble PD-1': 'PDCD1',
    'Vascular endothelial growth factor A (VEGF A)': 'VEGFA'}


def s072():
    p = pd.read_csv(P('row072_Padron2022_PRINCE', 'Olink_NatureMed_Olink_select_proteins.csv'))
    m = pd.read_csv(P('row072_Padron2022_PRINCE', 'Olink_NatureMed_Olink_metadata.csv')).set_index('sample.id')
    c = pd.read_csv(P('row072_Padron2022_PRINCE', 'Clinical_PICI0002_ph2_clinical.csv')).set_index('Deidentified.ID')
    p['marker'] = p['paper.name'].map(PRINCE_MAP).fillna(p['paper.name'])
    # IL-10, IL-6 and LAMP3 have two values per sample (adjacent rows; most likely the same protein on two Olink
    # panels, which the release does not label). They are NOT averaged: kept as two rows, numbered by file order.
    p['measurement'] = p.groupby(['sample.id', 'paper.name']).cumcount() + 1
    twice = set(p.loc[p['measurement'] > 1, 'marker'])
    p['key'] = [f'{m}|m{k}' if m in twice else m for m, k in zip(p['marker'], p['measurement'])]
    mat = p.pivot(index='key', columns='sample.id', values='value')
    inv = {v: k for k, v in PRINCE_MAP.items()}
    base = [k.split('|')[0] for k in mat.index]
    mk = pd.DataFrame({'original_name': base, 'uniprot': np.nan,
                       'panel': [f'measurement {k.split("|m")[1]} of 2 (panel not labelled in release)' if '|m' in k else 'Olink (selected proteins)' for k in mat.index],
                       'olink_id': np.nan, 'paper_name': [inv.get(b, b) for b in base],
                       'measurement': [int(k.split('|m')[1]) if '|m' in k else 1 for k in mat.index]}, index=mat.index)
    meta = m.reindex(mat.columns).join(c, on='Deidentified.ID').drop(columns=['ms'])
    meta.insert(0, 'sample_id', meta.index)
    meta.columns = snake(meta.columns)
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes='Only the 30 proteins released by the authors (not the full panel); values as provided. IL-10, IL-6 and LAMP3 have '
                      'two values per sample in the release (probably two Olink panels; not labelled) - kept as separate rows '
                      '"|measurement 1 of 2" / "|measurement 2 of 2" in file order, not averaged (their agreement is weak for '
                      'IL-10 and LAMP3, r ~0.3-0.5). Descriptive protein names '
                      'mapped to Olink assay names (paper_name kept in markers file). Clinical/outcomes (arm, BOR, PFS, OS) joined by patient.')


def s073():
    mk, mat, qc = read_export(P('row073_Liu2022_TNBC', '41467_2022_30569_MOESM5_ESM.xlsx'))
    mat = mat[[c for c in mat.columns if re.match(r'^\d{3}-\d$', c)]]
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    meta['patient_id'] = meta.index.str.split('-').str[0]
    meta['timepoint'] = meta.index.str.split('-').str[1].map({'1': 'pre-treatment', '2': 'post-treatment'})
    mk, mat = finalize(mk, mat)
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes='NPX (log2), Target 96 IO. Suffix -1 = pre-treatment, -2 = post-treatment (per the file header). '
                      'Response data are not in this supplementary file.')


def s078():
    out = []
    f = P('row078_Barroso-Sousa2021_TNBC', '41523_2021_287_MOESM1_ESM.xlsx')
    for s, lab in [('Plasma biomarker assoc with PFS', 'baseline'), ('Plasma biomarker change', 'change')]:
        d = pd.read_excel(f, sheet_name=s, header=2)
        d.columns = [str(c).lower() for c in d.columns]
        d = d.rename(columns={'biomarker': 'marker'}).dropna(subset=['marker']).set_index('marker')
        out.append(d.add_prefix(lab + '_'))
    res = pd.concat(out, axis=1)
    mat = res[['baseline_hazardratio', 'change_hazardratio']].apply(pd.to_numeric, errors='coerce')
    mk = simple_markers(res.index, 'Target 96 Immuno-Oncology').join(res)
    meta = pd.DataFrame({'sample_id': mat.columns, 'statistic': ['PFS hazard ratio, baseline level', 'PFS hazard ratio, on-treatment change']}).set_index('sample_id', drop=False)
    return dict(mat=mat, meta=meta, markers=mk, kind='summary',
                notes='SUMMARY ONLY: Cox hazard ratios for PFS per marker (baseline level and change); CIs and p-values in markers file. '
                      'No per-sample values released.')


def _s086(f, study):
    d = pd.read_csv(P('row086_Glass2024_CTCL', f))
    d = d[d['SampleID'].astype(str).str.match(r'^1[03]_\d+\.')]
    d['key'] = d['Assay'].astype(str) + '|' + d['OlinkID'].astype(str)
    mat = d.pivot_table(index='key', columns='SampleID', values='NPX', aggfunc='mean')
    mk = d.drop_duplicates('key').set_index('key').reindex(mat.index)[['Assay', 'UniProt', 'Panel', 'OlinkID']]
    mk.columns = ['original_name', 'uniprot', 'panel', 'olink_id']
    a = mk['original_name']
    nk = a.where(~a.duplicated(keep=False), a + '|' + mk['panel'].astype(str).str.replace('Olink ', ''))
    nk = nk.where(~nk.duplicated(keep=False), nk + '|' + mk['olink_id'].astype(str))
    mk.index = nk.values
    mat.index = nk.values
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    meta['subject'] = meta.index.str.split('.').str[0]
    meta['timepoint'] = meta.index.str.split('.').str[1]
    pm = pd.read_csv(P('row086_Glass2024_CTCL', 'patient_metadata.csv')).set_index('subject')
    meta = meta.join(pm, on='subject')
    qc = d.groupby('SampleID')['QC_Warning'].agg(lambda s: 'WARN' if s.astype(str).str.upper().str.startswith('WARN').any() else 'PASS')
    meta['qc_warning'] = qc.reindex(meta.index).values
    panels = '11 Target 96 panels' if study == 'CITN-10' else 'Olink Explore (Cardiometabolic, Inflammation, Oncology, Neurology)'
    return dict(mat=mat, meta=meta, markers=mk, kind='per-sample',
                notes=f'{study}: NPX from {panels}; the IO panel was NOT used. Controls removed. Response (R/NR, PR/SD/PD) and '
                      'disease (MF/SS) joined by subject. Assays measured on several panels are kept as separate rows.')


def s086a():
    return _s086('olink_citn_10.csv', 'CITN-10')


def s086b():
    return _s086('olink_citn_13.csv', 'CITN-13')


def _parker(folder, files, clin, note):
    d = pd.concat([pd.read_csv(P(folder, f)) for f in files])
    d['sample_id'] = d['manuscript.id'] + '_' + d['timepoint.id']
    mat = d.pivot_table(index='protein.name', columns='sample_id', values='npx', aggfunc='mean')
    sm = d.drop_duplicates('sample_id').set_index('sample_id')[['manuscript.id', 'cohort', 'treatment', 'timepoint.id']]
    cl = pd.read_csv(P(folder, clin)).set_index('manuscript.id')
    meta = sm.join(cl.drop(columns=['cohort', 'treatment'], errors='ignore'), on='manuscript.id')
    meta.insert(0, 'sample_id', meta.index)
    meta.columns = snake(meta.columns)
    meta = meta.reindex(mat.columns)
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology + Immune Response'), kind='per-sample', notes=note)


def s093():
    return _parker('row093_OReilly2026_REVOLUTION',
                   ['translational-data_REV_cohortA_olink_public.csv', 'translational-data_REV_cohortB_olink_public.csv'],
                   'clinical-data_REV_clinical_data_public.csv',
                   'NPX; cohorts A and B combined (cohort column). 172 proteins (Target 96 Immuno-Oncology + Immune Response panels). '
                   'Outcomes (BOR, PFS, OS, 1-yr OS status) joined by patient.')


def s105():
    return _parker('row105_Galsky2025_PORTER', ['translational-data_PORTER_olink_public.csv'],
                   'clinical-data_PORTER_clinical_data_public.csv',
                   'NPX; cohorts A-C. 172 proteins (Target 96 Immuno-Oncology + Immune Response panels). Outcomes (BOR, rPFS, OS, composite '
                   'responder, PSA/CTC response) joined by patient.')


def s096():
    rows = []
    x = pd.read_excel(P('row096_Fedak2026_MG', '41591_2025_4170_MOESM11_ESM.xlsx'), header=None)
    j = 0
    while j < x.shape[1]:
        name = x.iat[0, j]
        if isinstance(name, str) and name.strip():
            k = j + 1
            while k < x.shape[1] and isinstance(x.iat[0, k], str) and x.iat[0, k].strip():
                for v in pd.to_numeric(x.iloc[1:, k], errors='coerce').dropna():
                    rows.append((name.strip(), x.iat[0, k].replace('\xa0', ' ').strip(), 'Extended Data Fig. 3', float(v)))
                k += 1
            j = k
        else:
            j += 1
    s = pd.read_excel(P('row096_Fedak2026_MG', '41591_2025_4170_MOESM6_ESM.xlsx'), sheet_name='4N', header=None)
    for jj in range(1, s.shape[1]):
        for v in pd.to_numeric(s.iloc[1:, jj], errors='coerce').dropna():
            rows.append((str(s.iat[0, jj]).strip(), 'D1 vs SCR', 'Fig. 4N', float(v)))
    L = pd.DataFrame(rows, columns=['marker', 'group', 'source', 'delta_npx'])
    mat = L.groupby(['marker', 'group'])['delta_npx'].median().unstack()
    mk = simple_markers(mat.index, 'Target 96 Inflammation (selected)')
    mk['n_values'] = L.groupby('marker').size().reindex(mat.index).values
    meta = pd.DataFrame({'sample_id': mat.columns, 'statistic': 'median delta-NPX vs screening across patients'}).set_index('sample_id', drop=False)
    return dict(mat=mat, meta=meta, markers=mk, kind='summary', long=L,
                notes='SUMMARY ONLY: figure source data give delta-NPX values without patient IDs, so values are summarised as medians per '
                      'group/visit (columns). All individual values are in long_values.csv.')


def s103():
    d = pd.read_excel(P('row103_Yu2025_LGL', '41392_2025_2164_MOESM1_ESM.xlsx')).set_index('Protein')
    mat = d.apply(pd.to_numeric, errors='coerce')
    meta = pd.DataFrame(index=mat.columns)
    meta['sample_id'] = meta.index
    meta['group'] = np.where(meta.index.str.startswith('Health'), 'Healthy donor', 'LGL leukemia patient')
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Inflammation', uniprot=list(mat.index)), kind='per-sample',
                notes='Values are on a LINEAR scale as published (not log2 NPX); a log2 copy is also written. Markers are given as UniProt '
                      'accessions. Inflammation panel (not IO). Response data are in the paper, not in this table.')


def s062():
    f = P('row062_Coffey2024_Myeloma', 'jitc-2023-008110supp002_data_supplement.xlsx')
    d = pd.read_excel(f, sheet_name='Supplemental table 3', header=2)
    d = d.dropna(subset=['Patient ID'])
    d['sample_id'] = ['P%02d_%s' % (int(p), str(t).replace(' ', '')) for p, t in zip(d['Patient ID'], d['Time Point'])]
    x = d.drop(columns=['Patient ID', 'Time Point']).set_index('sample_id')
    mat = x.T.apply(pd.to_numeric, errors='coerce')
    t1 = pd.read_excel(f, sheet_name='Supplemental table 1', header=3).dropna(subset=['Patient ID'])
    t1.columns = ['patient', 'stage', 'cytogenetic_risk', 'fish', 'm_protein_immunophenotype', 'imwg_screening', 'imwg_day90',
                  'imwg_day365', 'imwg_best']
    t1['patient'] = t1['patient'].astype(int)
    meta = d[['sample_id', 'Patient ID', 'Time Point']].copy()
    meta.columns = ['sample_id', 'patient', 'timepoint']
    meta['patient'] = meta['patient'].astype(int)
    meta = meta.merge(t1, on='patient', how='left')
    col = {'Screening': 'imwg_screening', 'Day 90': 'imwg_day90', 'Day 365': 'imwg_day365'}
    meta['imwg_response_at_timepoint'] = [r[col[t]] if t in col else np.nan for t, (_, r) in zip(meta['timepoint'], meta.iterrows())]
    meta['patient_id'] = ['P%02d' % p for p in meta['patient']]
    meta = meta.drop(columns=['patient']).set_index('sample_id', drop=False).reindex(mat.columns)
    return dict(mat=mat, meta=meta, markers=simple_markers(mat.index, 'Target 96 Immuno-Oncology'), kind='per-sample',
                notes='Supplemental Table 3 ("unfiltered Olink normalized expression values"): NPX, Target 96 Immuno-Oncology, BONE MARROW '
                      'plasma (per the paper), 15 multiple myeloma patients x Screening / Day 90 / Day 365 after autologous stem cell '
                      'transplant with elotuzumab + lenalidomide. IMWG response per timepoint and best response from Supplemental '
                      'Table 1 (with ISS/DS stage, cytogenetic risk, FISH, M-protein). Per-patient PFS is not released (Supplemental '
                      'Table 5 has only Cox hazard ratios). 4 values are blank in the source and stay NaN.')


STUDIES = {
    'row010_Meng2026_GastricCa': s010, 'row012_Johnson2026_ccRCC': s012, 'row024_Shan2026_OvarianCAR-T': s024,
    'row062_Coffey2024_Myeloma': s062,
    'row027_Tikkanen2025_SolidTumors': s027, 'row032_Gao2025_NSCLC': s032, 'row033_Xiao2025_TNBC': s033,
    'row040_Plummer2025_TIMEPOINT': s040a, 'row040_Plummer2025_OUTREACH': s040b, 'row044_Braun2025_RCCvaccine': s044,
    'row052_Tsimberidou2024_AMADEUS': s052, 'row070_Nunez2023_irAE': s070, 'row072_Padron2022_PRINCE': s072,
    'row073_Liu2022_TNBC': s073, 'row078_Barroso-Sousa2021_TNBC': s078, 'row086_Glass2024_CTCL_CITN10': s086a,
    'row086_Glass2024_CTCL_CITN13': s086b, 'row093_OReilly2026_REVOLUTION': s093, 'row096_Fedak2026_MG': s096,
    'row103_Yu2025_LGL': s103, 'row105_Galsky2025_PORTER': s105}
