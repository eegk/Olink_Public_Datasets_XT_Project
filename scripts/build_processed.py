"""Write cleaned per-study outputs and the cross-study marker-overlap summary.

Run from anywhere:  python3 build_processed.py
Outputs go to Olink_IO_data/processed/<dataset>/ and Olink_IO_data/processed/_overlap/.
"""
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import load_studies as L  # noqa: E402

OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

# UniProt accessions for names that carry no accession in their source file
# (UniProt mnemonics used by the Parker Institute files, plus common aliases). Resolved via rest.uniprot.org.
MANUAL = {
    'AGRG1': 'Q9Y653', 'ANGP1': 'Q15389', 'ANGP2': 'O15123', 'ARGI1': 'P05089', 'BGAL': 'P16278', 'BT3A2': 'P78410',
    'CAH9': 'Q16790', 'CLC4A': 'Q9UMR7', 'CLC4C': 'Q8WTT0', 'CLC4D': 'Q8WXI8', 'CLC4G': 'Q6UXB4', 'CLC6A': 'Q6EIG7',
    'CLC7A': 'Q9BXN2', 'CNTP2': 'Q9UHC6', 'CXAR': 'P78310', 'CXL10': 'P02778', 'CXL11': 'O14625', 'CXL13': 'O43927',
    'DCBD2': 'Q96PD2', 'DHI1': 'P28845', 'GALT3': 'Q14435', 'GRAA': 'P12544', 'GRAB': 'P10144', 'GRAH': 'P20718',
    'GROA': 'P09341', 'HEXI1': 'O94992', 'I12R1': 'P42701', 'ICA69': 'Q05084', 'IF4G1': 'Q04637', 'IF5A1': 'P63241',
    'EIF5A1': 'P63241', 'IMA5': 'P52294', 'INLR1': 'Q8IU57', 'ITA11': 'Q9UKX5', 'ITA6': 'P23229', 'ITB6': 'P18564',
    'K1C19': 'P08727', 'KI3L1': 'P43629', 'KPCT': 'Q04759', 'LEG1': 'P09382', 'LEG9': 'O00182', 'LIRB4': 'Q8NHJ6',
    'MERL': 'P35240', 'NCTR1': 'O76036', 'NEB2': 'Q96SB3', 'NFAC3': 'Q12968', 'PD1L1': 'Q9NZQ7', 'PD1L2': 'Q9BQ51',
    'PGS2': 'P07585', 'PLXA4': 'Q9HCM2', 'RO52': 'P19474', 'SDF1': 'P48061', 'SH21A': 'O60880', 'SPY2': 'O43597',
    'TNF10': 'P50591', 'TNF12': 'O43508', 'TNF14': 'O43557', 'TNFA': 'P01375', 'TNFL6': 'P48023', 'TNR12': 'Q9NP84',
    'TNR21': 'O75509', 'TNR4': 'P43489', 'TNR5': 'P25942', 'TNR9': 'Q07011', 'TRYB1': 'Q15661', 'VGFR2': 'P35968',
    'X3CL1': 'P78423', 'ZBT16': 'Q05516', 'BCAP': 'Q6ZUJ8', 'ICOSL': 'O75144', 'PLGF': 'P49763', 'PIGF': 'P49763',
    'RIGI': 'O95786', 'TNFRSF5': 'P25942', '4-1BB': 'Q07011', 'RANKL': 'O14788', 'X4EBP1': 'Q13541', 't_4EBP1': 'Q13541',
    'MICA': 'Q29980,Q29983', 'MIC_AvB': 'Q29980,Q29983', 'IL12A': 'P29459,P29460'}

# Preferred display names (Olink Target 96 assay names) where sources disagree
PREFERRED = {'P10145': 'IL8', 'Q07011': 'TNFRSF9', 'Q02763': 'TIE2', 'P80098': 'MCP-3', 'P29965': 'CD40-L', 'P01583': 'IL-1 alpha',
             'Q16790': 'CAIX', 'P01579': 'IFN-gamma', 'P13500': 'MCP-1', 'P80075': 'MCP-2', 'Q99616': 'MCP-4',
             'Q29980,Q29983': 'MIC-A/B', 'P29459,P29460': 'IL12', 'P01137': 'LAP TGF-beta-1', 'P01127': 'PDGF subunit B',
             'Q9NZQ7': 'PD-L1', 'Q9BQ51': 'PD-L2', 'P35968': 'VEGFR-2', 'P09382': 'Gal-1', 'O00182': 'Gal-9',
             'P09601': 'HO-1', 'P50591': 'TRAIL', 'O43508': 'TWEAK', 'Q9Y653': 'ADGRG1', 'P28845': 'HSD11B1', 'O14788': 'TRANCE',
             'P01375': 'TNF', 'Q14790': 'CASP8', 'P09603': 'CSF1', 'P49763': 'PGF', 'Q13541': '4E-BP1', 'P05089': 'ARG1', 'Q15389': 'ANGPT1', 'O15123': 'ANGPT2',
             'P25942': 'CD40', 'O75144': 'ICOSLG', 'P48023': 'FASLG', 'P78423': 'CX3CL1', 'P48061': 'CXCL12',
             'P09341': 'CXCL1', 'P02778': 'CXCL10', 'O14625': 'CXCL11', 'O43927': 'CXCL13', 'P12544': 'GZMA',
             'P10144': 'GZMB', 'P20718': 'GZMH', 'P43489': 'TNFRSF4', 'Q9NP84': 'TNFRSF12A', 'O75509': 'TNFRSF21',
             'O43557': 'TNFSF14', 'P42701': 'IL12RB1', 'P43629': 'KIR3DL1', 'O76036': 'NCR1', 'P07585': 'DCN'}


def norm(s):
    s = str(s).upper().replace('ALPHA', 'A').replace('BETA', 'B').replace('GAMMA', 'G')
    return re.sub(r'[^A-Z0-9]', '', s)


# Older Olink annotations use superseded accessions for a few proteins; map them to the current entry
ACC_ALIAS = {'Q8WWJ7': 'P30203',    # CD6 (older Inflammation-panel annotation)
             'Q8NF90': 'P12034',    # FGF5
             'P0CG05': 'P0DOY2',    # IGLC2
             'NTproBNP': 'P16860'}  # NT-proBNP has no accession in Olink files; NPPB gene product


def canon_acc(u):
    if not isinstance(u, str) or u in ('', 'nan'):
        return None
    parts = [re.sub(r'-\d+$', '', p) for p in re.split(r'[,_;\s]+', u.strip()) if p]
    parts = [ACC_ALIAS.get(p, p) for p in parts]
    return ','.join(sorted(set(parts))) if parts else None


def main():
    R = {}
    for k, f in L.STUDIES.items():
        R[k] = f()
        print('loaded', k, R[k]['mat'].shape)

    # 1) name -> accession reference built from every file that carries accessions, plus the Inflammation panel list
    ref = {}
    for r in R.values():
        for nm, u in zip(r['markers']['original_name'], r['markers']['uniprot']):
            a = canon_acc(u)
            if a:
                ref.setdefault(norm(nm), Counter())[a] += 1
    pl = pd.read_csv(os.path.join(L.BASE, 'row040_Plummer2025_TIMEPOINT', 'OlinkTarget96Inflammation-2025-01-27.csv'))
    for v in pl.iloc[:, 0]:
        a = str(v).split(';')
        ref.setdefault(norm(a[1]), Counter())[canon_acc(a[0])] += 5
    for nm, a in MANUAL.items():
        ref[norm(nm)] = Counter({canon_acc(a): 100})
    best = {n: c.most_common(1)[0][0] for n, c in ref.items()}

    # 2) assign accession to every marker row
    unresolved = []
    for k, r in R.items():
        mk = r['markers']
        acc = []
        for nm, u in zip(mk['original_name'], mk['uniprot']):
            a = canon_acc(u) or best.get(norm(nm))
            if a is None:
                unresolved.append((k, nm))
            acc.append(a)
        mk['uniprot_harmonized'] = acc

    # 3) display name per accession: preferred, else most common original name
    names = {}
    for r in R.values():
        for nm, a in zip(r['markers']['original_name'], r['markers']['uniprot_harmonized']):
            if a:
                names.setdefault(a, Counter())[str(nm)] += 1
    disp = {a: PREFERRED.get(a, c.most_common(1)[0][0]) for a, c in names.items()}

    # 4) write per-study outputs
    inventory = []
    for k, r in R.items():
        d = os.path.join(OUT, k)
        os.makedirs(d, exist_ok=True)
        mk = r['markers'].copy()
        mk['marker'] = [disp.get(a, nm) if a else nm for a, nm in zip(mk['uniprot_harmonized'], mk['original_name'])]
        # unique row labels: harmonized name, suffixed with panel (then source key) if repeated within the study
        lab = mk['marker'].copy()
        orig = mk['original_name'].astype(str)
        dup = lab.duplicated(keep=False) & (orig != lab)
        lab[dup] = lab[dup] + '|' + orig[dup]
        dup = lab.duplicated(keep=False)
        pan = mk['panel'].astype(str).str.replace(r'Olink (Target 96 )?', '', regex=True).str.replace(r'\(.*\)', '', regex=True).str.strip()
        lab[dup] = lab[dup] + '|' + pan[dup]
        dup = lab.duplicated(keep=False)
        lab[dup] = lab[dup] + '|' + pd.Series(mk.index, index=mk.index)[dup].astype(str)
        mat = r['mat'].copy()
        mat.index = lab.loc[mat.index].values
        mat.index.name = 'marker'
        mk.insert(0, 'marker_row_label', lab.values)
        mk = mk.drop(columns=['marker'])
        mk.insert(1, 'marker', [disp.get(a, nm) if a else nm for a, nm in zip(mk['uniprot_harmonized'], mk['original_name'])])
        meta = r['meta'].copy()
        meta.index.name = None
        mat.to_csv(os.path.join(d, 'expression_matrix.csv'))
        if k.startswith('row103'):
            np.log2(mat.clip(lower=1e-9)).to_csv(os.path.join(d, 'expression_matrix_log2.csv'))
        if 'long' in r:
            r['long'].to_csv(os.path.join(d, 'long_values.csv'), index=False)
        meta.to_csv(os.path.join(d, 'sample_metadata.csv'), index=False)
        mk.to_csv(os.path.join(d, 'marker_annotation.csv'), index=False)
        with open(os.path.join(d, 'README.txt'), 'w') as fh:
            fh.write(f'{k}\nType: {r["kind"]}\nMatrix: {mat.shape[0]} markers (rows) x {mat.shape[1]} '
                     f'{"samples" if r["kind"] == "per-sample" else "summary columns"} (columns)\n\n{r["notes"]}\n\n'
                     'Files: expression_matrix.csv (rows = harmonized marker names; "|panel" suffix only when an assay appears on '
                     'several panels), sample_metadata.csv (one row per matrix column, same order), marker_annotation.csv '
                     '(original assay name, UniProt, panel, OlinkID).\n')
        oc = [c for c in meta.columns if re.search(r'resp|orr|bor|best_overall|pfs|os_|^os$|os_months|surviv|dead|death|relapse|tox|irae|disease_control|composite|benefit', c, re.I)]
        inventory.append(dict(dataset=k, type=r['kind'], n_markers=mat.shape[0], n_unique_proteins=mk['uniprot_harmonized'].nunique(),
                              n_columns=mat.shape[1], n_patients=(meta[[c for c in meta.columns if c in ('patient_id', 'subject_id', 'subject', 'manuscript_id', 'deidentified_id')][0]].nunique()
                                                                   if any(c in meta.columns for c in ('patient_id', 'subject_id', 'subject', 'manuscript_id', 'deidentified_id')) else np.nan),
                              outcome_columns='; '.join(oc), metadata_columns=meta.shape[1], notes=r['notes']))

    inv = pd.DataFrame(inventory)
    inv.to_csv(os.path.join(OUT, 'dataset_inventory.csv'), index=False)

    # 5) overlap
    od = os.path.join(OUT, '_overlap')
    os.makedirs(od, exist_ok=True)
    pres = {}
    for k, r in R.items():
        pres[k] = set(a for a in r['markers']['uniprot_harmonized'] if a)
    allacc = sorted(set().union(*pres.values()))
    P = pd.DataFrame({k: [int(a in s) for a in allacc] for k, s in pres.items()}, index=allacc)
    persample = [k for k in R if R[k]['kind'] == 'per-sample']
    P.insert(0, 'marker', [disp.get(a, a) for a in P.index])
    P.insert(1, 'n_datasets_all', P[list(R)].sum(axis=1))
    P.insert(2, 'n_datasets_per_sample', P[persample].sum(axis=1))
    P.index.name = 'uniprot'
    P = P.sort_values(['n_datasets_all', 'marker'], ascending=[False, True])
    P.to_csv(os.path.join(od, 'marker_presence_by_dataset.csv'))
    ks = list(R)
    pw = pd.DataFrame([[len(pres[a] & pres[b]) for b in ks] for a in ks], index=ks, columns=ks)
    pw.to_csv(os.path.join(od, 'pairwise_shared_markers.csv'))
    pd.DataFrame(unresolved, columns=['dataset', 'original_name']).to_csv(os.path.join(od, 'unresolved_marker_names.csv'), index=False)
    print('unresolved', unresolved)
    print(inv[['dataset', 'type', 'n_markers', 'n_unique_proteins', 'n_columns', 'n_patients', 'outcome_columns']].to_string())


if __name__ == '__main__':
    main()
