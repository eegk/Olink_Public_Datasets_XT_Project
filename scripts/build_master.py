"""Combine the 17 per-sample datasets in processed/ into one master object. NO normalization: every expression
value is copied unchanged from processed/<dataset>/expression_matrix.csv (which equals the raw source files).

Outputs (processed/_master/):
  master_long.csv.gz         every value of every assay (lossless): sample_uid, dataset, marker, uniprot, assay row, panel, value
  master_matrix.csv.gz       proteins (rows, UniProt-harmonized) x samples (columns); one assay per protein per dataset
  master_features.csv        one row per protein: name, UniProt, number of datasets/samples measured, IO-core flag, panels
  master_metadata.csv        one row per sample: harmonized core columns (study, disease, treatment, timepoint, baseline,
                             response, BOR, PFS/OS, irAE, sample type, panel, value scale ...)
  master_metadata_full.csv   core columns + every original metadata column (prefixed orig__)
  master_object.h5ad         AnnData: X = samples x proteins (raw values, NaN = not measured), obs = master metadata,
                             var = features, layers['measured'] mask
"""
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PROC = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(PROC, '_master')
os.makedirs(OUT, exist_ok=True)

SUMMARY_ONLY = ['row044_Braun2025_RCCvaccine', 'row078_Barroso-Sousa2021_TNBC', 'row096_Fedak2026_MG']
PANEL_PRIORITY = ['Immuno-Oncology', 'Inflammation', 'Oncology', 'Immune Response', 'selected', 'measurement 1']

# ------------------------------------------------------------------ study-level facts (from the papers / deposits)
STUDY = {
    'row010_Meng2026_GastricCa': dict(study='Meng 2026, Mol Biomed', doi='10.1186/s43556-026-00510-8', disease='Gastric/GEJ adenocarcinoma',
                                      treatment='Anti-PD-1 + chemotherapy', sample_type='Plasma', panels='Target 96 Immuno-Oncology', value_scale='NPX (log2)'),
    'row012_Johnson2026_ccRCC': dict(study='Johnson 2026, Nat Commun', doi='10.1038/s41467-026-73336-1', disease='Clear cell renal cell carcinoma',
                                     treatment='Pembrolizumab, then pembrolizumab + high-dose IL-2', sample_type='Plasma',
                                     panels='Target 96 Immuno-Oncology; Oncology II; Cardiometabolic', value_scale='NPX (log2)'),
    'row024_Shan2026_OvarianCAR-T': dict(study='Shan 2026, Cancer', doi='10.1002/cncr.70279', disease='Ovarian cancer (platinum-resistant)',
                                         treatment='Mesothelin TCR-like CAR-T', sample_type='Plasma', panels='Target 96 Oncology II', value_scale='NPX (log2)'),
    'row062_Coffey2024_Myeloma': dict(study='Coffey 2024, J Immunother Cancer', doi='10.1136/jitc-2023-008110', disease='Multiple myeloma',
                                      treatment='Elotuzumab + autologous stem cell transplant + lenalidomide', sample_type='Bone marrow plasma',
                                      panels='Target 96 Immuno-Oncology', value_scale='NPX (log2)'),
    'row027_Tikkanen2025_SolidTumors': dict(study='Tikkanen 2025, Sci Rep', doi='10.1038/s41598-025-31229-1', disease='Advanced solid tumours',
                                            treatment='Anti-PD-(L)1', sample_type='Plasma', panels='Target 96 Immuno-Oncology', value_scale='NPX (log2)'),
    'row032_Gao2025_NSCLC': dict(study='Gao 2025, Mol Biomed', doi='10.1186/s43556-025-00291-6', disease='NSCLC',
                                 treatment='Anti-PD-1/PD-L1 (+/- chemotherapy)', sample_type='Plasma', panels='Target 96 Immuno-Oncology',
                                 value_scale='Centred NPX as deposited (not raw NPX)'),
    'row033_Xiao2025_TNBC': dict(study='Xiao 2025, Cancer Biol Med', doi='10.20892/j.issn.2095-3941.2025.0038', disease='Triple-negative breast cancer',
                                 treatment='Immunotherapy or chemotherapy (see treatment_detail)', sample_type='Plasma', panels='Target 96 Immuno-Oncology', value_scale='NPX (log2)'),
    'row040_Plummer2025_TIMEPOINT': dict(study='Plummer 2025, Cell Rep Med (TIMEPOINT)', doi='10.1016/j.xcrm.2025.102041', disease='Advanced solid tumours',
                                         treatment='MTL-CEBPA + pembrolizumab', sample_type='Plasma', panels='Target 96 Inflammation', value_scale='NPX (log2)'),
    'row040_Plummer2025_OUTREACH': dict(study='Plummer 2025, Cell Rep Med (OUTREACH)', doi='10.1016/j.xcrm.2025.102041', disease='Advanced solid tumours',
                                        treatment='MTL-CEBPA +/- pembrolizumab', sample_type='Plasma', panels='Target 96 Inflammation', value_scale='NPX (log2)'),
    'row052_Tsimberidou2024_AMADEUS': dict(study='Tsimberidou 2024, J Exp Med (AMADEUS)', doi='10.1084/jem.20240152', disease='Advanced solid tumours',
                                           treatment='Nivolumab +/- ipilimumab', sample_type='Plasma', panels='Target 96 Immuno-Oncology; Immune Response', value_scale='NPX (log2)'),
    'row070_Nunez2023_irAE': dict(study='Nunez 2023, Med', doi='10.1016/j.medj.2022.12.007', disease='Melanoma or NSCLC (not given per patient)',
                                  treatment='Immune checkpoint inhibitors', sample_type='Serum', panels='Target 96 Inflammation', value_scale='NPX (log2)'),
    'row072_Padron2022_PRINCE': dict(study='Padron 2022, Nat Med (PRINCE)', doi='10.1038/s41591-022-01829-9', disease='Metastatic pancreatic adenocarcinoma',
                                     treatment='Gemcitabine/nab-paclitaxel + nivolumab and/or sotigalimab', sample_type='Plasma',
                                     panels='Olink (30 selected proteins released)', value_scale='NPX (log2)'),
    'row073_Liu2022_TNBC': dict(study='Liu 2022, Nat Commun', doi='10.1038/s41467-022-30569-0', disease='Triple-negative breast cancer',
                                treatment='Camrelizumab + apatinib + eribulin', sample_type='Plasma', panels='Target 96 Immuno-Oncology', value_scale='NPX (log2)'),
    'row086_Glass2024_CTCL_CITN10': dict(study='Glass 2024, Cell Rep Med (CITN-10)', doi='10.1016/j.xcrm.2024.101527', disease='Cutaneous T cell lymphoma',
                                         treatment='Pembrolizumab', sample_type='Serum', panels='11 Target 96 panels (no Immuno-Oncology)', value_scale='NPX (log2)'),
    'row086_Glass2024_CTCL_CITN13': dict(study='Glass 2024, Cell Rep Med (CITN-13)', doi='10.1016/j.xcrm.2024.101527', disease='Cutaneous T cell lymphoma',
                                         treatment='Pembrolizumab + interferon-gamma', sample_type='Serum',
                                         panels='Olink Explore (Cardiometabolic, Inflammation, Oncology, Neurology)', value_scale='NPX (log2)'),
    'row093_OReilly2026_REVOLUTION': dict(study="O'Reilly 2026, JITC (REVOLUTION)", doi='10.1136/jitc-2025-012864', disease='Metastatic pancreatic adenocarcinoma',
                                          treatment='Gemcitabine/nab-paclitaxel + ipilimumab + nivolumab or hydroxychloroquine', sample_type='Serum',
                                          panels='Target 96 Immuno-Oncology; Immune Response', value_scale='NPX (log2)'),
    'row103_Yu2025_LGL': dict(study='Yu 2025, Signal Transduct Target Ther', doi='10.1038/s41392-025-02164-4', disease='Large granular lymphocytic leukemia',
                              treatment='Thalidomide-based regimen', sample_type='Serum', panels='Target 96 Inflammation',
                              value_scale='LINEAR scale as published (not log2 NPX)'),
    'row105_Galsky2025_PORTER': dict(study='Galsky 2025, Clin Cancer Res (PORTER)', doi='10.1158/1078-0432.CCR-24-3693', disease='Metastatic castration-resistant prostate cancer',
                                     treatment='Combination immunotherapy (see treatment_detail)', sample_type='Serum',
                                     panels='Target 96 Immuno-Oncology; Immune Response', value_scale='NPX (log2)'),
}

BOR_MAP = {'COMPLETE RESPONSE': 'CR', 'CR': 'CR', 'PARTIAL RESPONSE': 'PR', 'PR': 'PR', 'CONFIRMED PARTIAL RESPONSE': 'PR',
           'UNCONFIRMED PARTIAL RESPONSE': 'uPR', 'STABLE DISEASE': 'SD', 'SD': 'SD', 'PROGRESSIVE DISEASE': 'PD', 'PD': 'PD',
           'NOT EVALUABLE': 'NE', 'NE': 'NE'}


def bor(x):
    if pd.isna(x):
        return np.nan
    return BOR_MAP.get(str(x).strip().upper(), np.nan)


def orr_from_bor(b):
    if pd.isna(b) or b == 'NE':
        return np.nan
    return 'R' if b in ('CR', 'PR') else 'NR'


def yn(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().upper()
    if s in ('1', '1.0', 'Y', 'YES', 'TRUE'):
        return 1
    if s in ('0', '0.0', 'N', 'NO', 'FALSE'):
        return 0
    return np.nan


def harmonize(ds, m):
    """Return DataFrame of harmonized core columns for dataset ds (index aligned with m)."""
    h = pd.DataFrame(index=m.index)
    g = lambda c: m[c] if c in m.columns else pd.Series(np.nan, index=m.index)  # noqa: E731
    h['patient_id_original'] = np.nan
    h['timepoint_original'] = np.nan
    h['timepoint_std'] = np.nan
    h['disease_detail'] = np.nan
    h['treatment_detail'] = np.nan
    h['group'] = 'patient'
    for c in ['response_std', 'response_definition', 'bor_std', 'orr_recist', 'pfs_months', 'pfs_event', 'pfs_type', 'os_months',
              'os_event', 'irae', 'age', 'sex']:
        h[c] = np.nan

    if ds == 'row010_Meng2026_GastricCa':
        h['patient_id_original'] = g('patient_id')
        h['timepoint_original'] = g('timepoint_label')
        h['timepoint_std'] = np.where(g('timepoint_index') == 1, 'baseline', 'on-treatment')
        h['response_std'] = g('response')
        h['response_definition'] = 'Authors\' responder grouping (R/NR)'
    elif ds == 'row012_Johnson2026_ccRCC':
        h['patient_id_original'] = g('patient_id')
        h['timepoint_original'] = g('visit_code')
        h['timepoint_std'] = np.where(g('visit_code') == 'B1W1', 'baseline', 'on-treatment')
    elif ds == 'row062_Coffey2024_Myeloma':
        h['patient_id_original'] = g('patient_id')
        h['timepoint_original'] = g('timepoint')
        h['timepoint_std'] = np.where(g('timepoint') == 'Screening', 'baseline', 'on/post-treatment')
        h['disease_detail'] = g('stage').astype(str) + '; cytogenetic risk: ' + g('cytogenetic_risk').astype(str)
        h['bor_std'] = 'IMWG ' + g('imwg_best').astype(str)
        h['response_std'] = g('imwg_best').map(lambda x: np.nan if pd.isna(x) else ('NR' if str(x).strip() in ('PD', 'SD', 'MR') else 'R'))
        h['response_definition'] = 'IMWG best response >= PR (sCR/CR/VGPR/PR) = R; myeloma criteria, not RECIST'
    elif ds == 'row024_Shan2026_OvarianCAR-T':
        h['timepoint_std'] = 'unknown'
    elif ds == 'row027_Tikkanen2025_SolidTumors':
        h['patient_id_original'] = g('sample_id')
        h['timepoint_original'] = 'pre-treatment'
        h['timepoint_std'] = 'baseline'
        h['disease_detail'] = g('tumor_type')
        h['response_std'] = g('orr').map({1: 'R', 0: 'NR'})
        h['response_definition'] = 'Objective response (RECIST), as coded ORR 1/0'
        h['sex'] = g('sex')
    elif ds == 'row032_Gao2025_NSCLC':
        h['patient_id_original'] = g('sample_id')
        coh2 = g('cohort_inferred').astype(str).str.startswith('Cohort 2')
        h['timepoint_original'] = np.where(coh2, 'pre-treatment (validation cohort)', 'not stated')
        h['timepoint_std'] = np.where(coh2, 'baseline', 'unknown')
        h['bor_std'] = g('response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = h['orr_recist']
        h['response_definition'] = 'Objective response from BOR (PR vs SD/PD)'
        h['pfs_months'], h['pfs_event'], h['pfs_type'] = g('pfs'), g('relapse').map(yn), 'PFS'
        h['os_months'], h['os_event'] = g('os'), g('dead').map(yn)
        h['age'] = g('age').map({0.0: 'age group 0 (coded)', 1.0: 'age group 1 (coded)'})
        h['sex'] = g('gender').map({1.0: 'coded 1', 2.0: 'coded 2'})
    elif ds == 'row033_Xiao2025_TNBC':
        h['patient_id_original'] = g('patient_id')
        h['timepoint_original'] = g('timepoint')
        h['timepoint_std'] = np.where(g('timepoint') == 'Pre', 'baseline', 'on/post-treatment')
        h['treatment_detail'] = g('group').map({'GroupA': 'Neoadjuvant', 'GroupB': 'Adjuvant', 'GroupC': 'Advanced'}) + ' - ' + g('treatment').astype(str)
        r = g('response').astype(str)
        h['bor_std'] = r.map(lambda x: x if x in ('CR', 'PR', 'SD', 'PD') else np.nan)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = r.map({'Response': 'R', 'nonResponse': 'NR', 'CR': 'R', 'PR': 'R', 'SD': 'NR', 'PD': 'NR'})
        h['response_definition'] = np.where(g('group') == 'GroupC', 'RECIST objective response (advanced)',
                                            np.where(g('group') == 'GroupA', 'Authors\' neoadjuvant response (Response/nonResponse)', 'not assessed'))
        h['age'] = g('age')
    elif ds in ('row040_Plummer2025_TIMEPOINT', 'row040_Plummer2025_OUTREACH'):
        h['patient_id_original'] = g('sample_id')
        h['timepoint_original'] = 'C1D1'
        h['timepoint_std'] = 'baseline'
        if ds.endswith('TIMEPOINT'):
            h['disease_detail'] = g('tumour_type')
            h['bor_std'] = g('response').map(bor)
        else:
            h['bor_std'] = g('best_overall_response').map(bor)
            h['treatment_detail'] = g('treatment')
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = h['orr_recist']
        h['response_definition'] = 'Objective response from BOR' + (' (TIMEPOINT reports only PD/SD)' if ds.endswith('TIMEPOINT') else '')
    elif ds == 'row052_Tsimberidou2024_AMADEUS':
        h['patient_id_original'] = g('subject_id')
        h['timepoint_original'] = g('timepoint_id')
        h['timepoint_std'] = np.where(g('timepoint_id') == 'C1D1', 'baseline', 'on-treatment')
        h['disease_detail'] = g('tumor_type')
        h['treatment_detail'] = g('treatment') + ' (' + g('arm') + ' arm)'
        h['bor_std'] = g('best_overall_response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = g('responder_flag_crpr').map({'Y': 'R', 'N': 'NR'})
        h['response_definition'] = 'Responder flag CR/PR (study file)'
        h['pfs_months'], h['pfs_event'], h['pfs_type'] = g('pfs_months'), g('pfs_event_flag').map(yn), 'PFS'
        h['os_months'], h['os_event'] = g('os_months'), g('os_event_flag').map(yn)
        h['age'], h['sex'] = g('age'), g('sex')
    elif ds == 'row070_Nunez2023_irAE':
        # same person listed twice in the source (series 1/2, identical T1/T2 values): one patient id; series kept in orig__ columns
        h['patient_id_original'] = g('patient_id').astype(str)
        h['timepoint_original'] = g('timepoint_label')
        h['timepoint_std'] = np.where(g('timepoint') == 1, 'baseline', 'on-treatment')
        h['response_std'] = g('response')
        h['response_definition'] = 'Authors\' response group (Tox R / No Tox R / No Tox NR)'
        h['irae'] = g('irae_toxicity').map({'Yes': 1, 'No': 0})
    elif ds == 'row072_Padron2022_PRINCE':
        h['patient_id_original'] = g('deidentified_id')
        h['timepoint_original'] = g('timepoint_id')
        h['timepoint_std'] = np.where(g('timepoint_id') == 'C1D1', 'baseline', 'on-treatment')
        h['treatment_detail'] = g('arm_description')
        h['bor_std'] = g('best_overall_response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = h['orr_recist']
        h['response_definition'] = 'Objective response from BOR'
        h['pfs_months'], h['pfs_event'], h['pfs_type'] = g('clinical_observation_pfs') / 30.4375, g('clinical_observation_pfs_event').map(yn), 'PFS'
        h['os_months'], h['os_event'] = g('clinical_observation_os') / 30.4375, g('clinical_observation_os_event').map(yn)
        h['age'], h['sex'] = g('age'), g('sex')
    elif ds == 'row073_Liu2022_TNBC':
        h['patient_id_original'] = g('patient_id').astype(str)
        h['timepoint_original'] = g('timepoint')
        h['timepoint_std'] = np.where(g('timepoint') == 'pre-treatment', 'baseline', 'on/post-treatment')
    elif ds in ('row086_Glass2024_CTCL_CITN10', 'row086_Glass2024_CTCL_CITN13'):
        h['patient_id_original'] = g('subject')
        h['timepoint_original'] = g('timepoint')
        tp = g('timepoint').astype(str)
        h['timepoint_std'] = np.where(tp.isin(['Pre', 'BL']), 'baseline', np.where(tp == 'EOT', 'end of treatment', 'on-treatment'))
        h['disease_detail'] = g('status').map({'MF': 'Mycosis fungoides', 'SS': 'Sezary syndrome'})
        h['bor_std'] = g('granular_response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = g('response')
        h['response_definition'] = 'Global response R (CR/PR) vs NR (study file)'
    elif ds == 'row093_OReilly2026_REVOLUTION':
        h['patient_id_original'] = g('manuscript_id')
        h['timepoint_original'] = g('timepoint_id')
        h['timepoint_std'] = np.where(g('timepoint_id') == 'C1D1', 'baseline', 'on-treatment')
        h['treatment_detail'] = 'Cohort ' + g('cohort') + ': ' + g('treatment')
        h['bor_std'] = g('best_overall_response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = h['orr_recist']
        h['response_definition'] = 'Objective response from BOR (confirmed PR = R; unconfirmed PR coded uPR and counted NR)'
        h['pfs_months'], h['pfs_event'], h['pfs_type'] = g('pfs_months'), g('pfs_event_flag').map(yn), 'PFS'
        h['os_months'], h['os_event'] = g('os_months'), g('os_event_flag').map(yn)
        h['age'], h['sex'] = g('age'), g('sex')
    elif ds == 'row103_Yu2025_LGL':
        h['group'] = np.where(g('group') == 'Healthy donor', 'healthy control', 'patient')
        h['patient_id_original'] = g('sample_id')
        h['timepoint_std'] = np.where(g('group') == 'Healthy donor', 'healthy control', 'unknown')
    elif ds == 'row105_Galsky2025_PORTER':
        h['patient_id_original'] = g('manuscript_id')
        h['timepoint_original'] = g('timepoint_id')
        tp = g('timepoint_id').astype(str)
        # cohort C has a lead-in treatment phase: Lead-in.Day1 is its baseline and Cycle1.Day1 is already on treatment
        pre = tp.isin(['BL', 'Lead-in.Day1']) | ((tp == 'Cycle1.Day1') & (g('cohort') != 'C'))
        h['timepoint_std'] = np.where(pre, 'baseline', 'on-treatment')
        h['treatment_detail'] = 'Cohort ' + g('cohort') + ': ' + g('treatment')
        h['bor_std'] = g('best_overall_response').map(bor)
        h['orr_recist'] = h['bor_std'].map(orr_from_bor)
        h['response_std'] = g('composite_endpoint_responder').map({'Y': 'R', 'N': 'NR'})
        h['response_definition'] = 'Composite endpoint responder (PSA/CTC/RECIST, study definition)'
        h['pfs_months'], h['pfs_event'], h['pfs_type'] = g('rpfs_months'), g('rpfs_event_flag').map(yn), 'rPFS'
        h['os_months'], h['os_event'] = g('os_months'), g('os_event_flag').map(yn)
        h['age'], h['sex'] = g('age'), g('sex')
    h['pfs_type'] = np.where(h['pfs_months'].notna(), h['pfs_type'], np.nan)
    return h


BASELINE_PRIORITY = {'BL': 0, 'Lead-in.Day1': 1, 'Pre': 0}


def main():
    datasets = [d for d in sorted(os.listdir(PROC)) if d.startswith('row') and d not in SUMMARY_ONLY]
    long_parts, meta_parts, full_parts = [], [], []
    for ds in datasets:
        mat = pd.read_csv(os.path.join(PROC, ds, 'expression_matrix.csv'), index_col=0)
        meta = pd.read_csv(os.path.join(PROC, ds, 'sample_metadata.csv'))
        ann = pd.read_csv(os.path.join(PROC, ds, 'marker_annotation.csv')).set_index('marker_row_label')
        assert list(mat.columns) == meta['sample_id'].astype(str).tolist(), ds
        prefix = ds.split('_')[0] + ('_' + ds.split('_')[-1] if ds.split('_')[-1] in ('TIMEPOINT', 'OUTREACH', 'CITN10', 'CITN13') else '')
        uid = [f'{prefix}__{s}' for s in mat.columns]
        # long table
        L = mat.copy()
        L.columns = uid
        L = L.reset_index().melt(id_vars='marker', var_name='sample_uid', value_name='value').dropna(subset=['value'])
        L = L.rename(columns={'marker': 'assay_row'})
        L['dataset'] = ds
        L['uniprot'] = L['assay_row'].map(ann['uniprot_harmonized'])
        L['marker'] = L['assay_row'].map(ann['marker'])
        L['original_assay_name'] = L['assay_row'].map(ann['original_name'])
        L['panel'] = L['assay_row'].map(ann['panel'])
        long_parts.append(L[['sample_uid', 'dataset', 'marker', 'uniprot', 'assay_row', 'original_assay_name', 'panel', 'value']])
        # metadata
        info = STUDY[ds]
        core = pd.DataFrame({'sample_uid': uid, 'dataset': ds, 'publication_sheet_row': int(ds[3:6]), 'study': info['study'],
                             'doi': info['doi'], 'sample_id_original': meta['sample_id'].astype(str).values})
        h = harmonize(ds, meta).reset_index(drop=True)
        core['patient_uid'] = [f'{prefix}__{p}' if pd.notna(p) else np.nan for p in h['patient_id_original']]
        core['patient_id_original'] = h['patient_id_original'].values
        core['group'] = h['group'].values
        core['disease'] = info['disease']
        core['disease_detail'] = h['disease_detail'].values
        core['treatment'] = info['treatment']
        core['treatment_detail'] = h['treatment_detail'].values
        core['sample_type'] = info['sample_type']
        core['olink_panels'] = info['panels']
        core['value_scale'] = info['value_scale']
        core['timepoint_original'] = h['timepoint_original'].values
        core['timepoint_std'] = h['timepoint_std'].values
        for c in ['response_std', 'response_definition', 'bor_std', 'orr_recist', 'pfs_months', 'pfs_event', 'pfs_type', 'os_months',
                  'os_event', 'irae', 'age', 'sex']:
            core[c] = h[c].values
        core['qc_warning'] = meta['qc_warning'].values if 'qc_warning' in meta.columns else np.nan
        core['n_assays_measured'] = mat.notna().sum(axis=0).values
        meta_parts.append(core)
        orig = meta.copy()
        orig.columns = ['orig__' + c for c in orig.columns]
        full_parts.append(pd.concat([core.reset_index(drop=True), orig.reset_index(drop=True)], axis=1))
        print(f'{ds:34s} samples={len(uid):4d} values={len(L):7d}')

    LONG = pd.concat(long_parts, ignore_index=True)
    META = pd.concat(meta_parts, ignore_index=True)
    FULL = pd.concat(full_parts, ignore_index=True)

    # one primary baseline sample per patient (earliest pre-treatment sample)
    META['is_baseline'] = META['timepoint_std'].eq('baseline')
    pr = META['timepoint_original'].map(BASELINE_PRIORITY).fillna(2)
    cand = META[META['is_baseline'] & META['patient_uid'].notna()].assign(_p=pr)
    first = cand.sort_values(['patient_uid', '_p', 'sample_uid']).drop_duplicates('patient_uid')['sample_uid']
    META['is_primary_baseline'] = META['sample_uid'].isin(first)
    FULL = FULL.merge(META[['sample_uid', 'is_baseline', 'is_primary_baseline']], on='sample_uid')

    # protein-level wide matrix: one assay per protein per dataset (panel priority, then file order)
    def prio(p):
        p = str(p)
        for i, k in enumerate(PANEL_PRIORITY):
            if k.lower() in p.lower():
                return i
        return len(PANEL_PRIORITY)
    ar = LONG[['dataset', 'uniprot', 'assay_row', 'panel']].drop_duplicates()
    ar['_prio'] = ar['panel'].map(prio)
    ar['_order'] = ar.groupby('dataset').cumcount()
    primary = ar.sort_values(['dataset', 'uniprot', '_prio', '_order']).drop_duplicates(['dataset', 'uniprot'])
    LONG = LONG.merge(primary[['dataset', 'assay_row']].assign(used_in_master_matrix=True), on=['dataset', 'assay_row'], how='left')
    LONG['used_in_master_matrix'] = LONG['used_in_master_matrix'].fillna(False).astype(bool)
    W = LONG[LONG['used_in_master_matrix']].pivot(index='uniprot', columns='sample_uid', values='value')
    W = W.reindex(columns=META['sample_uid'])

    # features
    names = LONG.groupby('uniprot')['marker'].agg(lambda s: s.value_counts().index[0])
    io_sets = [set(LONG.loc[LONG.dataset == d, 'uniprot']) for d in ['row010_Meng2026_GastricCa', 'row012_Johnson2026_ccRCC', 'row062_Coffey2024_Myeloma',
               'row027_Tikkanen2025_SolidTumors', 'row032_Gao2025_NSCLC', 'row033_Xiao2025_TNBC', 'row052_Tsimberidou2024_AMADEUS',
               'row073_Liu2022_TNBC', 'row093_OReilly2026_REVOLUTION', 'row105_Galsky2025_PORTER']]
    io_core = set.intersection(*io_sets)
    F = pd.DataFrame(index=W.index)
    F['marker'] = names.reindex(F.index)
    F['n_datasets'] = LONG.groupby('uniprot')['dataset'].nunique().reindex(F.index)
    F['n_samples_measured'] = W.notna().sum(axis=1)
    F['in_IO_core_91'] = F.index.isin(io_core)
    F['datasets'] = LONG.groupby('uniprot')['dataset'].agg(lambda s: '; '.join(sorted(set(s)))).reindex(F.index)
    F['panels'] = LONG.groupby('uniprot')['panel'].agg(lambda s: '; '.join(sorted(set(map(str, s))))).reindex(F.index)
    F['n_alternative_assays_not_in_matrix'] = LONG[~LONG.used_in_master_matrix].groupby('uniprot')['assay_row'].nunique().reindex(F.index).fillna(0).astype(int)
    F = F.sort_values(['n_datasets', 'marker'], ascending=[False, True])
    W = W.loc[F.index]
    F.index.name = 'uniprot'
    W.index.name = 'uniprot'

    META['n_proteins_in_master_matrix'] = W.notna().sum(axis=0).reindex(META['sample_uid']).values
    FULL['n_proteins_in_master_matrix'] = META['n_proteins_in_master_matrix'].values

    # write
    LONG.to_csv(os.path.join(OUT, 'master_long.csv.gz'), index=False, compression='gzip')
    W.to_csv(os.path.join(OUT, 'master_matrix.csv.gz'), compression='gzip')
    F.to_csv(os.path.join(OUT, 'master_features.csv'))
    META.to_csv(os.path.join(OUT, 'master_metadata.csv'), index=False)
    FULL.to_csv(os.path.join(OUT, 'master_metadata_full.csv'), index=False)

    try:
        import anndata as ad
        obs = META.set_index('sample_uid').copy()
        for c in obs.columns:
            if obs[c].dtype == object:
                obs[c] = obs[c].astype('string').fillna('NA').astype(str)
        var = F.copy()
        for c in var.columns:
            if var[c].dtype == object:
                var[c] = var[c].astype(str)
        X = W.T.loc[obs.index, var.index].to_numpy(dtype='float64')
        A = ad.AnnData(X=X, obs=obs, var=var)
        A.layers['measured'] = (~np.isnan(X)).astype('int8')
        A.uns['description'] = ('Olink master object. X = values exactly as published (NOT normalized; scales differ by dataset - '
                                'see obs.value_scale). NaN = protein not measured in that dataset. One assay per protein per dataset; '
                                'all assays are in master_long.csv.gz.')
        A.write_h5ad(os.path.join(OUT, 'master_object.h5ad'), compression='gzip')
        print('wrote master_object.h5ad', A.shape)
    except ImportError:
        print('anndata not installed - skipped .h5ad')

    print('samples', len(META), 'patients', META['patient_uid'].nunique(), 'proteins', W.shape[0], 'long rows', len(LONG),
          'values in matrix', int(W.notna().sum().sum()))


if __name__ == '__main__':
    main()
