"""Cross-study marker overlap summary (run after build_processed.py)."""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

OD = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '_overlap'))
P = pd.read_csv(os.path.join(OD, 'marker_presence_by_dataset.csv'), index_col=0)
DS = [c for c in P.columns if c.startswith('row')]
SUMMARY = ['row044_Braun2025_RCCvaccine', 'row078_Barroso-Sousa2021_TNBC', 'row096_Fedak2026_MG']
GROUPS = {
    'Immuno-Oncology panel (per-sample)': ['row010_Meng2026_GastricCa', 'row012_Johnson2026_ccRCC', 'row027_Tikkanen2025_SolidTumors',
                                           'row062_Coffey2024_Myeloma',
                                           'row032_Gao2025_NSCLC', 'row033_Xiao2025_TNBC', 'row052_Tsimberidou2024_AMADEUS',
                                           'row073_Liu2022_TNBC', 'row093_OReilly2026_REVOLUTION', 'row105_Galsky2025_PORTER'],
    'Inflammation panel (per-sample)': ['row040_Plummer2025_TIMEPOINT', 'row040_Plummer2025_OUTREACH', 'row070_Nunez2023_irAE',
                                        'row103_Yu2025_LGL'],
    'Oncology II panel': ['row024_Shan2026_OvarianCAR-T'],
}


def S(k):
    return set(P.index[P[k] == 1])


def names(s):
    return sorted(P.loc[list(s), 'marker'])


lines = []
io = set.intersection(*[S(k) for k in GROUPS['Immuno-Oncology panel (per-sample)']])
inf = set.intersection(*[S(k) for k in GROUPS['Inflammation panel (per-sample)']])
onc = S('row024_Shan2026_OvarianCAR-T')
persample = [k for k in DS if k not in SUMMARY]
core_all = set.intersection(*[S(k) for k in persample])
core_no_prince = set.intersection(*[S(k) for k in persample if k != 'row072_Padron2022_PRINCE'])
lines += [f'Unique proteins across all {len(DS)} datasets: {len(P)}',
          f'Proteins in >=2 datasets: {(P.n_datasets_all >= 2).sum()}; in >=10 datasets: {(P.n_datasets_all >= 10).sum()}',
          '',
          f'Shared by all {len(GROUPS["Immuno-Oncology panel (per-sample)"])} per-sample Immuno-Oncology datasets: {len(io)} (full IO panel minus Gal-1, which is absent from the REVOLUTION and PORTER releases)',
          f'Shared by all 4 per-sample Inflammation-panel datasets: {len(inf)}',
          f'IO core ∩ Inflammation core: {len(io & inf)}: ' + ', '.join(names(io & inf)),
          f'IO core ∩ Oncology II: {len(io & onc)}: ' + ', '.join(names(io & onc)),
          f'IO ∩ Inflammation ∩ Oncology II: {len(io & inf & onc)}: ' + ', '.join(names(io & inf & onc)),
          f'Measured in every per-sample dataset except PRINCE ({len(persample)-1}): {len(core_no_prince)}: ' + ', '.join(names(core_no_prince)),
          f'Measured in every per-sample dataset incl. PRINCE (30 selected proteins): {len(core_all)}: ' + ', '.join(names(core_all)),
          f'IO core covered by CITN-10 (11 Target 96 panels): {len(io & S("row086_Glass2024_CTCL_CITN10"))}/{len(io)}; by CITN-13 (Explore): {len(io & S("row086_Glass2024_CTCL_CITN13"))}/{len(io)}',
          f'PRINCE released proteins within IO core: {len(S("row072_Padron2022_PRINCE") & io)}/30',
          '',
          'Markers per dataset (unique UniProt):']
for k in DS:
    lines.append(f'  {k}: {int(P[k].sum())}' + ('  [summary only]' if k in SUMMARY else ''))
open(os.path.join(OD, 'overlap_summary.txt'), 'w').write('\n'.join(lines) + '\n')

core = P.loc[sorted(io, key=lambda a: P.loc[a, 'marker'].upper()), ['marker']].copy()
core['also_on_inflammation_panel'] = core.index.isin(inf)
core['also_on_oncology_II'] = core.index.isin(onc)
core['n_datasets_all'] = P.loc[core.index, 'n_datasets_all']
core.to_csv(os.path.join(OD, 'IO_core_91_markers.csv'))

# heatmap of pairwise shared proteins
pw = pd.read_csv(os.path.join(OD, 'pairwise_shared_markers.csv'), index_col=0)
order = sum(GROUPS.values(), []) + ['row072_Padron2022_PRINCE', 'row086_Glass2024_CTCL_CITN10', 'row086_Glass2024_CTCL_CITN13'] + SUMMARY
pw = pw.loc[order, order]
lab = [k.split('_', 1)[1] + (' *' if k in SUMMARY else '') for k in order]
fig, ax = plt.subplots(figsize=(13, 11))
v = np.log10(pw.values.astype(float) + 1)
im = ax.imshow(v, cmap='Blues')
for i in range(len(order)):
    for j in range(len(order)):
        n = pw.values[i, j]
        ax.text(j, i, str(n), ha='center', va='center', fontsize=7, color='white' if v[i, j] > v.max() * 0.6 else '#222')
ax.set_xticks(range(len(order)))
ax.set_yticks(range(len(order)))
ax.set_xticklabels(lab, rotation=60, ha='right', fontsize=8)
ax.set_yticklabels(lab, fontsize=8)
for b in [10, 14, 15, 18]:
    ax.axhline(b - 0.5, color='#555', lw=0.8)
    ax.axvline(b - 0.5, color='#555', lw=0.8)
ax.set_title('Proteins shared between datasets (UniProt-harmonized)\n'
             'blocks: IO panel | Inflammation | Oncology II | PRINCE + CITN broad panels | * summary-only', fontsize=11)
cb = fig.colorbar(im, ax=ax, shrink=0.6)
cb.set_label('log10(shared proteins + 1)')
fig.tight_layout()
fig.savefig(os.path.join(OD, 'pairwise_shared_markers_heatmap.png'), dpi=160)
print('\n'.join(lines))
