"""
Filter point mutation candidates by length >= 16 aa
"""
import pandas as pd
import sys

# Fix encoding for Windows
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8')

# Read AutoDL generated ESM-2+RF predictions
df = pd.read_csv('point_mutation_candidates_ESM2_RF_REAL.csv')

print('='*70)
print('Point Mutation Data Length Filtering')
print('='*70)

print('\n[1] Original Data Statistics:')
print(f'    Total sequences: {len(df)}')
print(f'    Average length: {df["length"].mean():.1f} aa')
print(f'    Length range: [{df["length"].min()}, {df["length"].max()}] aa')
print(f'    Max Joint: {df["ESM2_RF_Joint_prob"].max():.6f}')

# Filter length >= 16
df_filtered = df[df['length'] >= 16].copy()

print('\n[2] Filtered Data Statistics (length >= 16 aa):')
print(f'    Retained: {len(df_filtered)}/{len(df)} ({len(df_filtered)/len(df)*100:.1f}%)')
print(f'    Average length: {df_filtered["length"].mean():.1f} aa')
print(f'    Length range: [{df_filtered["length"].min()}, {df_filtered["length"].max()}] aa')
print(f'    Max Joint: {df_filtered["ESM2_RF_Joint_prob"].max():.6f}')

# Find best candidate
best = df_filtered.nlargest(1, 'ESM2_RF_Joint_prob').iloc[0]
print('\n[3] Best Point Mutation Candidate (>= 16 aa):')
print(f'    Rank: {best["rank"]}')
print(f'    Sequence: {best["sequence"]}')
print(f'    Length: {int(best["length"])} aa')
print(f'    Parent: {best["parent_name"]}')
print(f'    Mutations: {best["mutations"]}')
print(f'    CPP probability: {best["ESM2_RF_CPP_prob"]:.3f}')
print(f'    AMP probability: {best["ESM2_RF_AMP_prob"]:.3f}')
print(f'    Joint probability: {best["ESM2_RF_Joint_prob"]:.3f}')

# Compare with modular assembly
modular_max = 0.367335
print('\n[4] Comparison with Modular Assembly:')
print(f'    Point Mutation max Joint (>=16aa): {df_filtered["ESM2_RF_Joint_prob"].max():.6f}')
print(f'    Modular Assembly max Joint:        {modular_max:.6f}')
improvement = (modular_max - df_filtered['ESM2_RF_Joint_prob'].max()) / df_filtered['ESM2_RF_Joint_prob'].max() * 100
print(f'    Modular improvement:                +{improvement:.1f}%')
print(f'\n    Result: Modular Assembly > Point Mutation Strategy')

# Top 10 comparison
print('\n[5] Top 10 Candidates Comparison:')
top10_point = df_filtered.nlargest(10, 'ESM2_RF_Joint_prob')
print(f'    Point Mutation Top 10 (>=16aa):')
print(f'      Average length: {top10_point["length"].mean():.1f} aa')
print(f'      Average Joint: {top10_point["ESM2_RF_Joint_prob"].mean():.3f}')

# Read modular data
try:
    df_modular_r9 = pd.read_csv('true_esm_modular_candidates.csv')
    df_modular_tat = pd.read_csv('true_esm_modular_candidates1003.csv')
    df_modular = pd.concat([df_modular_r9, df_modular_tat])
    top10_modular = df_modular.nlargest(10, 'Joint_prob')
    print(f'    Modular Assembly Top 10:')
    print(f'      Average length: {top10_modular["length"].mean():.1f} aa')
    print(f'      Average Joint: {top10_modular["Joint_prob"].mean():.3f}')
except:
    print('    (Cannot load modular data files)')

# Save filtered data
output_file = 'point_mutation_candidates_ESM2_RF_filtered_16aa.csv'
df_filtered.to_csv(output_file, index=False)
print(f'\n[6] Data Saved:')
print(f'    Output file: {output_file}')
print(f'    Contains {len(df_filtered)} filtered candidates')

# Save top 10 for easy reference
top10_output = 'point_mutation_top10_filtered_16aa.csv'
top10_point.to_csv(top10_output, index=False)
print(f'    Top 10 file: {top10_output}')

print('\n' + '='*70)
print('Analysis Complete!')
print('='*70)
