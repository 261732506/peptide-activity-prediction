"""
Generate Supplementary Figure S1: ESM-2 Length Bias Analysis
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# Read data
df = pd.read_csv('point_mutation_candidates_ESM2_RF_REAL.csv')

# Create figure with 2 subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Scatter plot with regression line
ax1 = axes[0]
x = df['length']
y = df['ESM2_RF_Joint_prob']

# Scatter plot
ax1.scatter(x, y, alpha=0.5, s=30, c='#3498db', edgecolors='white', linewidth=0.5)

# Regression line
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
line_x = np.array([x.min(), x.max()])
line_y = slope * line_x + intercept
ax1.plot(line_x, line_y, 'r--', linewidth=2, label=f'Linear fit: r={r_value:.3f}')

# Add vertical line at 16 aa threshold
ax1.axvline(x=16, color='green', linestyle='--', linewidth=2, alpha=0.7,
            label='Length filter (≥16 aa)')

# Annotations
ax1.text(0.05, 0.95, f'Pearson r = {r_value:.3f}\np < 0.0001',
         transform=ax1.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax1.set_xlabel('Sequence Length (amino acids)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Joint Probability (CPP × AMP)', fontsize=11, fontweight='bold')
ax1.set_title('A. ESM-2 Length Bias in Point Mutations', fontsize=12, fontweight='bold')
ax1.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
ax1.grid(True, alpha=0.3)

# Panel B: Box plot comparison
ax2 = axes[1]

# Group by length
length_groups = []
joint_probs = []
colors_box = []

for length in sorted(df['length'].unique()):
    subset = df[df['length'] == length]
    length_groups.append(f'{length} aa')
    joint_probs.append(subset['ESM2_RF_Joint_prob'].values)
    colors_box.append('#e74c3c' if length < 16 else '#27ae60')

# Box plot
bp = ax2.boxplot(joint_probs, labels=length_groups, patch_artist=True,
                  widths=0.6, showfliers=True,
                  boxprops=dict(linewidth=1.5),
                  whiskerprops=dict(linewidth=1.5),
                  capprops=dict(linewidth=1.5),
                  medianprops=dict(color='darkblue', linewidth=2))

# Color boxes
for patch, color in zip(bp['boxes'], colors_box):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

# Add modular assembly reference line
modular_max = 0.367
ax2.axhline(y=modular_max, color='purple', linestyle='--', linewidth=2,
            label=f'Modular max = {modular_max:.3f}')

# Add filtered point mutation max
filtered_max = df[df['length'] >= 16]['ESM2_RF_Joint_prob'].max()
ax2.axhline(y=filtered_max, color='green', linestyle=':', linewidth=2,
            label=f'Filtered point mutation max = {filtered_max:.3f}')

ax2.set_xlabel('Sequence Length', fontsize=11, fontweight='bold')
ax2.set_ylabel('Joint Probability (CPP × AMP)', fontsize=11, fontweight='bold')
ax2.set_title('B. Joint Probability Distribution by Length', fontsize=12, fontweight='bold')
ax2.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
ax2.grid(True, alpha=0.3, axis='y')
ax2.set_xticklabels(length_groups, rotation=45, ha='right')

# Add color legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#e74c3c', alpha=0.6, label='Excluded (<16 aa)'),
                   Patch(facecolor='#27ae60', alpha=0.6, label='Retained (≥16 aa)')]
ax2.legend(handles=legend_elements + [bp['medians'][0]],
           labels=['Excluded (<16 aa)', 'Retained (≥16 aa)',
                   f'Modular max = {modular_max:.3f}',
                   f'Filtered max = {filtered_max:.3f}'],
           loc='upper right', frameon=True, fancybox=True, shadow=True)

plt.tight_layout()

# Save figure
output_file = 'ESM2_length_bias_analysis.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f'Supplementary Figure S1 saved: {output_file}')

# Also save as PDF for publication
plt.savefig('ESM2_length_bias_analysis.pdf', dpi=300, bbox_inches='tight')
print(f'PDF version saved: ESM2_length_bias_analysis.pdf')

# Generate statistics summary
print('\n' + '='*70)
print('Statistical Summary for Supplementary Figure S1')
print('='*70)

print('\nLength-Probability Correlation:')
print(f'  Pearson r: {r_value:.3f}')
print(f'  p-value: {p_value:.2e}')
print(f'  95% CI: [{r_value - 1.96*std_err:.3f}, {r_value + 1.96*std_err:.3f}]')

print('\nMean Joint Probability by Length Group:')
for length in sorted(df['length'].unique()):
    subset = df[df['length'] == length]
    mean_joint = subset['ESM2_RF_Joint_prob'].mean()
    std_joint = subset['ESM2_RF_Joint_prob'].std()
    n = len(subset)
    status = 'EXCLUDED' if length < 16 else 'RETAINED'
    print(f'  {length} aa: {mean_joint:.3f} ± {std_joint:.3f} (n={n}) [{status}]')

print('\nComparison:')
excluded = df[df['length'] < 16]['ESM2_RF_Joint_prob']
retained = df[df['length'] >= 16]['ESM2_RF_Joint_prob']
print(f'  Excluded (<16 aa): mean={excluded.mean():.3f}, max={excluded.max():.3f}, n={len(excluded)}')
print(f'  Retained (≥16 aa): mean={retained.mean():.3f}, max={retained.max():.3f}, n={len(retained)}')
print(f'  Modular Assembly: max={modular_max:.3f}')
print(f'  Modular > Filtered Point Mutation: +{(modular_max-retained.max())/retained.max()*100:.1f}%')

print('\n' + '='*70)
