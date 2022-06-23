import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
sns.set(font="Arial", font_scale=1.2, style='ticks')
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap4 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#191970", "#ffdab9", "#8B0000"])
%matplotlib
%autoindent

# Epimutation Burden-related analyses

dmr_met = pd.read_table("DMR_abs10_smooth.txt", index_col=0)
dmr_met.columns = list(map(lambda x: 'X'+x, dmr_met.columns))
dmr_met = dmr_met*100

dmr = sc.AnnData(dmr_met.iloc[:, 84:].T) # Only Tumor
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X

dmr.obs = clinic_info.iloc[84:, :].copy()

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=83, zero_center=True)
#sc.pl.pca(dmr, color='TN', palette={'Normal':'midnightblue', 'Tumor':'darkred'}, annotate_var_explained=True, size=100)
sc.pl.pca(dmr, color=['EBV'], palette={'Negative':'midnightblue', 'Positive':'darkred'}, annotate_var_explained=True, size=300)
sns.despine()
g = sc.pl.pca(dmr, color=['EpiBurden'], color_map=cmap, annotate_var_explained=True, size=300, show=False)
g.set_title("Epimutation Burden")
sns.despine()

PC1 = pd.DataFrame(list(map(lambda x: dmr.obsm['X_pca'][x][0], list(range(0,84)))), columns=['PC1'], index=dmr.obs.index)
df = pd.concat([PC1, dmr.obs['EpiBurden']], axis=1)
g = sns.lmplot(data=df, x="PC1", y="EpiBurden")
g.set_ylabels("Epimutation Burden")
stats.pearsonr(df['EpiBurden'], df['PC1'])

# TCGA subtypes
tcga_info = pd.read_csv('/data/Projects/phenomata/01.Projects/12.Thesis/Part2/tcga_subtypes.csv', index_col='Sample')
dmr.obs['TCGA Subtypes'] = tcga_info['Type']

df = pd.concat([dmr.obs['EpiBurden'], dmr.obs['TCGA Subtypes']], axis=1)
df.rename(columns={"EpiBurden": "Epimutation Burden"}, inplace=True)

p = sns.boxplot(data=df, x='TCGA Subtypes', y='Epimutation Burden', palette={'EBV': '#DC143C', 'MSI': '#4169E1', 'GS/CIN': '#9370DB'}, width=0.8, showfliers=True, order=['EBV', 'MSI', 'GS/CIN'])
p = sns.stripplot(data=df, x='TCGA Subtypes', y='Epimutation Burden', jitter=True, marker='o', color='black', size=2.0, alpha=0.5, order=['EBV', 'MSI', 'GS/CIN'])
p.set_xticklabels(['EBV (N=4)', 'MSI (N=10)', 'GS/CIN (N=70)'])
plt.tight_layout()
sns.despine()

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'EBV']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'MSI']['Epimutation Burden'])
MannwhitneyuResult(statistic=27.0, pvalue=0.3736263736263737)

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'EBV']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'GS/CIN']['Epimutation Burden'])
MannwhitneyuResult(statistic=199.0, pvalue=0.16814499237806202)

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'MSI']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'GS/CIN']['Epimutation Burden'])
MannwhitneyuResult(statistic=356.0, pvalue=0.9362267365375121)
