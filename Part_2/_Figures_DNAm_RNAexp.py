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

clinic_info = pd.read_csv("2022_WC300_clinical_information_Xadded_ver2.0.csv", index_col=0)

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

# Super-enhancer analysis
# Super Enhancer Overlapped Genes
se_genes = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/SE_DMR_abs15_hyper_overlapped_genes.txt")
se_genes = list(se_genes.values.flatten()) # 100개

deg_tn[deg_tn.index.isin(se_genes)] # ==> 83개
deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] > 0)].index
deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] < 0)].index

# deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ==> 12개
# deg_tn_downlist[deg_tn_downlist.isin(se_genes)] ==> 24개
se_deg = list( deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ) + list( deg_tn_downlist[deg_tn_downlist.isin(se_genes)] )

col_colors1 = ['#C0C0C0']*84 + ['#000000']*84
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057B8", "#000000", "#ffd700"])
cmap_rkg = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#00FF00", "#000000", "#ff0000"])   #"R"ed blac"k" "G"reen

g = sns.clustermap(gene_vst[gene_vst.index.isin( se_deg )],
                   col_cluster=False,
                    method='ward',
                   metric='euclidean',
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap_rkg,
                    col_colors=[col_colors1],
                   xticklabels=False,
                    yticklabels=True, vmin=-2.5, vmax=2.5)
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')

se_deg_rna_roworder = g.dendrogram_row.reordered_ind

#  Array of codes for making SE_DEG_ALL.txt (refer to /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/DMR/DMR_min55_new/Enrichment/Stomach_SE)
reorder = list(gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index)

se_deg_met = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/SE_DEG_ALL.txt", index_col=0)
se_deg_met.columns = list(map(lambda x: 'X'+x, se_deg_met.columns))

a = pd.DataFrame(list(map(lambda x: x.split('/')[0], se_deg_met.index)), columns=['Region'], index=se_deg_met.index)
b = pd.DataFrame(list(map(lambda x: x.split('/')[1], se_deg_met.index)), columns=['GeneID'], index=se_deg_met.index)
c = pd.DataFrame(list(map(lambda x: x.split('/')[2], se_deg_met.index)), columns=['CpG'], index=se_deg_met.index)
se_deg_met_info = pd.concat([a,b,c], axis=1)
del a, b, c

reorder_iloc = list()

for i in reorder:
    reorder_iloc.append(list(se_deg_met_info['GeneID'].values).index(i))

spearmanr = list()
pvalues = list()

for i in list(range(36)):
    spearmanr.append(stats.spearmanr(se_deg_met.iloc[reorder_iloc].iloc[i], gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].iloc[i])[0])

for i in list(range(36)):
    pvalues.append(stats.spearmanr(se_deg_met.iloc[reorder_iloc].iloc[i], gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].iloc[i])[1])

from statsmodels.stats.multitest import multipletests
bh_corrected_pval = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1]

spearmanr_df = pd.DataFrame(spearmanr, index=gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index, columns=['Spearman_Rho'])
bh_corrected_pval_df = pd.DataFrame(bh_corrected_pval, index=gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index, columns=['FDR_BH_Pval'])
df_final = pd.concat([spearmanr_df, bh_corrected_pval_df], axis=1)
df_final = df_final[df_final['FDR_BH_Pval'] < 0.01].sort_values(by='Spearman_Rho', ascending=True)

# (Again) Heatmap of SE overlapped gene expression 
g = sns.clustermap(gene_vst[gene_vst.index.isin(df_final.index)].reindex(df_final.index),
                   col_cluster=False,
                   row_cluster=False,
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap_rkg,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=True, vmin=-2.5, vmax=2.5)
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')

sorter = list(df_final.index) # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
sorterIndex = dict(zip(sorter, range(len(sorter))))
new_se_deg_met_info = se_deg_met_info[se_deg_met_info['GeneID'].isin(df_final.index)]
new_se_deg_met_info['New_order'] = new_se_deg_met_info['GeneID'].map(sorterIndex)
new_se_deg_met_info.sort_values(by='New_order', inplace=True)


# DNA methylation of Super enhancers overlapped with DEG
g = sns.clustermap(se_deg_met[se_deg_met.index.isin(new_se_deg_met_info.index)].reindex(new_se_deg_met_info.index),
                   col_cluster=False,
                   row_cluster=False,
                   cmap=cmap,
                   z_score=None,
                   standard_scale=0,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=False)

sns.clustermap(df_final['Spearman_Rho'], col_cluster=False, row_cluster=False, vmin=-0.65, vmax=0.22)

# Homeobox (HOX) cluster DNA methylation
hox = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/HOX_Clusters_ALL.txt", index_col=0)
hox_fillna0 = hox.fillna(0)
hoxa = hox_fillna0.iloc[100:210].copy()
hoxb = hox_fillna0.iloc[330:].copy()
hoxc = hox_fillna0.iloc[210:330].copy()
hoxd = hox_fillna0.iloc[:100].copy()
hox_fillna0_new = pd.concat([hoxa, hoxb, hoxc, hoxd])
sns.clustermap(hoxa,
                   method='ward',
                   metric='euclidean',
                   col_cluster=False,
                   row_cluster=False,
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=None,
                   row_colors=None,
                   cbar_kws={'label': 'DNA methylation'},
                   vmin=-3, vmax=2) # Original cmap='gnuplot2'

