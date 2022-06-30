# cm03 wc300_ver2 environment
# ipython --profile=thesis

import matplotlib_venn
from matplotlib_venn import venn3
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy
import scipy.stats as stats
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import chi2_contingency
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests

plt.rcParams['figure.figsize'] = (7,7)
sns.set(font="Arial", font_scale=1.2, style='ticks')
%autoindent
%matplotlib

# Version Info
print(matplotlib_venn.__version__) # 0.11.7
print(sns.__version__) # 0.11.2
print(matplotlib.__version__) # 3.4.3
print(pd.__version__) # 1.3.5
print(scipy.__version__) # 1.8.1

def get_asterisks_for_pval(p_val):
    """Receives the p-value and returns asterisks string."""
    if p_val > 0.05:
        p_text = "ns"  # above threshold => not significant
    elif p_val < 1e-4:  
        p_text = '****'
    elif p_val < 1e-3:
        p_text = '***'
    elif p_val < 1e-2:
        p_text = '**'
    else:
        p_text = '*'
    
    return p_text

def chisq_and_posthoc_corrected(df):
    """Receives a dataframe and performs chi2 test and then post hoc.
    Prints the p-values and corrected p-values (after FDR correction)"""
    # start by running chi2 test on the matrix
    chi2, p, dof, ex = chi2_contingency(df, correction=True)
    print(f"Chi2 result of the contingency table: {chi2}, p-value: {p}")
    
    # post-hoc
    all_combinations = list(combinations(df.index, 2))  # gathering all combinations for post-hoc chi2
    p_vals = []
    print("Significance results:")
    for comb in all_combinations:
        new_df = df[(df.index == comb[0]) | (df.index == comb[1])]
        chi2, p, dof, ex = chi2_contingency(new_df, correction=True)
        p_vals.append(p)

    # checking significance
    # correction for multiple testing
    reject_list, corrected_p_vals = multipletests(p_vals, method='fdr_bh')[:2]
    for p_val, corr_p_val, reject, comb in zip(p_vals, corrected_p_vals, reject_list, all_combinations):
        print(f"{comb}: p_value: {p_val}; corrected: {corr_p_val} ({get_asterisks_for_pval(p_val)}) reject: {reject}")

sohn_emtindex = pd.read_table("EMT_Index_HGSOC.txt", index_col=0, sep='\t')
sohn_emtindex['Cluster'].replace(to_replace={"B":"Cluster B", "A": "Cluster A"}, inplace=True)
sohn_emtindex.rename(columns={'Cluster':'HGSOC Cluster'}, inplace=True)
sohn_tcga_emtindex = pd.read_table("TCGA_EMT-index.txt", index_col=0, sep='\t')
sohn_tcga_emtindex.rename(columns={'EMT_class':'Defined by EMT index'}, inplace=True)
sohn_tcga_emtindex['Defined by EMT index'].replace(to_replace={"EMT_high":"EMT-high", "EMT_low": "EMT-low"}, inplace=True)

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2]
    return lst3 

emt_sohn = list(map(lambda x: x.strip(), open("EMT_Genes.txt", 'r').readlines()))[1:]           
emt_cristescu = list(map(lambda x: x.strip(), open("EMT_Genes_Cristescu.txt", 'r').readlines()))[1:]
emt_byers = list(map(lambda x: x.strip(), open("EMT_Genes_Byers.txt", 'r').readlines()))[1:]

# Venn Diagram (cm03 wc300_ver2 environment)
vd = venn3(subsets = (35, 70, 0, 141, 2, 5, 1), set_labels = ("Sohn's Thesis", "Byers et al., 2013", "Cristescu et al., 2015"), alpha = 0.7)
#plt.title("# of Overlapped Genes")
#plt.tight_layout()

# Cristescu gene signatures on TPM (SNUH)
df_sohn = pd.read_table("HGOC_TPM_new_GeneLevel_geneSymbol_onlyTumor.txt", sep="\t", index_col=0)
"""
# Checking for missing gene 
for i in emt_cristescu: 
    if i not in list(df_sohn[df_sohn.index.isin(emt_cristescu)].index):
        print(i) 

"""

df_sohn_cristescu = df_sohn[df_sohn.index.isin(emt_cristescu)]
"""
df_sohn_cristescu.to_csv("EMT_Genes_Cristescu_SNUH.txt", sep="\t")
Then, Run _calculate_geometric_mean_from_TPM.py ==> GM_EMT_Genes_Cristescu_SNUH.txt
"""
emt_index_from_cristescu = pd.read_table("GM_EMT_Genes_Cristescu_SNUH.txt", index_col=0, sep="\t")
sohn_ws76['EMT-index_Cristescu'] = emt_index_from_cristescu.loc['EMT_index'].values # From 76 gene sigantures table

label_orders = ['Cluster A', 'Cluster B']
g = sns.lmplot(data=sohn_ws76, x='EMT-index', y='EMT-index_Cristescu', hue='HGSOC Cluster', palette={'Cluster A':'midnightblue', 'Cluster B':'darkred'}, fit_reg=False, facet_kws={'legend_out': True})
for text, label in zip(g._legend.texts, label_orders):
    text.set_text(label)

sns.regplot(data=sohn_ws76, x='EMT-index', y='EMT-index_Cristescu', scatter=False, color='black', ax=g.axes[0,0])
g.set_xlabels("EMT index")
g.set_ylabels("EMT index from Cristescu et al.")
stats.pearsonr(sohn_ws76['EMT-index_Cristescu'], sohn_ws76['EMT-index'])

# 76 gene signatures on TPM (SNUH)
df_sohn = pd.read_table("HGOC_TPM_new_GeneLevel_geneSymbol_onlyTumor.txt", sep="\t", index_col=0)
"""
# Checking for missing gene 
for i in emt_byers: 
    if i not in list(df_sohn[df_sohn.index.isin(emt_byers)].index):
        print(i) 

"""

df_sohn_byers = df_sohn[df_sohn.index.isin(emt_byers)]
cdh1 = df_sohn_byers.apply(lambda x: stats.pearsonr(x.values, df_sohn_byers.loc['CDH1'].values)[0], axis=1)
df_sohn_byers['PCC'] = cdh1.copy() # Warnings https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
ws_emt = pd.DataFrame(list(map(lambda x: (df_sohn_byers.iloc[: , x] *  df_sohn_byers.iloc[:, -1]).sum(), list(range(0,20)))), index=df_sohn_byers.iloc[:,:-1].columns, columns=['WS_EMT']) # Weighted Sum of EMT
ws_emt['Norm_WS_EMT'] = (ws_emt - ws_emt.mean()).values.flatten() # Normalized Weighted Sum of EMT
sohn_ws76 = pd.concat([ws_emt, sohn_emtindex], axis=1)

label_orders = ['Cluster A', 'Cluster B']
g = sns.lmplot(data=sohn_ws76, x='EMT-index', y='Norm_WS_EMT', hue='HGSOC Cluster', palette={'Cluster A':'midnightblue', 'Cluster B':'darkred'}, fit_reg=False, facet_kws={'legend_out': True})
for text, label in zip(g._legend.texts, label_orders):
    text.set_text(label)

sns.regplot(data=sohn_ws76, x='EMT-index', y='Norm_WS_EMT', scatter=False, color='black', ax=g.axes[0,0])
g.set_xlabels("EMT index")
g.set_ylabels("EMT score from Guo et al.")
stats.pearsonr(sohn_ws76['Norm_WS_EMT'], sohn_ws76['EMT-index'])

# 76 gene signatures on log2TPM (SNUH)
df_sohn_log2 = pd.read_table("HGOC_TPM_new_GeneLevel_geneSymbol_onlyTumor_log2.txt", sep="\t", index_col=0)
df_sohn_log2_byers = df_sohn_log2[df_sohn_log2.index.isin(emt_byers)]
cdh1_log2 = df_sohn_log2_byers.apply(lambda x: stats.pearsonr(x.values, df_sohn_log2_byers.loc['CDH1'].values)[0], axis=1)
df_sohn_log2_byers['PCC'] = cdh1_log2.copy() # Warnings https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
ws_emt_log2 = pd.DataFrame(list(map(lambda x: (df_sohn_log2_byers.iloc[: , x] *  df_sohn_log2_byers.iloc[:, -1]).sum(), list(range(0,20)))), index=df_sohn_log2_byers.iloc[:,:-1].columns, columns=['WS_EMT']) # Weighted Sum of EMT
ws_emt_log2['Norm_WS_EMT'] = (ws_emt_log2 - ws_emt_log2.mean()).values.flatten() # Normalized Weighted Sum of EMT
sohn_ws76 = pd.concat([ws_emt_log2, sohn_emtindex], axis=1)

g = sns.lmplot(data=sohn_ws76, x='EMT-index', y='Norm_WS_EMT')
g.set_xlabels("EMT index")
g.set_ylabels("EMT score from Guo et al.")
stats.pearsonr(sohn_ws76['Norm_WS_EMT'], sohn_ws76['EMT-index'])

# 76 gene signatures on TPM (TCGA)

df_tcga = pd.read_table("TCGA-OV_TPM_GeneLevel_geneID.txt", sep="\t", index_col=0)
emt_byers = list(map(lambda x: x.strip(), open("EMT_Genes_Byers_forTCGA.txt", 'r').readlines()))[1:]
df_tcga_byers = df_tcga[df_tcga.index.isin(emt_byers)]
"""
ADGRF1  GPR110
ADGRG1  GPR56
PATJ    INADL
"""
cdh1_tcga = df_tcga_byers.apply(lambda x: stats.pearsonr(x.values, df_tcga_byers.loc['CDH1'].values)[0], axis=1)
df_tcga_byers['PCC'] = cdh1_tcga.copy() # Warnings https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
ws_emt_tcga = pd.DataFrame(list(map(lambda x: (df_tcga_byers.iloc[: , x] *  df_tcga_byers.iloc[:, -1]).sum(), list(range(0,379)))), index=df_tcga_byers.iloc[:,:-1].columns, columns=['WS_EMT']) # Weighted Sum of EMT
ws_emt_tcga['Norm_WS_EMT'] = (ws_emt_tcga - ws_emt_tcga.mean()).values.flatten() # Normalized Weighted Sum of EMT
tcga_ws76 = pd.concat([ws_emt_tcga, sohn_tcga_emtindex], axis=1)

label_orders = ['EMT-high', 'EMT-low']
g = sns.lmplot(data=tcga_ws76, x='EMT-index', y='Norm_WS_EMT', hue='Defined by EMT index', palette={'EMT-low':'midnightblue', 'EMT-high':'darkred'}, fit_reg=False, facet_kws={'legend_out':True})
for text, label in zip(g._legend.texts, label_orders):
    text.set_text(label)

sns.regplot(data=tcga_ws76, x='EMT-index', y='Norm_WS_EMT', scatter=False, color='black', ax=g.axes[0,0])
g.set_xlabels("EMT index")
g.set_ylabels("EMT score from Guo et al.")
stats.pearsonr(tcga_ws76['Norm_WS_EMT'], tcga_ws76['EMT-index'])

plt.rcParams['figure.figsize'] = (10,10)
ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
ax3 = plt.subplot(2,2,3)
ax4 = plt.subplot(2,2,4)
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['CDH1'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1), x='Norm_WS_EMT', y='CDH1', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax1) # E-cadherin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['CDH2'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1), x='Norm_WS_EMT', y='CDH2', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax2) # N-cadherin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['VIM'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1), x='Norm_WS_EMT', y='VIM', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax3) # Vimentin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['TGFB1'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1), x='Norm_WS_EMT', y='TGFB1', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax4) # TGFB

ax1.set_xlabel("EMT score from Guo et al.")
ax1.set_ylabel("Log2[TPM]")
ax1.set_title("E-cadherin ($\it{CDH1}$)")

ax2.set_xlabel("EMT score from Guo et al.")
ax2.set_ylabel("Log2[TPM]")
ax2.set_title("N-cadherin ($\it{CDH2}$)")

ax3.set_xlabel("EMT score from Guo et al.")
ax3.set_ylabel("Log2[TPM]")
ax3.set_title("Vimentin ($\it{VIM}$)")

ax4.set_xlabel("EMT score from Guo et al.")
ax4.set_ylabel("Log2[TPM]")
ax4.set_title(r"TGF$\beta$ ($\it{TGFB1}$)")

sns.despine()
plt.tight_layout()

data = pd.concat([np.log2(df_tcga+1).loc['CDH1'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1)
stats.pearsonr(data['CDH1'], data['Norm_WS_EMT'])
data = pd.concat([np.log2(df_tcga+1).loc['CDH2'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1)
stats.pearsonr(data['CDH2'], data['Norm_WS_EMT'])
data = pd.concat([np.log2(df_tcga+1).loc['VIM'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1)
stats.pearsonr(data['VIM'], data['Norm_WS_EMT'])
data = pd.concat([np.log2(df_tcga+1).loc['TGFB1'], tcga_ws76.loc[:, 'Norm_WS_EMT']], axis=1)
stats.pearsonr(data['TGFB1'], data['Norm_WS_EMT'])

plt.rcParams['figure.figsize'] = (7,7)

# Cristescu et al (TCGA)

df_tcga = pd.read_table("TCGA-OV_TPM_GeneLevel_geneID.txt", sep="\t", index_col=0)
emt_byers = list(map(lambda x: x.strip(), open("EMT_Genes_Byers_forTCGA.txt", 'r').readlines()))[1:]
df_tcga_byers = df_tcga[df_tcga.index.isin(emt_byers)]
"""
C10orf38	FAM171A1
C10orf56	ZCCHC24
C16orf45	BMERB1
C9orf19	GLIPR2
CTGF	CCN2
DFNA5	GSDME
FAM101B	RFLNB
GLT25D2	COLGALT2
KIRREL	KIRREL1
LEPREL1	P3H2
LHFP	LHFPL6
PTRF	CAVIN1
TMSL8	TMSB15A
ZNF788	ZNF788P
python _pick_specific_genes_from_list.py EMT_Genes_Cristescu_forTCGA.txt
Then, Run python _calculate_geometric_mean_from_TPM.py EMT_Genes_Cristescu_forTCGA_TCGA-OV.txt ==> GM_EMT_Genes_Cristescu_forTCGA_TCGA-OV.txt
"""
emt_index_from_cristescu_tcga = pd.read_table("GM_EMT_Genes_Cristescu_forTCGA_TCGA-OV.txt", index_col=0, sep="\t")
tcga_ws76['EMT-index_Cristescu'] = emt_index_from_cristescu_tcga.loc['EMT_index'].values # From 76 gene sigantures table

label_orders = ['EMT-high', 'EMT-low']
g = sns.lmplot(data=tcga_ws76, x='EMT-index', y='EMT-index_Cristescu', hue='Defined by EMT index', palette={'EMT-low':'midnightblue', 'EMT-high':'darkred'}, fit_reg=False, facet_kws={'legend_out':True})
for text, label in zip(g._legend.texts, label_orders):
    text.set_text(label)

sns.regplot(data=tcga_ws76, x='EMT-index', y='EMT-index_Cristescu', scatter=False, color='black', ax=g.axes[0,0])
g.set_xlabels("EMT index")
g.set_ylabels("EMT index from Cristescu et al.")
stats.pearsonr(tcga_ws76['EMT-index_Cristescu'], tcga_ws76['EMT-index'])

plt.rcParams['figure.figsize'] = (10,10)
ax1 = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
ax3 = plt.subplot(2,2,3)
ax4 = plt.subplot(2,2,4)
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['CDH1'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1), x='EMT-index_Cristescu', y='CDH1', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax1) # E-cadherin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['CDH2'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1), x='EMT-index_Cristescu', y='CDH2', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax2) # N-cadherin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['VIM'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1), x='EMT-index_Cristescu', y='VIM', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax3) # Vimentin
sns.regplot(data=pd.concat([np.log2(df_tcga+1).loc['TGFB1'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1), x='EMT-index_Cristescu', y='TGFB1', scatter_kws={'color':'black'}, line_kws={'color': 'black'}, ax=ax4) # TGFB

ax1.set_xlabel("EMT score from Cristescu et al.")
ax1.set_ylabel("Log2[TPM]")
ax1.set_title("E-cadherin ($\it{CDH1}$)")

ax2.set_xlabel("EMT score from Cristescu et al.")
ax2.set_ylabel("Log2[TPM]")
ax2.set_title("N-cadherin ($\it{CDH2}$)")

ax3.set_xlabel("EMT score from Cristescu et al.")
ax3.set_ylabel("Log2[TPM]")
ax3.set_title("Vimentin ($\it{VIM}$)")

ax4.set_xlabel("EMT score from Cristescu et al.")
ax4.set_ylabel("Log2[TPM]")
ax4.set_title(r"TGF$\beta$ ($\it{TGFB1}$)")

sns.despine()
plt.tight_layout()

data = pd.concat([np.log2(df_tcga+1).loc['CDH1'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1)
stats.pearsonr(data['CDH1'], data['EMT-index_Cristescu'])
data = pd.concat([np.log2(df_tcga+1).loc['CDH2'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1)
stats.pearsonr(data['CDH2'], data['EMT-index_Cristescu'])
data = pd.concat([np.log2(df_tcga+1).loc['VIM'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1)
stats.pearsonr(data['VIM'], data['EMT-index_Cristescu'])
data = pd.concat([np.log2(df_tcga+1).loc['TGFB1'], tcga_ws76.loc[:, 'EMT-index_Cristescu']], axis=1)
stats.pearsonr(data['TGFB1'], data['EMT-index_Cristescu'])

plt.rcParams['figure.figsize'] = (7,7)

# For Guo et al. and Cristescu et al.
label_orders = ['EMT-high', 'EMT-low']
g = sns.lmplot(data=tcga_ws76, x='EMT-index_Cristescu', y='Norm_WS_EMT', hue='Defined by EMT index', palette={'EMT-low':'midnightblue', 'EMT-high':'darkred'}, fit_reg=False, facet_kws={'legend_out':True})
for text, label in zip(g._legend.texts, label_orders):
    text.set_text(label)

sns.regplot(data=tcga_ws76, x='EMT-index_Cristescu', y='Norm_WS_EMT', scatter=False, color='black', ax=g.axes[0,0])
g.set_xlabels("EMT index from Cristescu et al.")
g.set_ylabels("EMT index from Guo et al.")
stats.pearsonr(tcga_ws76['EMT-index_Cristescu'], tcga_ws76['Norm_WS_EMT'])


# Assigning groups for survival analysis
conditions_guo = [(tcga_ws76['Norm_WS_EMT'] > 0), (tcga_ws76['Norm_WS_EMT'] < 0)]
choices_guo = ['Epithelial', 'Mesenchymal']
tcga_ws76['EMT Groups (Guo et al.)'] = np.select(conditions_guo, choices_guo, default='NA')
# sns.boxplot(data=tcga_ws76, x='EMT Groups (Guo et al.)', y='Norm_WS_EMT', palette={'Epithelial':'midnightblue', 'Mesenchymal':'darkred'}) # sanity check

conditions_cristescu = [(tcga_ws76['EMT-index_Cristescu'] > tcga_ws76['EMT-index_Cristescu'].median()), (tcga_ws76['EMT-index_Cristescu'] <= tcga_ws76['EMT-index_Cristescu'].median())]
choices_cristescu = ['EMT(+)', 'EMT(-)']
tcga_ws76['EMT Groups (Cristescu et al.)'] = np.select(conditions_cristescu, choices_cristescu, default='NA')

# TCGA subtypes
tcga_subtypes = pd.read_csv("TCGA_OV_EMThighlow_TCGAsubtypes.csv", index_col=0)

tcga_ws76_reindex = list(map(lambda x: x[:-1], list(tcga_ws76.index))) # reindexing tcga_ws76 to match id with that of tcga_subtypes
tcga_ws76.index = tcga_ws76_reindex

new_tcga_ws76 = pd.concat([tcga_ws76, pd.DataFrame(tcga_subtypes.loc[:, 'gene_expression_subtype'])], axis=1)

type_palette = {'proliferative': '#4B0082', 'differentiated': '#CD853F', 'immunoreactive': '#228B22', 'mesenchymal': '#800000'}
dot_palette_sohn = {'EMT-high': '#FF0000', 'EMT-low': '#0000FF'}
dot_palette_cristescu = {'EMT(+)': '#FF0000', 'EMT(-)': '#0000FF'}
dot_palette_guo = {'Mesenchymal': '#FF0000', 'Epithelial': '#0000FF'}

# Sohn
p = sns.violinplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette=type_palette)
p = sns.swarmplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".5", size=3, hue="Defined by EMT index", palette=dot_palette_sohn)
sns.despine()

# Cristescu et al.
p = sns.violinplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index_Cristescu", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette=type_palette)
p = sns.swarmplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index_Cristescu", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".5", size=3, hue="EMT Groups (Cristescu et al.)", palette=dot_palette_cristescu)
p.set_xticklabels(['Proliferative', 'Differentiated', 'Immunoreactive', 'Mesenchymal'])
p.legend(title='Group defined by method\nfrom Cristescu et al.', loc='upper center', frameon=False)
p.set_xlabel("TCGA HGSOC subtypes")
p.set_ylabel("EMT index from Cristescu et al.")
sns.despine()
plt.tight_layout()

dat = new_tcga_ws76.copy()
dat = dat[~dat['gene_expression_subtype'].isna()]
dat.rename(columns={'EMT-index_Cristescu': 'EMT_index_Cristescu'}, inplace=True)
dat2 = dat.reset_index()

model = ols('EMT_index_Cristescu ~ C(gene_expression_subtype)', dat).fit() # One-way ANOVA
anova_lm(model)
#                               df       sum_sq     mean_sq          F        PR(>F)
#C(gene_expression_subtype)    3.0  2786.875126  928.958375  76.430113  1.400172e-35
#Residual                    263.0  3196.594143   12.154350        NaN           NaN

posthoc = pairwise_tukeyhsd(dat2['EMT_index_Cristescu'], dat2['gene_expression_subtype'], alpha=0.05)
print(posthoc)
#        Multiple Comparison of Means - Tukey HSD, FWER=0.05         
#====================================================================
#    group1         group2     meandiff p-adj   lower   upper  reject
#--------------------------------------------------------------------
#differentiated immunoreactive  -0.0196    1.0 -1.5716  1.5323  False
#differentiated    mesenchymal   8.1726   -0.0  6.5332  9.8121   True
#differentiated  proliferative   1.6942 0.0276  0.1325   3.256   True
#immunoreactive    mesenchymal   8.1923   -0.0   6.619  9.7655   True
#immunoreactive  proliferative   1.7139 0.0171  0.2217   3.206   True
#   mesenchymal  proliferative  -6.4784   -0.0 -8.0613 -4.8955   True

contingency_dat = dat.groupby(['gene_expression_subtype', 'EMT Groups (Cristescu et al.)']).size().unstack()
print(chisq_and_posthoc_corrected(contingency_dat))
#Chi2 result of the contingency table: 51.02165774042108, p-value: 4.8403848915596475e-11
#Significance results:
#('differentiated', 'immunoreactive'): p_value: 0.28694231985919705; corrected: 0.3282628746122071 (ns) reject: False
#('differentiated', 'mesenchymal'): p_value: 1.222374921665085e-08; corrected: 3.6671247649952555e-08 (****) reject: True
#('differentiated', 'proliferative'): p_value: 0.3282628746122071; corrected: 0.3282628746122071 (ns) reject: False
#('immunoreactive', 'mesenchymal'): p_value: 8.50421565753058e-12; corrected: 5.102529394518348e-11 (****) reject: True
#('immunoreactive', 'proliferative'): p_value: 0.020741076521377178; corrected: 0.03111161478206577 (*) reject: True
#('mesenchymal', 'proliferative'): p_value: 1.1232482012425821e-06; corrected: 2.2464964024851642e-06 (****) reject: True

# Guo et al.
p = sns.violinplot(data=new_tcga_ws76, x="gene_expression_subtype", y="Norm_WS_EMT", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette=type_palette)
p = sns.swarmplot(data=new_tcga_ws76, x="gene_expression_subtype", y="Norm_WS_EMT", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".5", size=3, hue="EMT Groups (Guo et al.)", palette=dot_palette_guo)
p.set_xticklabels(['Proliferative', 'Differentiated', 'Immunoreactive', 'Mesenchymal'])
p.legend(title='Group defined by method\nfrom Guo et al.', loc='lower center', frameon=False)
p.set_xlabel("TCGA HGSOC subtypes")
p.set_ylabel("EMT score from Guo et al.")
sns.despine()
plt.tight_layout()

model = ols('Norm_WS_EMT ~ C(gene_expression_subtype)', dat).fit() # One-way ANOVA
anova_lm(model)
#                               df        sum_sq        mean_sq         F    PR(>F)
#C(gene_expression_subtype)    3.0  1.384702e+06  461567.499537  2.666874  0.048205
#Residual                    263.0  4.551855e+07  173074.343414       NaN       NaN

posthoc = pairwise_tukeyhsd(dat2['Norm_WS_EMT'], dat2['gene_expression_subtype'], alpha=0.05)
print(posthoc)
#          Multiple Comparison of Means - Tukey HSD, FWER=0.05           
#========================================================================
#    group1         group2      meandiff p-adj    lower    upper   reject
#------------------------------------------------------------------------
#differentiated immunoreactive -172.2362 0.0787 -357.4289  12.9565  False
#differentiated    mesenchymal -139.3145 0.2564 -334.9452  56.3162  False
#differentiated  proliferative -182.8411 0.0567 -369.2027   3.5206  False
#immunoreactive    mesenchymal   32.9217 0.9689 -154.8154 220.6588  False
#immunoreactive  proliferative  -10.6049 0.9987 -188.6625 167.4527  False
#   mesenchymal  proliferative  -43.5266 0.9333 -232.4169 145.3637  False
#------------------------------------------------------------------------

contingency_dat = dat.groupby(['gene_expression_subtype', 'EMT Groups (Guo et al.)']).size().unstack()
print(chisq_and_posthoc_corrected(contingency_dat))
#Chi2 result of the contingency table: 8.706193232309863, p-value: 0.03346333724992227
#Significance results:
#('differentiated', 'immunoreactive'): p_value: 0.010681446593179485; corrected: 0.06408867955907692 (*) reject: False
#('differentiated', 'mesenchymal'): p_value: 0.24433070837670517; corrected: 0.38356677746108875 (ns) reject: False
#('differentiated', 'proliferative'): p_value: 0.02640650739286577; corrected: 0.07921952217859732 (*) reject: False
#('immunoreactive', 'mesenchymal'): p_value: 0.25571118497405915; corrected: 0.38356677746108875 (ns) reject: False
#('immunoreactive', 'proliferative'): p_value: 0.8611818382885372; corrected: 0.8611818382885372 (ns) reject: False
#('mesenchymal', 'proliferative'): p_value: 0.4190575275621842; corrected: 0.5028690330746209 (ns) reject: False
