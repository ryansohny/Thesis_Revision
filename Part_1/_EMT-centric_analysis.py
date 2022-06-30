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

# Cristescu
p = sns.violinplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index_Cristescu", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette=type_palette)
p = sns.swarmplot(data=new_tcga_ws76, x="gene_expression_subtype", y="EMT-index_Cristescu", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".5", size=3, hue="EMT Groups (Cristescu et al.)", palette=dot_palette_cristescu)
p.set_xticklabels(['Proliferative', 'Differentiated', 'Immunoreactive', 'Mesenchymal'])
p.legend(title='Group defined by method\nfrom Cristescu et al.', loc='upper center', frameon=False)
p.set_xlabel("TCGA HGSOC subtypes")
p.set_ylabel("EMT index from Cristescu et al.")
sns.despine()
plt.tight_layout()

# Guo
p = sns.violinplot(data=new_tcga_ws76, x="gene_expression_subtype", y="Norm_WS_EMT", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette=type_palette)
p = sns.swarmplot(data=new_tcga_ws76, x="gene_expression_subtype", y="Norm_WS_EMT", order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".5", size=3, hue="EMT Groups (Guo et al.)", palette=dot_palette_guo)
p.set_xticklabels(['Proliferative', 'Differentiated', 'Immunoreactive', 'Mesenchymal'])
p.legend(title='Group defined by method\nfrom Guo et al.', loc='lower center', frameon=False)
p.set_xlabel("TCGA HGSOC subtypes")
p.set_ylabel("EMT index from Guo et al.")
sns.despine()
plt.tight_layout()
