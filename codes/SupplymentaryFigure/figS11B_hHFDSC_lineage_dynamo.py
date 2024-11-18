import scanpy as sc ###
import dynamo as dyn
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt 
from scipy.sparse import csr_matrix
import squidpy as sq ###
import sys
import time
import os
import re
from pandas.core.frame import DataFrame

id = sys.argv[1]  # hair follicle selection id
sample = id.split('_SCT_selection')[0]
workpath="../data/figure2/Spatial/04.RCTD/dynamo/"+sample+'/'+id #working path
inh5ad=id+"_all_assays.h5ad" #input path of h5ad file with annotation
outpath="../data/figure2/dynamo/"+sample+'/'+id+"/HFDSC"
root_state=0 #"cell_in_HFSC" - 0,   "cell_in_HFDSC" - 1
root_state = int(root_state)

print(f"id: {id}")
print(f"workpath: {workpath}")
print(f"inh5ad: {inh5ad}")
print(f"outpath: {outpath}")
print(f"root_state: {root_state}")

if not os.path.exists(workpath):
    os.makedirs(workpath,exist_ok=True)
    print("file created")
else:
    print("file exist")

os.chdir(workpath)
anndata = sc.read(inh5ad)

# import data
anndata.obsm['X_spatial'] = anndata.obsm['spatial'].copy()

color_df = pd.read_csv('../data/RCTD.color_LastlEdition240612.txt', sep='\t')
co_all = dict(zip(color_df['order'], color_df['Color']))


# subset cell types
cell_in_HFSC = ["Active_HFSC",'Quiescent_HFSC','ORS_Suprabasal','ORS',"TAC_2","TAC_3","Medulla_Cortex_Cuticle","IRS_cuticle","HF_IRS"]
cell_in_HFDSC = ["hfDSC","Dermal_papilla","Dermal_sheath"]
cell_in_Epidermis = ["HF_basal","IFE_basal","IFE_granular","IFE_Spinous","IFE_Suprabasal"]
cell_list = [cell_in_HFSC,cell_in_HFDSC,cell_in_Epidermis]

anndata = anndata[anndata.obs['celltype'].isin(cell_list[root_state])]
root_cell = cell_list[root_state][0]
co = {celltype: co_all[celltype] for celltype in cell_list[root_state] if celltype in co_all} 
if not os.path.exists(outpath):
    os.makedirs(outpath,exist_ok=True)
    print("file created")
else:
    print("file exist")

os.chdir(outpath)


# QC
sc.pp.calculate_qc_metrics(anndata, expr_type='unspliced', layer = "unspliced", inplace = True)
sc.pp.calculate_qc_metrics(anndata, expr_type='spliced', layer = "spliced", inplace = True)
sc.pl.violin(anndata, ["n_genes_by_spliced", "n_genes_by_unspliced"])
plt.savefig('violin.n_genes.png')

sc.pl.violin(anndata, ["log1p_total_spliced",  "log1p_total_unspliced"])
plt.savefig('violin.log1p_total.png')
sc.pl.violin(anndata, ["total_spliced",  "total_unspliced"])
plt.savefig('violin.total.png')

print("qc plot have read")

# select genes for analysis
expressed_genes = (anndata.layers['spliced'].A + anndata.layers['unspliced'].A).sum(0) > n_gene 
expressed_genes = anndata.var_names[expressed_genes]

len(expressed_genes)

sq.gr.spatial_neighbors(
     anndata,
     n_neighs= 8
)

sq.gr.spatial_autocorr(
    anndata,
    mode = 'moran',
    genes=list(expressed_genes),
    n_perms=100,
    n_jobs=30,
)

moran_top_gene = anndata.uns["moranI"].head(2000).index
moran_top_gene = expressed_genes.intersection(moran_top_gene)
final_gene = (set(moran_top_gene) & set(expressed_genes))


len(final_gene)

# run dynamo
preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True,force_gene_list=list(expressed_genes))
preprocessor.preprocess_adata(anndata, recipe='monocle')

# preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
# preprocessor.preprocess_adata(anndata, recipe='monocle')
dyn.tl.reduceDimension(anndata, enforce = True)
dyn.pl.umap(anndata, color='celltype',show_legend= 'on data',alpha=1,background='black',pointsize=0.25,figsize=(28, 12), color_key=co) 
plt.savefig('uamp.svg')
print("Preprocessor have read!!!!!!!")

dyn.tl.dynamics(anndata, cores=10)
dyn.tl.gene_wise_confidence(anndata, group='celltype',lineage_dict = { root_cell :["Dermal_sheath","Dermal_papilla"]} )

os.chdir(outpath)
print("gene_wise_confidence have read")

dyn.pl.phase_portraits(anndata, genes=anndata.var_names[anndata.var.use_for_dynamics][:4], figsize=(6, 4), color='celltype')
plt.savefig('umap_phase_portraits_gene.svg')  


dyn.tl.reduceDimension(anndata, enforce = True)
dyn.tl.cell_velocities(anndata, method="fp", basis='spatial',enforce=True,  transition_genes = list(anndata.var_names[anndata.var.use_for_pca]))
dyn.tl.cell_wise_confidence(anndata)
dyn.tl.confident_cell_velocities(anndata, group='celltype', lineage_dict={ root_cell:["Dermal_sheath","Dermal_papilla"]},only_transition_genes=True,basis='spatial')



dyn.pl.cell_wise_vectors(anndata, color=['celltype'], basis='spatial', show_legend='on data', quiver_length=12, quiver_size=8, pointsize=0.1, show_arrowed_spines=False,figsize=(28, 12),frontier = True,color_key=co,calpha=1) #background='black'
plt.savefig('cell_wise_vectors.spatial.svg')


dyn.pl.streamline_plot(anndata, color=['celltype'], basis='spatial', show_legend='on data',quiver_length=12, quiver_size=8, pointsize=0.1, color_key=co, show_arrowed_spines=True,background='black',figsize=(28, 12),calpha=1) 
plt.savefig('cell_wise_plot.spatial.svg')


dyn.vf.VectorField(anndata,
                   basis='spatial',
                   M=100)

dyn.vf.rank_velocity_genes(anndata,
                           groups='celltype',
                           vkey="velocity_S")

dyn.vf.acceleration(anndata, basis='spatial')
dyn.ext.ddhodge(anndata, basis='spatial')
transition = DataFrame(anndata.var_names[anndata.var.use_for_transition])
transition.columns = ['gene']
transition_genes = transition['gene'].tolist()

dyn.pl.kinetic_heatmap(anndata,
                       genes=transition_genes,
                       tkey='spatial_ddhodge_potential',
                       gene_order_method='maximum',
                       mode='pseudotime',
                       color_map='inferno',
                       yticklabels=False,basis='spatial'
                      )
plt.savefig('kinetic_heatmap.spatial.pdf')

dyn.pl.kinetic_heatmap(anndata,
                       genes=transition_genes,
                       tkey='spatial_ddhodge_potential',
                       gene_order_method='maximum',
                       mode='pseudotime',
                       color_map='viridis',
                       yticklabels=True,basis='spatial',
                       figsize = (11.5, 100)
                      )
plt.savefig('kinetic_heatmap.spatial.leng.png')

heatmap = anndata.uns['kinetics_heatmap']
heatmap.to_csv('kinetics_heatmap.umap.csv')

dyn.pl.streamline_plot(anndata, color=['celltype','spatial_ddhodge_potential'], 
                       basis='spatial', 
                       show_legend='None',
                       figsize=(28, 12),
                       color_key=co,
                       quiver_length=12, quiver_size=20,
                       pointsize=0.1,
                       background='black',
                       cmap ='viridis',
                       show_arrowed_spines=True,
                       calpha=1,
                       dpi=300,
                       linewidth=1,
                       s_kwargs_dict = {"alpha": 0.4},
                       streamline_color = "#ffffff",
                       streamline_kwargs = {"arrowsize": 40}                       
                      )
plt.savefig('ddhodge_potential.celltype.spatial.svg')

dyn.cleanup(anndata)
anndata.write_h5ad(id+'_dynamo.h5ad')
