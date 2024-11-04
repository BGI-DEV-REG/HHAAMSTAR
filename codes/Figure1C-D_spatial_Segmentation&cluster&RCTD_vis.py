import sys
sys.path.append('../')
from Tools.Segmentation import Segmentation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import os
import warnings
warnings.filterwarnings('ignore')
import squidpy as sq
import time
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2
import dynamo as dyn
sys.path.append('../')
from Tools.Spatial import *

# 01.Segmentation
inputfile = f'../data/figure1/Spatial/00.data/{id}/'
outputfile = f'../data/figure1/Spatial/01.Segmentation/{id}/'
bs = 9
ot = 0.03
md = 7
ed = 2

bs = int(bs)
ot = float(ot)
md = int(md)
ed = int(ed)

if not os.path.exists(outputfile):
    os.mkdir(outputfile)

sobj = Segmentation()
sobj.load(
    img_path = f'{inputfile}/{id}_matched_ssDNA.png',
    mRNA_path = f'{inputfile}/{id}.gem.gz',
)

plt.figure(figsize=(16,16))
plt.imshow(sobj.raw_img, 'gray')

path = f'{outputfile}/{id}_matched_ssDNA.png'
cv2.imwrite(path, sobj.raw_img , [cv2.IMWRITE_PNG_COMPRESSION, 0])

## pre processing   sobj.raw_img --> sobj.img
sobj.pre_process(threshold='auto')

path = f'{outputfile}/{id}_process_matched_ssDNA.png'
cv2.imwrite(path, sobj.img)

## watershed
sobj.watershed(
    block_size=bs,
    offset=ot,
    min_distance=md,
    expand_distance=ed
)

path = f'{outputfile}/{id}_label_ssDNA.png'
cv2.imwrite(path, sobj.label)

## save
save_path = os.path.join(outputfile)
sobj.save_scGEM(save_path = save_path, name = id)


# 03.Spatial_cluster
inh5ad = f'../data/figure1/Spatial/02.Cluster/{id}/{id}_SCT.h5ad'
inputpath = f'../data/figure1/Spatial/01.Segmentation/{id}/'
outpath = f'../data/figure1/Spatial/03.Spatial_cluster/{id}/'
id = sys.argv[4]
e_neigh = 30
s_neigh = 8
resolution = 1

e_neigh = int(e_neigh)
s_neigh = int(s_neigh)
resolution = float(resolution)
data = sc.read(inh5ad)
data.var_names_make_unique()

sc.pp.pca(data,n_comps=30)

## Consider Spatial distance Clustering
def spatially_constrained_cluster(adata, e_neigh=30, s_neigh=6, binsize=None,resolution = 5):
    sc.pp.neighbors(adata, n_neighbors=e_neigh)
    sq.gr.spatial_neighbors(adata, n_neighs=s_neigh)
    conn = adata.obsp['connectivities'].copy()
    conn.data[conn.data > 0] = 1
    adj = conn + adata.obsp['spatial_connectivities']
    adj.data[adj.data  > 0] = 1
    sc.tl.leiden(adata, adjacency=adj, key_added='spatial_leiden_e' + str(e_neigh) + '_s' + str(s_neigh),resolution = resolution)

spatially_constrained_cluster(data,e_neigh=e_neigh, s_neigh=s_neigh, binsize=None,resolution = resolution)

## Spatial plot

a = np.load(inputpath + '/' + str(id) + '_mask.npy')
#a= a[~(a==0).all(1)]

data.uns[str(id)]={}
data.uns[str(id)]['seg_cell']=a
data.obs['Batch'] = str(id)

#name.columns = 'CellID'
data.obs['cell_id'] = data.obs['CellID'].apply(lambda x:int(x.split('.')[1])).values
data.obs['cell_id']= data.obs['cell_id'].astype('category')

data.uns['angle_dict']={}
data.uns['angle_dict'][str(id)]=0

data.write(os.path.join(outpath,str(id)+'.Spatial.distance.cluster.h5ad'))


# RCTD result vis
inh5ad = f'{outpath}/{id}.Spatial.distance.cluster.h5ad'
inputmeta = f'../data/figure1/Spatial/04.RCTD/{id}/{id}.SCT.metadata.csv'
outpath = f'../data/figure1/Spatial/04.RCTD/{id}/'

os.chdir(outpath)
data = sc.read(inh5ad)
data.var_names_make_unique()

Seu = pd.read_csv(inputmeta)
Seu = Seu[['cellid', 'RCTD']]
new_column_names = {'cellid': 'CellID', 'RCTD': 'RCTD'}
Seu = Seu.rename(columns=new_column_names)
obs = pd.merge(data.obs, Seu, left_on='CellID', right_on='CellID', how='left')
obs = obs.set_index('CellID') 
obs.index.name = None


obs.RCTD.value_counts()
data.obs = obs

# remove unknown
cells = obs[obs['RCTD'] != 'unknown'].index
data = data[cells,]

color_df = pd.read_csv('../data/RCTD.color_LastlEdition240612.txt', sep='\t')

featureplot_slices_discrete(obj = data,
    feature = 'RCTD',
    fname = os.path.join(str(id)+'.RCTD_unrm_sep0.Spatial.pdf'),
    show = False,
    scale = True,
    legend_size = 6,
    dpi= 300,
    sep = 0,
    order = color_df['abbreviation'].tolist(),
    colors = color_df['Color'].tolist(),
    slices = None,
    angle_dict = data.uns['angle_dict'],
    nrow = 1,
    ncol = 1,
    compress_factor=False,
    raw=False)



