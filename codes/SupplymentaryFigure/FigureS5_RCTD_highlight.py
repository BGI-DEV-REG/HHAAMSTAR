import scanpy as sc
import matplotlib.pyplot as plt
# import squidpy as sq
import numpy as np
import pandas as pd
import sys
import time
import os
sys.path.append('../')
from Tools.Spatial import featureplot_slices_discrete

inh5ad =  f'../data/figure1/Spatial/03.Spatial_cluster/{id}/{id}.Spatial.distance.cluster.h5ad'
outpath = f'../data/figure1/Spatial/04.RCTD/{id}/'
meta = f'../data/figure1/Spatial/04.RCTD/{id}/{id}.SCT.metadata.csv'

os.chdir(outpath)
data = sc.read(inh5ad)
data.var_names_make_unique()

Seu = pd.read_csv(meta) #, index_col = 0)
Seu = Seu[['cellid', 'celltype']]
new_column_names = {'cellid': 'CellID', 'celltype': 'celltype'}
Seu = Seu.rename(columns=new_column_names)
obs = pd.merge(data.obs, Seu, left_on='CellID', right_on='CellID', how='left')
obs = obs.set_index('CellID')
obs.index.name = None

# remove unknown
cells = obs[obs['celltype'] != 'unknown'].index
data = data[cells,]

color_df = pd.read_csv('../data/RCTD.color_LastlEdition240612.txt', sep='\t')
cluster = color_df.loc[celltype_index, 'abbreviation']
Color = color_df.loc[celltype_index, 'Color']

featureplot_slices_discrete(obj = data,
    feature = 'celltype',
    fname = os.path.join(str(id)+'.'+cluster+'.RCTD.Spatial.highlight_sep0.pdf'),
    show = False,
    scale = True,
    legend_size = 10,
    dpi= 300,
    sep = 0,
    order = [cluster],
    colors = [Color],
    blank_color=(48, 48, 48),
    slices = None,
    nrow = 1,
    ncol = 1,
    compress_factor=False,
    raw=False)


