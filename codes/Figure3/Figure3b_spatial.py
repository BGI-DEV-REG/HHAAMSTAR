# ============================================================
# Figure 3b - Spatial Plot
# ============================================================

import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import os
os.environ["OPENCV_SHRT_MAX"] = str(pow(2, 40))
import cv2
sys.path.append('/data/work/ST_Analysis/RYTools')
from Tools.Spatial import *

DATA_DIR = "data"
OUT_DIR  = "output"
os.makedirs(OUT_DIR, exist_ok=True)

num = 51

sampleID = pd.read_csv(os.path.join(DATA_DIR, 'Figure3b_sample.csv'), index_col=0)
color_df = pd.DataFrame({
    'abbreviation': ["hHFDSC", "aHFSC", "TAC1", "TAC2", 'Other'],
    'Color': ['#00FF00', '#FFFF00', '#FF34FF', '#1E90FF', '#303030']
})

args_id = sampleID['ID'].iloc[num]
args_sample = sampleID['Sample'].iloc[num]
print(f"Processing sample: {args_id} (num={num})")

inh5ad = "/zfsms3/ST_STOMICS/STOmics_cloud/odms/prd/sz_history/P21Z10200N0171/Project/Hair_follicle/03.Spatial_Cluster/" \
        f"{args_sample}/{args_id}/{args_id}.Spatial.distance.cluster.h5ad"
inputmeta = os.path.join(DATA_DIR, "Figure3b_meta.csv")

data = sc.read(inh5ad)
data.var_names_make_unique()
data.obs['CellID'] = data.obs.index

Seu = pd.read_csv(inputmeta)
Seu = Seu.set_index('cellid')
Seu['cellid'] = Seu.index
Seu = Seu[['cellid', 'Figure3I_cell']]

obs = pd.merge(data.obs, Seu, left_on='CellID', right_index=True, how='left')
obs = obs.set_index('CellID')
obs.index.name = None
print(f"{args_id} cell counts:")
print(obs.Figure3I_cell.value_counts())
data.obs = obs

cells = obs[obs['Figure3I_cell'] != 'unknown'].index
data = data[cells, ]

featureplot_slices_discrete(
    obj=data, feature='Figure3I_cell',
    fname=os.path.join(OUT_DIR, f"Figure3b_{args_id}.cell_unrm.Spatial.pdf"),
    show=False, scale=True, legend_size=6,
    order=color_df['abbreviation'].tolist(),
    colors=color_df['Color'].tolist(),
    slices=None, angle_dict=data.uns['angle_dict'],
    nrow=1, ncol=1, sep=0,
    compress_factor=False, raw=False
)

print("Done!")
