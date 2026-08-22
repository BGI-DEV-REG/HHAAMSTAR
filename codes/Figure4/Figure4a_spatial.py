# ============================================================
# Figure 4a - Spatial Plot
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

sampleID = pd.read_csv(os.path.join(DATA_DIR, 'Figure4a_sample.csv'), index_col=0)
color_df = pd.read_csv(os.path.join(DATA_DIR, 'RCTD.color1_251210.csv'))

for num in range(0, 77):
    try:
        args_id = sampleID['ID'].iloc[num]
        args_sample = sampleID['Sample'].iloc[num]
        print(f"Processing sample: {args_id} (num={num})")

        inh5ad = "/zfsms3/ST_STOMICS/STOmics_cloud/odms/prd/sz_history/P21Z10200N0171/Project/Hair_follicle/03.Spatial_Cluster/" \
                f"{args_sample}/{args_id}/{args_id}.Spatial.distance.cluster.h5ad"
        inputmeta = os.path.join(DATA_DIR, "Figure4a_meta.csv")

        data = sc.read(inh5ad)
        data.var_names_make_unique()

        Seu = pd.read_csv(inputmeta)
        Seu = Seu[Seu['Sample'] == args_id].copy()
        Seu = Seu[['cellid', 'RCTD']]

        obs = pd.merge(data.obs, Seu, left_on='CellID', right_on='cellid', how='left')
        obs = obs.set_index('CellID')
        obs.index.name = None
        print(f"{args_id} RCTD counts:")
        print(obs.RCTD.value_counts())
        data.obs = obs

        cells = obs[obs['RCTD'] != 'unknown'].index
        data = data[cells, ]

        featureplot_slices_discrete(
            obj=data, feature='RCTD',
            fname=os.path.join(OUT_DIR, f"Figure4a_{args_id}.RCTD_unrm.Spatial.pdf"),
            show=False, scale=True, legend_size=6,
            order=color_df['abbreviation'].tolist(),
            colors=color_df['Color'].tolist(),
            slices=None, angle_dict=data.uns['angle_dict'],
            nrow=1, ncol=1, sep=0,
            compress_factor=False, raw=False
        )
    except Exception as e:
        print(f"Error processing num={num}: {str(e)}")
        continue

print("All samples done!")
