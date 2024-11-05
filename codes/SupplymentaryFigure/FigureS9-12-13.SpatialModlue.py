import scanpy as sc
### cluster viss sc
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import pandas as pd
import sys
import time
import os
os.environ["OPENCV_SHRT_MAX"] = str(pow(2,40))
import cv2
import dynamo as dyn
sys.path.append('../')
from Tools.Spatial import *

sample = id.split('_SCT_selection')[0]
inh5ad = f'../data/figure1/Spatial/03.Spatial_cluster/{sample}/{id}/{id}.Spatial.distance.cluster.h5ad'
outpath = f'/../data/figure4/GO/Spatial_MarkerModule/{sample}/{id}/'

if not os.path.exists(outpath):
    os.makedirs(outpath,exist_ok=True)

os.chdir(outpath)

data = sc.read(inh5ad)
data.var_names_make_unique()

meta = pd.read_csv(meta)
meta = meta.set_index('cellid')
meta

data.obs[feature] = meta[feature]
data.obs

featureplot_slices_continuous(obj=data,
                            cmap='viridis', #调色板
                            feature=feature,
                            fname=feature+'.Spatial.pdf',
                            sep = 0,
                            show=False,
                            title=f"{feature}",
                            scale=True,
                            nrow=1,
                            ncol=1,
                            compress_factor=2,
                            raw=False)
