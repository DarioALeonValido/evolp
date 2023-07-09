import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from pathlib import Path

datapath = Path('../databases_generated/HPA_exp')
output_path = Path('output')

normal = np.loadtxt(datapath / 'normalHPA.dat')
tumor = np.loadtxt(datapath / 'tumorHPA.dat')

ref = normal.mean(axis=0)
data = np.concatenate([normal, tumor]) - ref
U, S, VT = np.linalg.svd(data, full_matrices=0)
data_transform = data @ VT[:3].T

normal = data_transform[:normal.shape[0]]
tumor = data_transform[normal.shape[0]:]

fig_3d = plt.figure(figsize=(6.5, 5))
ax_3d = fig_3d.add_subplot(projection='3d')

ax_3d.scatter(normal[:, 0], normal[:, 1], normal[:, 2], label='Normal', c='b')
ax_3d.scatter(tumor[:, 0], tumor[:, 1], tumor[:, 2], label='Tumor', c='r', marker='s')

ax_3d.set_xlabel('PC1', fontsize=11)
ax_3d.set_ylabel('PC2', fontsize=11)
ax_3d.set_zlabel('PC3', fontsize=11)
ax_3d.tick_params(labelsize=11)
ax_3d.legend(fontsize=11)

fig_3d.tight_layout()
fig_3d.savefig(output_path / 'fig4.pdf')
fig_3d.savefig(output_path / 'fig4.png')

plt.show()
