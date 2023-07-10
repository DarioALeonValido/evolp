import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from pathlib import Path
from scipy.stats.mstats import gmean


# Set output path
output_path = Path('output')

# Set path to data
datapath = Path(r'..\databases_generated\TCGA_exp\meanvalues')

# Set tissue types
tissues = np.array(['BLCA', 'BRCA', 'COAD', 'ESCA', 'HNSC', 'KIRC', 'KIRP',
                    'LIHC', 'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD', 'THCA',
                    'UCEC'])

# Load data
data = np.zeros((30, 60483))
for i, tissue in enumerate(tissues):
    data[i] = np.loadtxt(datapath / f'normal{tissue}.dat')
    data[i + 15] = np.loadtxt(datapath / f'tumor{tissue}.dat')

# Calculate reference (normal center)
ref = gmean(data[:15])

# Normalize data
data = np.log2(data/ref)

U, S, VT = np.linalg.svd(data, full_matrices=0)
data_transform = data @ VT[:3].T

normal = data_transform[:15]
tumor = data_transform[15:]
normal[:, 0] = -normal[:, 0]
tumor[:, 0] = -tumor[:, 0]

point1 = tumor.mean(axis=0)
point2 = normal.mean(axis=0)
vec = point1 - point2

d1 = - point1.dot(vec)
d2 = - point2.dot(vec)

pad = 0
xx1, yy1 = np.meshgrid(
    range(
        int(np.floor(tumor[:,0].min()))-pad-30, int(np.ceil(tumor[:,0].max()))+pad+10
        ),
    range(
        int(np.floor(tumor[:,1].min()))-pad, int(np.ceil(tumor[:,1].max()))+pad-10
        ),
    )
xx2, yy2 = np.meshgrid(
    range(
        int(np.floor(normal[:,0].min()))-pad-50, int(np.ceil(normal[:,0].max()))+pad+20
        ),
    range(
        int(np.floor(normal[:,1].min()))-pad + 10, int(np.ceil(normal[:,1].max()))+pad - 90
        ),
    )

x10 = 140
y10 = 160
z10 = 360

point1[2] += 5
zz1 = np.sqrt((point1[0] - x10)**2 + (point1[1] - y10)**2 + (point1[2] + z10)**2\
    - (xx1 - x10)**2 - (yy1 - y10)**2) - z10

x20 = 100
y20 = 130
z20 = 350

point2[2] -= 5
zz2 = np.sqrt((point2[0] - x20)**2 + (point2[1] - y20)**2 + (point2[2] + z20)**2 -\
    (xx2 - x20)**2 - (yy2 - y20)**2) - z20 

zz1[zz1>100] = 100 #point1[2]

# Fig2b
fig_3d = plt.figure(figsize=(6.5, 5))
ax_3d = fig_3d.add_subplot(projection='3d')

ax_3d.scatter(normal[:, 0], normal[:, 1], normal[:, 2], label='Normal', c='b')
ax_3d.scatter(tumor[:, 0], tumor[:, 1], tumor[:, 2], label='Tumor', c='r', marker='s')

ax_3d.set_xlabel('PC1', fontsize=11)
ax_3d.set_ylabel('PC2', fontsize=11)
ax_3d.set_zlabel('PC3', fontsize=11)
ax_3d.tick_params(labelsize=11)
ax_3d.legend(fontsize=11)

ax_3d.set(xlim=(-175, 125))
ax_3d.set(zlim=( -75, 225))
ax_3d.set(zlim=(-125, 105))
ax_3d.view_init(18, -76)

fig_3d.tight_layout()
fig_3d.savefig(output_path / 'fig2a.pdf')
fig_3d.savefig(output_path / 'fig2a.png')

# Scheme
#========
fig_esq = plt.figure(figsize=(6.5, 5))
ax_esq = fig_esq.add_subplot(projection='3d')

ax_esq.plot_surface(xx2, yy2, zz2, alpha=0.6, color='b')
ax_esq.plot_surface(xx1, yy1, zz1, alpha=0.6, color='r')

ax_esq.set(xlim=(-175, 125))
ax_esq.set(zlim=( -75, 225))
ax_esq.set(zlim=(-125, 105))
ax_esq.view_init(18, -76)

# Set x, y, z axes label
ax_esq.set_xlabel('PC1', fontsize=11)
ax_esq.set_ylabel('PC2', fontsize=11)
ax_esq.set_zlabel('PC3', fontsize=11)
legend_elements = [Patch(facecolor='b', alpha=0.6, lw=2, label='Normal'),
                   Patch(facecolor='r', alpha=0.6, lw=2, label='Tumor')]
ax_esq.legend(handles=legend_elements, fontsize=11)
ax_esq.tick_params(labelsize=11)

x2, y2, _ = proj3d.proj_transform(point1[0], point1[1], point1[2], ax_esq.get_proj())
x3, y3, _ = proj3d.proj_transform(point2[0]-20, point2[1]+25, point2[2]-20, ax_esq.get_proj())


label = ax_esq.annotate(
    '', 
    xy=(x2,y2), xytext=(x3, y3),
    arrowprops=dict(arrowstyle='->',
        shrinkA=0, shrinkB=0 ))

angle = 90 + 180/np.pi*np.arctan((y2 - y3)/(x2 - x3))
text_axis = proj3d.proj_transform(-38, 0, -8, ax_esq.get_proj())
ax_esq.annotate('Cancer\nprogression\naxis',
                xy=(0,0), xytext=(text_axis[0], text_axis[1]),
                rotation=angle,
                ha='right',
                va='center',
                fontweight='bold',
                fontsize=10)

text_axis = proj3d.proj_transform(-48, 100, 65, ax_esq.get_proj())
ax_esq.annotate('Cancer\nmanifold',
                xy=(0,0), xytext=(text_axis[0], text_axis[1]),
                fontweight='bold',
                fontsize=10)

text_axis = proj3d.proj_transform(30, 20, -30, ax_esq.get_proj())
ax_esq.annotate('Nomal\ntissue\nmanifold',
                xy=(0,0), xytext=(text_axis[0], text_axis[1]),
                fontweight='bold',
                fontsize=10)

fig_esq.tight_layout()
fig_esq.savefig(output_path / 'fig2b.pdf')
fig_esq.savefig(output_path / 'fig2b.png')

plt.show()
