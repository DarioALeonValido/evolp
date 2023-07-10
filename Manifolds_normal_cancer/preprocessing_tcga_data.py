import numpy as np
import pandas as pd

from pathlib import Path
from scipy.stats.mstats import gmean
from tqdm import tqdm

tissues = np.array(['BLCA', 'BRCA', 'COAD', 'ESCA', 'HNSC', 'KIRC', 'KIRP',
                    'LIHC', 'LUAD', 'LUSC', 'PRAD', 'READ', 'STAD', 'THCA',
                    'UCEC'])
databases_path = Path(r'..\databases_external\TCGA')
output_path = Path('output')

for tissue in tissues:
    data_path = databases_path / f'data{tissue}' 

    sample = pd.read_excel(databases_path / f'sample{tissue}.xls')
    normal = sample['Sample Type'] == 'Solid Tissue Normal'
    tumor = sample['Sample Type'] != 'Solid Tissue Normal'

    fpkm = pd.read_table(data_path / sample['File Name'][0],
                         names=['gene', 'value'])
    data = np.zeros([sample.shape[0], fpkm.shape[0]])
    for i, filename in enumerate(tqdm(sample['File Name'])):
        data[i] = pd.read_table(data_path / filename,
                                names=['gene', 'value'])['value']

    np.savetxt(output_path/f'normal{tissue}.dat', gmean(data[normal] + 0.1))
    np.savetxt(output_path/f'tumor{tissue}.dat', gmean(data[tumor] + 0.1))
