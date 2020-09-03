import numpy as np
import xlrd as xl
from scipy.stats.mstats import gmean
from scipy.sparse.linalg import eigsh

wb = xl.open_workbook("/nfs/home/icimaf06/3. TCGA-KIRC/sample.xls")
#ws = wb["gdc_sample_sheet.2019-04-16"]
ws = wb.sheet_by_index(0)

directory = '/nfs/home/icimaf06/3. TCGA-KIRC/datos3/'

sample_type = []
#data = []
sanos = []
tumor = []
datos = []
#pos_sano = []
#pos_tumor = []

for i in range(1, ws.nrows):
    
#    data.append([])
    nombre = ws.cell_value(i,0)
    sample_type.append(ws.cell_value(i,3))

    file_object = open(file = directory + nombre, mode = "r")
    lines = file_object.readlines()

    if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
        sanos.append(i - 1)
    else:
        tumor.append(i - 1)

    datos.append([])
    for line in lines:
        datos[-1].append(line.split()[1])


datos = np.array(datos, dtype='double') + 0.1

referencia = gmean(datos[sanos])

t = np.log2(datos/referencia)
#t_reference = np.mean(t, axis=0)

#a = t - t_reference

covariance = np.dot(t.T, t)/t.shape[0]

vals, vect = eigsh(covariance, k = 100)

vals_normalizados = vals/vals.sum()

proyeccion = np.dot(vect.T,t.T)


np.savetxt('results2/PC', proyeccion.T, fmt='%f')
np.savetxt('results2/vec', vect.T, fmt='%f')
np.savetxt('results2/resultados', vals_normalizados, fmt='%f')
np.savetxt('results2/autovalores', vals, fmt='%f')
np.savetxt('results2/vec2', vect, fmt='%f')

principal = vect[:, vals.argmax()]
indices = np.argpartition(-np.abs(principal), 20)[:20]
componentes = principal[indices]

np.savetxt('results2/20_indices.txt', indices, fmt='%i')
np.savetxt('results2/20_componentes.txt', componentes, fmt='%f')
np.savetxt('results2/sanos.txt', sanos, fmt='%i')
np.savetxt('results2/tumor.txt', tumor, fmt='%i')