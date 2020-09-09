import numpy as np
import xlrd as xl
import PCA_core as pca

# Names ----------------------------------------------------------

datapath="../external_databases/TCGA/"
tissue_id="KIRC"
#possible targets:
#tissue_id=["BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC",
#                  "PRAD","STAD","THCA","UCEC","BLCA","ESCA","READ"]

# Functions ------------------------------------------------------

def read_data(datapath,tissue_id):
    wb = xl.open_workbook(datapath+"samples"+tissue_id+".xls")
    ws = wb.sheet_by_index(0)
    print(ws.nrows,"individual files were identified.")

    sample_type = []
    normal = []
    tumor = []
    data = []

    #reading first file to evaluate the cost
    name = ws.cell_value(1,0)
    sample_type.append(ws.cell_value(1,3))
    file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
    lines = file_object.readlines()
    print("Each file contains",len(lines),"genes.")
    print("Loading large databases...")
    if(ws.cell_value(1,3) == 'Solid Tissue Normal'):
        normal.append(0)
    else:
        tumor.append(0)
    data.append([])
    for line in lines:
        data[-1].append(line.split()[1])        
    file_object.close()

    #reading the rest of the files
    for i in range(2, ws.nrows):        
        name = ws.cell_value(i,0)
        sample_type.append(ws.cell_value(i,3))
        file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
        lines = file_object.readlines()
        if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
            normal.append(i - 1)
        else:
            tumor.append(i - 1)
        data.append([])
        for line in lines:
            data[-1].append(line.split()[1])        
        file_object.close()

    return([np.array(data, dtype='double') + 0.1, normal, tumor]) 
        ## Es realmente necesario usar 'double'? Tampoco se entiende el + 0.1

#def calcc(data, normal):
#    #ref = gmean(data[normal])  ## esto no me funciono asi
#    data_n = []
#    for i in normal:
#        data_n.append(data[i])
#    ref = gmean(data_n)#
#
#    t = np.log2(data/ref)
#    covariance = np.dot(t.T, t)/t.shape[0]
#    #eigenvalues, eigenvectors = eigsh(covariance, k = 100) ## no es necesario para la version lite
#    eigenvalues, eigenvectors = eigh(covariance)
#    eigenvalues_normalized = eigenvalues/eigenvalues.sum()
#    projection = np.dot(eigenvalues.T,t.T)
#    
#    return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])


# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
#if __name__ == "__main__":
data, normal, tumor = read_data(datapath,tissue_id)
print("Data successfully loaded!")
print("normal =", len(normal))
print("tumor =", len(tumor))
print("total =", len(data))

# reducing the data by 30 genes for lite version
print("Lite version: reducing the number of genes by", 30)
data_red = []
for i in range(len(data)):
    data_red.append([])
    for j in range(30):
        data_red[i].append(data[i][j])

print("Computing PCA components...")
#eigenvalues, eigenvectors, eigenvalues_normalized, projection = pca.calc(data, normal) #Full version
eigenvalues, eigenvectors, eigenvalues_normalized, projection = pca.calc(data_red, normal)
principal = eigenvectors[:, eigenvalues.argmax()]
index = np.argpartition(-np.abs(principal), 20)[:20]
components = principal[index]
print("Done!")


# Output data and figures --------------------------------------------------

print("Exporting data...")
np.savetxt('PC_dat.dat', projection.T, fmt='%f')
np.savetxt('eigenvectors_dat.dat', eigenvectors.T, fmt='%f')
np.savetxt('eigenvalues_normalized_dat.dat', eigenvalues_normalized, fmt='%f')
np.savetxt('eigenvalues_dat.dat', eigenvalues, fmt='%f')
np.savetxt('eigenvectorsT_dat.dat', eigenvectors, fmt='%f')
np.savetxt('20_index_dat.dat', index, fmt='%i')
np.savetxt('20_component_dat.dat', components, fmt='%f')
np.savetxt('normal_dat.dat', normal, fmt='%i')
np.savetxt('tumor_dat.dat', tumor, fmt='%i')    
print("Done!")

