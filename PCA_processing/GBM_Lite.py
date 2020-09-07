import numpy as np
import xlrd as xl
from scipy.stats.mstats import gmean
from scipy.sparse.linalg import eigsh
from datetime import datetime

def read_data():
    wb = xl.open_workbook("sample.xls")
    ws = wb.sheet_by_index(0)

    directory = '.\\datos6\\'

    sample_type = []
    sanos = []
    tumor = []
    datos = []

    for i in range(1, ws.nrows):
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

        file_object.close()

    datos = np.array(datos, dtype='double') + 0.1

    referencia = gmean(datos[sanos])

    t = np.log2(datos/referencia)

    np.savetxt('results/sanos_dat.dat', sanos, fmt='%i')
    np.savetxt('results/tumor_dat.dat', tumor, fmt='%i')

    return(t)

def CovMat(t, name = "upper_triangular_dat.npy"):
    file = open(name, "wb")
    for i in range(t.T.shape[0]):
        np.save(file, np.dot(t.T[i],t)[i:].astype("f2")/t.shape[0])
    file.close()

def CovMat2(t, m = 5000, name = "upper_triangular_dat.npy"):
    inicio = datetime.today()
    t0 = datetime.today()
    t1 = datetime.today()

    file = open(name, "wb")
    cant = 0
    buffer = []
    for i in range(t.T.shape[0]):
        buffer.append(np.dot(t.T[i],t)[i:].astype("f2")/t.shape[0])
        if(np.mod(i+1, m) == 0):
            t1 = datetime.today()
            print("\ti: ", i, "\t T: ", t1 - inicio, "\tt: ", t1 - t0)
            t0 = t1

            np.save(file, np.concatenate(buffer))
            cant += 1
            buffer = []
    if(buffer != []):
        np.save(file, np.concatenate(buffer))
        cant += 1
    np.savetxt("cant.txt", [cant])
    file.close()

def CovMatCluster(t):
    covariance = np.dot(t.T, t)/t.shape[0]
    return(covariance)

#####################################################
def MultMV(vec, name = "upper_triangular.npy"):
    r = np.zeros_like(vec)
    file = open(name, "rb")
    inicio = datetime.today()
    t0 = datetime.today()
    t1 = datetime.today()
    for i in range(np.alen(vec)):
        if(np.mod(i,100) == 0):
            t1 = datetime.today()
            print("\t\tfila: ", i, "/", np.alen(vec) ,"\t T: ", t1 - inicio, "\tt: ", t1 - t0, "\t <t>: ", (t1-inicio)/(i+1))
            t0 = t1
        buff = np.load(file)
        for j in range(np.alen(buff)):
            r[i] += vec[j + i] * buff[j]
            if(i != i + j):
                r[j + i] += vec[i] * buff[j]
    file.close()
    return(r)

def MultMV2(vec, name = "upper_triangular_dat.npy"):
    r = np.zeros_like(vec)
    file = open(name, "rb")
    c = np.alen(vec)
    i = 0
    j = 0
    while(i < c):
        buff = np.load(file)
        l= np.alen(buff)
        cumulative = 0
        while(cumulative < l):
            r[i] += np.dot(buff[cumulative:cumulative + c - i],vec[i:])
            r[i+1:] += vec[i] * buff[cumulative + 1:cumulative + c - i]
            cumulative += c - i
            i += 1
    file.close()
    return(r)

#####################################################
def comp(vec, w, i):
    for j in range(i + 2):
        buff = np.abs(np.dot(w,vec[j]))
        if(buff > 1e-7):
            return([False, buff])
    return([True, 0])

def ortogonalize(vec, w, i, buff, dim):
    if(buff < 0.05):
        q = w
    else:
        while(buff >= 0.05):
            q = np.random.rand(dim)
            _, buff = comp(vec, q, i)
    # q = w
    for j in range(i + 2):
        q = q  - np.dot(w, vec[j]) * vec[j]

    return(q)

def lanczos(iter, dim):
    #dim = np.alen(m)
    vec = np.zeros((iter + 1, dim))
    vec[1] = np.random.rand(dim)
    vec[1] = vec[1]/np.linalg.norm(vec[1])

    alpha = np.zeros(iter)
    beta = np.zeros(iter)

    inicio = datetime.today()
    t0 = datetime.today()
    t1 = datetime.today()

    for i in range(iter - 1):
        t1 = datetime.today()
        print("\ti: ", i+1, "/", iter ,"\t T: ", t1 - inicio, "\tt: ", t1 - t0)
        t0 = t1
        # w = np.dot(m, vec[i + 1])
        w = MultMV2(vec[i + 1])
        alpha[i] = np.dot(w, vec[i + 1])
        w = w - alpha[i] * vec[i + 1] - beta[i] * vec[i]
        QOrtonormal, valor = comp(vec, w, i)
        if(QOrtonormal == False):
            w = ortogonalize(vec, w, i, valor, dim)
        beta[i + 1] = np.linalg.norm(w)
        vec[i + 2] = w/beta[i + 1]

    # w = np.dot(m, vec[-1])
    w = MultMV2(vec[-1])
    alpha[-1] = np.dot(w, vec[-1])

    return([alpha, beta[1:], vec[1:]])

#####################################################
def calc(t, name = "upper_triangle_dat.npy"):
    #vals, vect = eigsh(covariance, k = 100)

    iter = 100
    a, b, vect_l = lanczos(iter, t.shape[1])

    vals, vect_red = eigsh(np.diag(a) + np.diag(b, -1) + np. diag(b, 1), k = iter)

    vals_norm = vals/vals.sum()

    projection = np.dot(vect.T,t.T)

    return([vals, vect, vals_norm, projection])

def calcCluster(t, covariance):
    vals, vect = eigsh(covariance, k = 100)

    vals_norm = vals/vals.sum()

    projection = np.dot(vect.T,t.T)

    return([vals, vect, vals_norm, projection])


if __name__ == "__main__":
    inicio = datetime.today()
    t0 = 0
    t1 = 0

    print("reading data\t time: ", inicio)
    t = read_data()
    print("constructing covariance matrix\t T: ", datetime.today() - inicio, "\t t: ", datetime.today() - inicio)
    t0 = datetime.today()
    # Comment the next line for working in cluster
    ## CovMat(t)
    #CovMat2(t)

    t1 = datetime.today()
    print("Lanczos\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1

    iter = 101
    a, b, vect_l = lanczos(iter, t.shape[1])

    t1 = datetime.today()
    print("vals, vect_red\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1
    vals, vect_red = eigsh(np.diag(a) + np.diag(b, -1) + np.diag(b, 1), k = iter-1)

    t1 = datetime.today()
    print("Sorting and normalizing\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1

    vals = vals[np.argsort(-np.abs(vals))]
    vect_red = vect_red[np.argsort(-np.abs(vals))]
    vals_norm = vals/vals.sum()

    t1 = datetime.today()
    print("Eighenvectors\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1
    vect = np.dot(vect_red.T, vect_l)
    t1 = datetime.today()
    print("projecting\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1
    projection = np.dot(vect.T,t.T)

    #vals, vect, vals_norm, projection = calc(t)

    # Uncomment the next line for working in cluster
    # covariance = CovMatCluster(t)
    # vals, vect, vals_norm, projection = calc(t, covariance)


    t1 = datetime.today()
    print("saving data\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1
    #las proyecciones van a ir en vectores columna, de forma tal que la primera fila es la proyección de  los ptos en el primer vector
    np.savetxt('results/PC_dat.dat', projection, fmt='%f')
    #los vectores van a estar en vectores columna
    np.savetxt('results/vec_dat.dat', vect, fmt='%f')
    np.savetxt('results/results_dat.dat', vals_norm, fmt='%f')
    np.savetxt('results/eigenvalues_dat.dat', vals, fmt='%f')
    np.savetxt('results/vec2_dat.dat', vect, fmt='%f')

    principal = vect[:, vals.argmax()]
    indices = np.argpartition(-np.abs(principal), 20)[:20]
    componentes = principal[indices]

    np.savetxt('results/20_index_dat.dat', indices, fmt='%i')
    np.savetxt('results/20_components_dat.dat', componentes, fmt='%f')

    t1 = datetime.today()
    print("Ending\t T: ", t1 - inicio, "\tt: ", t1 - t0)
    t0 = t1