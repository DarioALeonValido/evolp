"""
This script file computes the entropy and overlap of PCA. 
The procedure and main results are described in paper 
arXiv:2005.02271
For more details on this file see author(s):
FEQ, DALV
"""

import numpy as np
import math, pylab
import xlrd as xl
import matplotlib.pyplot as plt
import sys,os
import random
from entropy import read_PC,compute_entropy

# General variables ------------------------------------------------------

Npc = 20
#tissues_id = ['BRCA','COAD','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','STAD','THCA','UCEC']
t_id='BRCA'
N_random = 1000
Nsamp=10000
N_samples=[20, 30, 50, 100, 200, 300, 400, 500] #cannot exceed Nsamp

srcpc = "../databases_generated/TCGA_pca/"
srcsampl = "../databases_external/TCGA/"
sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'


# Functions ------------------------------------------------

def compute_PCranges(pc_dataset):

    coreM=pc_dataset.T
    cov_mat=np.cov(coreM,ddof=1)
    mean=np.mean(coreM,axis=1)
    stdev=np.std(coreM,axis=1,ddof=1)
    minimum=np.min(coreM,axis=1)
    maximum=np.max(coreM,axis=1)

    return mean,minimum,maximum,stdev,cov_mat


def compute_gaussianS(mean,cov_mat):

    u,ss,vh =np.linalg.svd(cov_mat)
    s=np.diag(ss)
    cov_inv=np.matmul(np.matmul(vh.T,np.linalg.inv(s)),u.T)
    sum_log=sum(math.log(s[i][i]) for i in range(Npc))
    entropy=0.5*(sum_log + Npc*(1 + math.log(2*math.pi)))

    return entropy


def gaussian_ratio(vec0,vec1,cov_mat):

    u,ss,vh =np.linalg.svd(cov_mat)
    s=np.diag(ss)
    cov_inv=np.matmul(np.matmul(vh.T,np.linalg.inv(s)),u.T)
    
    exp0=np.matmul(cov_inv,vec0)
    exp0=np.matmul(vec0.T,exp0)
    exp1=np.matmul(cov_inv,vec1)
    exp1=np.matmul(vec1.T,exp1)

    arg = exp1-exp0
    frac=math.exp(exp1-exp0)

    return frac, arg


def MonteCarlo_simulation(Nsamp,scale,mean,dev,cov_mat):

    data_generated=[]

    dim=len(dev)
    V0=np.random.rand(dim)
    V0=(scale*V0-scale/2)*dev+mean
    Ngen=0
    for i in range(Nsamp):
        V1=np.random.rand(dim)
        V1=(scale*V1-scale/2)*dev + mean
        ratio,arg=gaussian_ratio(V0-mean,V1-mean,cov_mat)
        accept=np.random.rand()
        if(arg>math.log(accept)):
        #if(ratio>accept):
            data_generated.append(V1)
            V0=V1
            Ngen=Ngen+1
        else:
            V0=np.random.rand(dim)
            V0=(scale*V0-scale/2)*dev + mean
            data_generated.append(V0)

    return Ngen/Nsamp, data_generated


def compute_overlap(Npc,mean_n,mean_t,cov_mat_n,cov_mat_t,cov_inv_n,cov_inv_t):

    A, B = Npc/2 , Npc/4
    inv_Vc = cov_inv_n+cov_inv_t
    Vc = np.linalg.inv(inv_Vc)
        
    u_Vc, s_Vc, vh_Vc = np.linalg.svd(Vc)
    u_inv_Vc, s_inv_Vc, vh_inv_Vc = np.linalg.svd(inv_Vc)
    det_Vc, det_inv_Vc = np.prod(s_Vc), np.prod(s_inv_Vc)
        
    u_inv_Vn, s_inv_Vn, vh_inv_Vn = np.linalg.svd(cov_inv_n)
    u_inv_Ve, s_inv_Ve, vh_inv_Ve = np.linalg.svd(cov_inv_t)        
        
    log_sum=sum(-math.log(s_inv_Vn[i])-math.log(s_inv_Ve[i])+math.log(s_inv_Vc[i]) for i in range(Npc))        
        
    eta_n = np.matmul(cov_inv_n,mean_n)
    eta_e = np.matmul(cov_inv_t,mean_t)
    eta_c = eta_n + eta_e
        
    mixed = np.matmul(np.transpose(eta_n),np.matmul(cov_mat_n,eta_n))\
           +np.matmul(np.transpose(eta_e),np.matmul(cov_mat_t,eta_e))\
           -np.matmul(np.transpose(eta_c),np.matmul(Vc,eta_c))\
        
    arg= -0.25*(Npc*math.log(2*math.pi) + log_sum + mixed)        
        
    #overlap = math.pow(2,A)*math.pow(2*math.pi,B)*math.exp(arg)*math.pow(det_inv_Vc,-0.25)
    log_overlap = -A*math.log(2)-B*math.log(2*math.pi)-arg +0.25*math.log(det_inv_Vc)

    return log_overlap


# The main program -----------------------------------------------------------------------------------------
if __name__ == "__main__":
    #random.seed(a=1234567899)
    random.seed()
    print('\nMonte Carlo detailed balance for Gaussian-like entropy\n')

#   Monte Carlo detailed balance --------------------------------
    pc_normal,pc_tumor=read_PC(sample_path,PC_path,t_id,Npc)
    print('normal/tumor:',len(pc_normal),'/',len(pc_tumor))
    print('PCs:',str(len(pc_normal[0]))+'x'+str(len(pc_tumor[0])))

    mean_n,min_n,max_n,stdev_n,cov_mat_n=compute_PCranges(pc_normal)
    mean_t,min_t,max_t,stdev_t,cov_mat_t=compute_PCranges(pc_tumor )
    S_normal=compute_gaussianS(mean_n,cov_mat_n)
    S_tumor =compute_gaussianS(mean_t,cov_mat_t)

    print('asymptotic references Snormal/Stumor:',S_normal,'/',S_tumor)

    np.savetxt('PCranges_dat.dat', np.transpose([min_n,max_n,mean_n,stdev_n,min_t,max_t,mean_t,stdev_t]),
    header='min_n\t\t\t  max_n\t\t  mean_n\t\t\t  sdev_n\t\t  min_t\t\t  max_t\t\t  mean_t\t\t  sdev_t',
    comments='', delimiter="\t  ", fmt="%s")

#   generating data
    print('generating',Nsamp,'normal and tumor samples...')
    scale_n=1.90; scale_t=3.2
    rate_n, normal_generated = MonteCarlo_simulation(Nsamp,scale_n,mean_n,stdev_n,cov_mat_n)
    rate_t, tumor_generated  = MonteCarlo_simulation(Nsamp,scale_t,mean_t,stdev_t,cov_mat_t)

    #normal_generated=normal_generated-np.mean(normal_generated,axis=0)
    #tumor_generated =tumor_generated + 

    np.savetxt('normal_gen_dat.dat', normal_generated,
    header='Generated normal samples: Npc='+str(Npc),
    comments='', delimiter="\t  ", fmt="%s")

    np.savetxt('tumor_gen_dat.dat', tumor_generated,
    header='Generated tumor samples: Npc='+str(Npc),
    comments='', delimiter="\t  ", fmt="%s")

    print('Acceptance rate of norma/tumor samples:',rate_n,'/',rate_t)
#   ----------------------------------------------------------

    S_normal_samp_mean,S_normal_samp_std,S_tumor_samp_mean,S_tumor_samp_std=[],[],[],[]
    lnI_samp_mean,lnI_samp_std=[],[]

    normal_generated = np.array(normal_generated, dtype=np.float64)
    tumor_generated  = np.array(tumor_generated , dtype=np.float64)

    print('Data subsampling with',N_random,'sets of N normal and tumor samples:')

    for N_min in N_samples:
        Sn_i,St_i,lnI_i=[],[],[]
        for i in range(N_random):
            pc_normal_samp=normal_generated[random.sample(range(Nsamp),N_min)]
            pc_tumor_samp=  tumor_generated[random.sample(range(Nsamp),N_min)]
            S_n,mean_n_samp,cov_mat_n_samp,cov_inv_n_samp=compute_entropy(pc_normal_samp)
            S_t,mean_t_samp,cov_mat_t_samp,cov_inv_t_samp=compute_entropy(pc_tumor_samp)
            Sn_i.append(S_n)
            St_i.append(S_t)
            logI_id=compute_overlap(Npc,mean_n_samp,mean_t_samp,cov_mat_n_samp,cov_mat_t_samp,
                                        cov_inv_n_samp,cov_inv_t_samp)
            lnI_i.append(logI_id)
        S_normal_samp_mean.append(np.mean(Sn_i)); S_normal_samp_std.append(np.std(Sn_i,ddof=1))
        S_tumor_samp_mean.append(np.mean(St_i));  S_tumor_samp_std.append(np.std(St_i,ddof=1))
        lnI_samp_mean.append(np.mean(lnI_i));     lnI_samp_std.append(np.std(lnI_i,ddof=1))
        print('N =',N_min)
    print('last Snormal/Stumor:',S_normal_samp_mean[-1],'/',S_tumor_samp_mean[-1])

    print('Data subsampling with',N_random,'sets of N normal samples:')
    lnI_samp_mean_n,lnI_samp_std_n=[],[]
    S_t,mean_t,cov_mat_t,cov_inv_t=compute_entropy(pc_tumor)

    for N_min in N_samples:
        lnI_i=[]
        for i in range(N_random):
            pc_normal_samp=normal_generated[random.sample(range(Nsamp),N_min)]
            S_n,mean_n_samp,cov_mat_n_samp,cov_inv_n_samp=compute_entropy(pc_normal_samp)
            logI_id=compute_overlap(Npc,mean_n_samp,mean_t,cov_mat_n_samp,cov_mat_t,
                                        cov_inv_n_samp,cov_inv_t)
            lnI_i.append(logI_id)
        lnI_samp_mean_n.append(np.mean(lnI_i));lnI_samp_std_n.append(np.std(lnI_i,ddof=1))
        print('N =',N_min)
#   --------------------------------

    print('Plotting figures')
    fig1, ax1 = plt.subplots(3,1,sharex=True)
    fig2, ax2 = plt.subplots()
    f=0.95  
    fig1.set_size_inches(f*1*5, f*3*3)
    fig2.set_size_inches(f*1*5, f*1*3)
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['ytick.labelsize'] = 15

    ax1[0].plot([N_samples[0],N_samples[-1]],[S_normal, S_normal],c='b',ls='--')
    ax1[0].plot([N_samples[0],N_samples[-1]],[S_tumor , S_tumor ],c='r',ls='--')
    ax1[0].plot(N_samples,S_tumor_samp_mean ,c='r', markersize=5, marker='o',label='tumor' )
    ax1[0].plot(N_samples,S_normal_samp_mean,c='b', markersize=5, marker='o',label='normal')

    S_normal_samp_mean=np.array(S_normal_samp_mean);S_normal_samp_std=np.array(S_normal_samp_std)
    S_tumor_samp_mean= np.array(S_tumor_samp_mean );S_tumor_samp_std =np.array(S_tumor_samp_std )

    ax1[1].errorbar(N_samples,S_tumor_samp_mean-S_normal_samp_mean,yerr=S_normal_samp_std+
    S_tumor_samp_std,c='k',ecolor='gray',ls='-',marker='o',fmt='o',capsize=3)
    ax1[2].errorbar(N_samples,S_tumor-S_normal_samp_mean,yerr=S_normal_samp_std,
    c='k',ecolor='gray',ls='-',marker='o',fmt='o',capsize=3)

    ax1[0].legend(edgecolor='k', borderaxespad=0.2)
    ax1[0].set_ylabel('Entropy'), ax1[1].set_ylabel(r'$\Delta S$') 
    ax1[2].set_xlabel(r'$N_{samples}$'), ax1[2].set_ylabel(r'$S^{\infty}_{tumor}-S^N_{normal}$')
    fig1.tight_layout()
    fig1.savefig(t_id+'_MC_sim_entropy_fig.pdf')

    ax2.errorbar(N_samples[1:],lnI_samp_mean[1:],yerr=lnI_samp_std[1:],label=r'$N_t=N_n$',
    c='orange',ls='-',marker='o',fmt='o',capsize=3)
    ax2.errorbar(N_samples[1:],lnI_samp_mean_n[1:],yerr=lnI_samp_std_n[1:],label=r'$N_t=\infty$',
    c='green',ls='-',marker='o',fmt='o',capsize=3)
    
    ax2.set_ylim(8,18)
    ax2.set_xlabel(r'$N_{samples}$'), ax2.set_ylabel(r'$-ln(I)$')
    ax2.legend(edgecolor='k', borderaxespad=0.2)
    fig2.tight_layout()
    fig2.savefig(t_id+'_MC_sim_overlap_fig.pdf')

print("Done!\n")
print("Check out the outputs!\n")
