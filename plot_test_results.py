#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################### RUN IMPORTS ############################
# from asyncio import coroutines
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib
from matplotlib import rcParams
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 10
import os, sys
import matplotlib.patches as patches
from matplotlib.ticker import NullFormatter, LinearLocator
from scipy import stats
import json
from scipy import stats
import pandas
from numpy import unravel_index
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

####################### SET PARAMS ############################################
# input directory, contains all the outputs from the MCMC code

Layer_Hist      = True
Layers_Aniso_Tr = True
Posterior       = True
Posterior2      = True
Sigmad          = True
Dispersion      = True
Convergence     = True
PsPp_Fit        = True
Corr_Hist      = True



# python plot_test_results.py ./OUT_TEST2

num_args=len(sys.argv)

if Corr_Hist:
    if num_args < 8:
        print('Required :       Script.py <dir> <av_int> <dis> <p1> <av_u_dep1> <p2> <av_u_dep2>')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     RJ_MCMC_Tests/XYZ_test/OUT_TEST')
        print('Options [2] :     Averaging interval, 50,100,200')
        print('Options [3] :     hist discretisation, 50,100,200')
        print('Options [4] :     Parameter 1, vph, vsv')
        print('Options [5] :     Parameter 1, av_u_dep: 100, 200 ....')
        print('Options [6] :     Parameter 2, vsv, xi')
        print('Options [7] :     Parameter 2, av_u_dep: 100, 200 ....')
        print('Number of arguments (' + str(num_args) +') too low... exit')
        exit('exiting....')
    
    av_int = int(sys.argv[2])
    if av_int!=50 and av_int!=100 and av_int!=200:
        print('Bad av_int supplied...')
        print('Options [2] :     Averaging interval, 50,100,200')
        exit('exiting....')

    dh=int(sys.argv[3])
    if dh!=50 and dh!=100 and dh!=200:
        print('Bad histogram discretization supplied...')
        print('Options [3] :     hist discretisation, 50,100,200')
        exit('exiting....')

    p1=str(sys.argv[4])
    if p1!='vph' and p1!='vsv':
        print('Bad Parameter 1 supplied...')
        print('Options [4] :     Parameter 1, vph, vsv')
        exit('exiting....')
    p1_u_dep=int(sys.argv[5])

    p2=str(sys.argv[6])
    if p2!='vsv' and p2!='xi':
        print('Bad Parameter 2 supplied...')
        print('Options [6] :     Parameter 2, vsv, xi')
        exit('exiting....')
    p2_u_dep=int(sys.argv[7])



else:
    if num_args < 2:
        print('Required :       Script.py test_directory')
        print('Arguments:      ',sys.argv)
        print('Options [1] :     RJ_MCMC_Tests/XYZ_test/OUT_TEST')
        print('Number of arguments (' + str(num_args) +') too low... exit')
        exit('exiting....')

directory = str(sys.argv[1])
print('Plotting results for: ' +str(directory))

if not os.path.exists(directory+'/PLOTS/'):
    os.mkdir(directory+'/PLOTS/')
maxnlay = 80

fname_pre=str(os.getcwd().split('/')[-1])+'_'


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average, math.sqrt(variance))



#########################################################
# Get Input models to inversion.
#########################################################


def get_inversion_model(filename='RM_ak135.txt'):
    '''
    reads the input inverse models into a dictionary.
    contains: 
        npt_true: number of points
        depth: depth of interfaces
        vsv, xi, vph

    Parameters
    ----------
    filename : str
        file containing the models saved from run_POC_synth.sh.
        E.g., RM_ak135.txt RM_craton.txt RM_craton_meta.txt
    
    Returns
    -------
    model : dict 
        contains the data of the inversion model.

    '''
    model           = {}

    file=filename
    # read file
    print(file)

    f=open(file,'r')
    lines=f.readlines()
    f.close()

    npt_true=len(lines)
    # depths AK135Vp AK135Vs Vp dVp VsH dVsH VsV dVsV  Vs_PR06 dVs_PR06 Xi_PR06 dXi_PR06 Vs_FW10 dVs_FW10 Xi_FW10 dXi_FW10 (calculate w.r.t. nul value, PR06=1, FW10=0)

    d=np.zeros(npt_true)
    xi=np.ones(npt_true)
    vsv=np.zeros(npt_true)
    vph=np.zeros(npt_true)


    for i in range(npt_true):
        data=lines[i].split()
        d[i]=float(data[0])
        xi[i]=float(data[11])
        vsv[i]=float(data[7])*1000.0
        vph[i]=float(data[3])*1000.0

    model['npt_true']=npt_true
    model['depth']=d
    model['vsv']=vsv
    model['xi']=xi
    model['vph']=vph

    return model




########## histogram of number of layers
if Layer_Hist:
    if os.path.isfile(directory+'/'+'NB_layers.out'):
            
        f=open(directory+'/'+'NB_layers.out','r')
        lines=f.readlines()
        f.close()

        n=[]
        for line in lines:
            n.append(float(line))
        
        plt.figure('nlayers_histogram')
        plt.title('nlayers histogram')
        plt.plot(n)
        plt.xlim([0, maxnlay])
        plt.xlabel('Number of layers')
        plt.ylabel('Frequency')
        
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Layer_hist.png',dpi=200)
        plt.close()

########## matrix containing number of layers and number of anisotropic layers
if Layers_Aniso_Tr:
    if os.path.isfile(directory+'/'+'Tradeoff.out'):
        f=open(directory+'/'+'Tradeoff.out','r')
        lines=f.readlines()
        f.close()

        nlay=int(lines[0])
        l=np.zeros((nlay+1,nlay))

        i=0
        j=0
        for line in lines[1:]:
            l[j,i]=float(line.split()[0])#/(i+1)
            j+=1
            if j==nlay+1:
                i+=1
                j=0

        plt.figure('Layers_Aniso_Tr')
        plt.title('Number of Layers versus Anisotropic Layers')
        plt.xlim([0, nlay])
        plt.ylim([0, nlay])
        plt.xlabel('Number of Layers')
        plt.ylabel('Number of Anisotropic Layers')
        plt.viridis()
        plt.pcolor(l)
        plt.colorbar()
        plt.plot([0, nlay], [0, nlay], '--r', linewidth=1)
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Layers_Aniso_Tr.png',dpi=200)
        plt.close()

############ Posterior #######################
if Posterior:

    if os.path.isfile(directory+'/'+'Proc_Posterior.out'):
        # Use the Output of Process_MCMC_output.py from transcale.
        file=open(directory+'/'+'Proc_Posterior.out','r')
    else:
        # Use with Original Posterior output from Fortran
        file=open(directory+'/'+'Posterior.out','r')

    lines=file.readlines()

    file.close()

    nlines = len(lines)
    data0=lines[0].split()
    prof=float(data0[0])
    disd=int(data0[1]) # Depth discretisation
    dmax=float(data0[2])

    burn_in=int(float(data0[3]))
    nsample=int(float(data0[4]))
    thinning=int(float(data0[5]))
    cores=int(float(data0[6]))


    [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[1].split()]
    disv=int(disv) # Velocity/aniso discretisation

    vsvd=np.zeros((disd,disv))
    xid=np.zeros((disd,disv))
    vpd=np.zeros((disd,disv))

    vsvs=np.linspace(vref_min,vref_max,disv)
    xis=np.linspace(xi_min,xi_max,disv)
    vps=np.linspace(vp_min,vp_max,disv)
    depths=np.linspace(0,prof,disd)

    i=0
    j=0

    for line in lines[2:]:
        [vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]
        j+=1
        if j==disv:
            j=0
            i+=1

    for i in range(np.shape(vsvd)[0]): # normalisation by depth, not strictly necessary
        s=np.amax(vsvd[i,:])
        vsvd[i,:]=vsvd[i,:]/s
    vsvd[-1,:]=0.0

    for i in range(np.shape(vpd)[0]): # normalisation by depth, not strictly necessary
        s=np.amax(vpd[i,:])
        vpd[i,:]=vpd[i,:]/s
    vpd[-1,:]=0.0
    
    for i in range(np.shape(xid)[0]): # normalisation by depth, not strictly necessary
        s=np.amax(xid[i,:])
        xid[i,:]=xid[i,:]/s
    xid[-1,:]=0.0

    ################# Average distribution ABOVE depth

    dep_val=200
    ind=np.where(depths<dep_val)[0]
    splice_vsvd_a=vsvd[ind,:]
    splice_xid_a=xid[ind,:]
    splice_vpd_a=vpd[ind,:]

    above_splice_vsvd=np.sum(splice_vsvd_a, axis=0)/len(ind) # Normalise by num layers above
    above_splice_xid=np.sum(splice_xid_a, axis=0)/len(ind)
    above_splice_vpd=np.sum(splice_vpd_a, axis=0)/len(ind)

    ################# Average distribution BELOW depth

    ind2=np.where(depths>dep_val)[0]
    splice_vsvd_b=vsvd[ind2,:]
    splice_xid_b=xid[ind2,:]
    splice_vpd_b=vpd[ind2,:]

    below_splice_vsvd=np.sum(splice_vsvd_b, axis=0)/len(ind2) # Normalise by num layers below
    below_splice_xid=np.sum(splice_xid_b, axis=0)/len(ind2) 
    below_splice_vpd=np.sum(splice_vpd_b, axis=0)/len(ind2)

    ################# Average distribution AT depth

    ind3=np.where(abs(depths-dep_val)==np.min(np.abs(depths-dep_val)))[0][0]
    splice_vsvd=vsvd[ind3,:]
    splice_xid=xid[ind3,:]
    splice_vpd=vpd[ind3,:]

    ### STATS ###

    dist_vsvs=[]
    dist_vsvs_b=[]
    dist_vsvs_a=[]
    for i in range(len(vsvs)):
        for j in range(int(1000*splice_vsvd[i])):
            dist_vsvs.append(vsvs[i])
        for k in range(int(1000*below_splice_vsvd[i])):
            dist_vsvs_b.append(vsvs[i])
        for l in range(int(1000*above_splice_vsvd[i])):
            dist_vsvs_a.append(vsvs[i])



    dist_xis=[]
    dist_xis_a=[]
    dist_xis_b=[]
    for i in range(len(xis)):
        for j in range(int(1000*splice_xid[i])):
            dist_xis.append(xis[i])
        for k in range(int(1000*above_splice_xid[i])):
            dist_xis_a.append(xis[i])
        for l in range(int(1000*below_splice_xid[i])):
            dist_xis_b.append(xis[i])



    dist_vps=[]
    dist_vps_a=[]
    dist_vps_b=[]
    for i in range(len(vps)):
        for j in range(int(1000*splice_vpd[i])):
            dist_vps.append(vps[i])
        for k in range(int(1000*below_splice_vpd[i])):
            dist_vps_b.append(vps[i])
        for l in range(int(1000*above_splice_vpd[i])):
            dist_vps_a.append(vps[i])

    std_vsvs=np.round_(np.std(dist_vsvs,axis=0),4)
    mean_vsvs=np.round_(np.mean(dist_vsvs,axis=0),4)
    median_vsvs=np.round_(np.median(dist_vsvs),4)
    vsvs_975=np.round_(np.quantile(dist_vsvs, 0.975),4)
    vsvs_825=np.round_(np.quantile(dist_vsvs, 0.825),4)
    vsvs_175=np.round_(np.quantile(dist_vsvs, 0.175),4)
    vsvs_025=np.round_(np.quantile(dist_vsvs, 0.025),4)

    std_vsvs_a=np.round_(np.std(dist_vsvs_a,axis=0),4)
    mean_vsvs_a=np.round_(np.mean(dist_vsvs_a,axis=0),4)
    median_vsvs_a=np.round_(np.median(dist_vsvs_a),4)
    vsvs_975_a=np.round_(np.quantile(dist_vsvs_a, 0.975),4)
    vsvs_825_a=np.round_(np.quantile(dist_vsvs_a, 0.825),4)
    vsvs_175_a=np.round_(np.quantile(dist_vsvs_a, 0.175),4)
    vsvs_025_a=np.round_(np.quantile(dist_vsvs_a, 0.025),4)

    std_vsvs_b=np.round_(np.std(dist_vsvs_b,axis=0),4)
    mean_vsvs_b=np.round_(np.mean(dist_vsvs_b,axis=0),4)
    median_vsvs_b=np.round_(np.median(dist_vsvs_b),4)
    vsvs_975_b=np.round_(np.quantile(dist_vsvs_b, 0.975),4)
    vsvs_825_b=np.round_(np.quantile(dist_vsvs_b, 0.825),4)
    vsvs_175_b=np.round_(np.quantile(dist_vsvs_b, 0.175),4)
    vsvs_025_b=np.round_(np.quantile(dist_vsvs_b, 0.025),4)

    std_xis=np.round_(np.std(dist_xis,axis=0),4)
    mean_xis=np.round_(np.mean(dist_xis,axis=0),4) # Mean includes zeros in xi
    median_xis=np.round_(np.median(dist_xis),4)
    xis_975=np.round_(np.quantile(dist_xis, 0.975),4)
    xis_825=np.round_(np.quantile(dist_xis, 0.825),4)
    xis_175=np.round_(np.quantile(dist_xis, 0.175),4)
    xis_025=np.round_(np.quantile(dist_xis, 0.025),4)

    std_xis_a=np.round_(np.std(dist_xis_a,axis=0),4) # Xi=1 removed later.
    mean_xis_a=np.round_(np.mean(dist_xis_a,axis=0),4)
    median_xis_a=np.round_(np.median(dist_xis_a),4)
    xis_975_a=np.round_(np.quantile(dist_xis_a, 0.975),4)
    xis_825_a=np.round_(np.quantile(dist_xis_a, 0.825),4)
    xis_175_a=np.round_(np.quantile(dist_xis_a, 0.175),4)
    xis_025_a=np.round_(np.quantile(dist_xis_a, 0.025),4)

    std_xis_b=np.round_(np.std(dist_xis_b,axis=0),4)
    mean_xis_b=np.round_(np.mean(dist_xis_b,axis=0),4)
    median_xis_b=np.round_(np.median(dist_xis_b),4)
    xis_975_b=np.round_(np.quantile(dist_xis_b, 0.975),4)
    xis_825_b=np.round_(np.quantile(dist_xis_b, 0.825),4)
    xis_175_b=np.round_(np.quantile(dist_xis_b, 0.175),4)
    xis_025_b=np.round_(np.quantile(dist_xis_b, 0.025),4)

    std_vps=np.round_(np.std(dist_vps,axis=0),4)
    mean_vps=np.round_(np.mean(dist_vps,axis=0),4)
    median_vps=np.round_(np.median(dist_vps),4)
    vps_975=np.round_(np.quantile(dist_vps, 0.975),4)
    vps_825=np.round_(np.quantile(dist_vps, 0.825),4)
    vps_175=np.round_(np.quantile(dist_vps, 0.175),4)
    vps_025=np.round_(np.quantile(dist_vps, 0.025),4)

    std_vps_b=np.round_(np.std(dist_vps_b,axis=0),4)
    mean_vps_b=np.round_(np.mean(dist_vps_b,axis=0),4)
    median_vps_b=np.round_(np.median(dist_vps_b),4)
    vps_975_b=np.round_(np.quantile(dist_vps_b, 0.975),4)
    vps_825_b=np.round_(np.quantile(dist_vps_b, 0.825),4)
    vps_175_b=np.round_(np.quantile(dist_vps_b, 0.175),4)
    vps_025_b=np.round_(np.quantile(dist_vps_b, 0.025),4)

    std_vps_a=np.round_(np.std(dist_vps_a,axis=0),4)
    mean_vps_a=np.round_(np.mean(dist_vps_a,axis=0),4)
    median_vps_a=np.round_(np.median(dist_vps_a),4)
    vps_975_a=np.round_(np.quantile(dist_vps_a, 0.975),4)
    vps_825_a=np.round_(np.quantile(dist_vps_a, 0.825),4)
    vps_175_a=np.round_(np.quantile(dist_vps_a, 0.175),4)
    vps_025_a=np.round_(np.quantile(dist_vps_a, 0.025),4)


    # print(mean_vsvs_a,std_vsvs_a,mean_vsvs,std_vsvs,mean_vsvs_b,std_vsvs_b)
    # print(mean_xis_a,std_xis_a,mean_xis,std_xis,mean_xis_b,std_xis_b)
    # print(mean_vps_a,std_vps_a,mean_vps,std_vps,mean_vps_b,std_vps_b)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, figsize=(9, 9), sharey=True) #

    if os.path.isfile(directory+'/'+'Proc_Posterior.out'):
        # When using Proc_Post the isotropic lies in the next index over
        xid[:,disv//2]=xid[:,disv//2-1] # else the isotropic layers dominate the anisotropy density plot
        above_splice_xid[disv//2]=above_splice_xid[disv//2-1]
        splice_xid[disv//2]=splice_xid[disv//2-1]
        below_splice_xid[disv//2]=below_splice_xid[disv//2-1]

        above_splice_xid=above_splice_xid/np.max(above_splice_xid)
        below_splice_xid=below_splice_xid/np.max(below_splice_xid)
        splice_xid=splice_xid/np.max(splice_xid)

    else:
        # xid[:,disv//2-1]=0 # else the isotropic layers dominate the anisotropy density plot
        xid[:,disv//2-1]=xid[:,disv//2] # else the isotropic layers dominate the anisotropy density plot
        above_splice_xid[disv//2-1]=above_splice_xid[disv//2]
        splice_xid[disv//2-1]=splice_xid[disv//2]
        below_splice_xid[disv//2-1]=below_splice_xid[disv//2]

    # print(np.max(above_splice_xid),np.max(splice_xid),np.max(below_splice_xid))

    ax1.plot(vsvs,above_splice_vsvd, c='r', linewidth=1,alpha=1)
    ax2.plot(xis,above_splice_xid, c='b', linewidth=1,alpha=1)
    ax3.plot(vps,above_splice_vpd, c='darkgreen', linewidth=1,alpha=1)

    ax4.plot(vsvs,splice_vsvd, c='r', linewidth=1,alpha=1)
    ax5.plot(xis,splice_xid, c='b', linewidth=1,alpha=1)
    ax6.plot(vps,splice_vpd, c='darkgreen', linewidth=1,alpha=1)

    ax7.plot(vsvs,below_splice_vsvd, c='r', linewidth=1,alpha=1)
    ax8.plot(xis,below_splice_xid, c='b', linewidth=1,alpha=1)
    ax9.plot(vps,below_splice_vpd, c='darkgreen', linewidth=1,alpha=1)

    ax1.fill_between(vsvs,above_splice_vsvd,0, where=above_splice_vsvd >= 0,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vsvs_025_a<=vsvs;     b=vsvs_975_a>=vsvs;     c=a==b
    ax1.fill_between(vsvs,above_splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vsvs_175_a<=vsvs;     b=vsvs_825_a>=vsvs;     c=a==b
    ax1.fill_between(vsvs,above_splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.2)
    #########################################################
    ax4.fill_between(vsvs,splice_vsvd,0, where=splice_vsvd >= 0,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vsvs_025<=vsvs;     b=vsvs_975>=vsvs;     c=a==b
    ax4.fill_between(vsvs,splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vsvs_175<=vsvs;     b=vsvs_825>=vsvs;     c=a==b
    ax4.fill_between(vsvs,splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.2)
    ########################################################
    ax7.fill_between(vsvs,below_splice_vsvd,0, where=below_splice_vsvd >= 0,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vsvs_025_b<=vsvs;     b=vsvs_975_b>=vsvs;     c=a==b
    ax7.fill_between(vsvs,below_splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vsvs_175_b<=vsvs;     b=vsvs_825_b>=vsvs;     c=a==b
    ax7.fill_between(vsvs,below_splice_vsvd,0, where=c,facecolor=[1.0, 0., 0.0], rasterized=False, alpha=0.2)


    #########################################################
    ax2.fill_between(xis,above_splice_xid,0, where=above_splice_xid >= 0,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=xis_025_a<=xis;     b=xis_975_a>=xis;     c=a==b
    ax2.fill_between(xis,above_splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=xis_175_a<=xis;     b=xis_825_a>=xis;     c=a==b
    ax2.fill_between(xis,above_splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.2)

    #########################################################
    ax5.fill_between(xis,splice_xid,0, where=splice_xid >= 0,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=xis_025<=xis;     b=xis_975>=xis;     c=a==b
    ax5.fill_between(xis,splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=xis_175<=xis;     b=xis_825>=xis;     c=a==b
    ax5.fill_between(xis,splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.2)

    #########################################################
    ax8.fill_between(xis,below_splice_xid,0, where=below_splice_xid >= 0,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.1)
    # For 95% quantile
    a=xis_025_b<=xis;     b=xis_975_b>=xis;     c=a==b
    ax8.fill_between(xis,below_splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.15)
    # For 65% quantile
    a=xis_175_b<=xis;     b=xis_825_b>=xis;     c=a==b
    ax8.fill_between(xis,below_splice_xid,0, where=c,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.2)

    #########################################################
    ax3.fill_between(vps,above_splice_vpd,0, where=above_splice_vpd >= 0,facecolor='darkgreen', rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vps_025_a<=vps;     b=vps_975_a>=vps;     c=a==b
    ax3.fill_between(vps,above_splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vps_175_a<=vps;     b=vps_825_a>=vps;     c=a==b
    ax3.fill_between(vps,above_splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.2)

    #########################################################
    ax6.fill_between(vps,splice_vpd,0, where=splice_vpd >= 0,facecolor='darkgreen', rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vps_025<=vps;     b=vps_975>=vps;     c=a==b
    ax6.fill_between(vps,splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vps_175<=vps;     b=vps_825>=vps;     c=a==b
    ax6.fill_between(vps,splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.2)

    #########################################################
    ax9.fill_between(vps,below_splice_vpd,0, where=below_splice_vpd >= 0,facecolor='darkgreen', rasterized=False, alpha=0.1)
    # For 95% quantile
    a=vps_025_b<=vps;     b=vps_975_b>=vps;     c=a==b
    ax9.fill_between(vps,below_splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.15)
    # For 65% quantile
    a=vps_175_b<=vps;     b=vps_825_b>=vps;     c=a==b
    ax9.fill_between(vps,below_splice_vpd,0, where=c,facecolor='darkgreen', rasterized=False, alpha=0.2)


    ax1.plot([mean_vsvs_a,mean_vsvs_a],[0,1],c='red',linewidth=1,alpha=1)
    ax2.plot([mean_xis_a,mean_xis_a],[0,1],c='red',linewidth=1,alpha=1)
    ax3.plot([mean_vps_a,mean_vps_a],[0,1],c='red',linewidth=1,alpha=1)

    ax1.plot([median_vsvs_a,median_vsvs_a],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax2.plot([median_xis_a,median_xis_a],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax3.plot([median_vps_a,median_vps_a],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)

    ax4.plot([mean_vsvs,mean_vsvs],[0,1],c='red',linewidth=1,alpha=1)
    ax5.plot([mean_xis,mean_xis],[0,1],c='red',linewidth=1,alpha=1)
    ax6.plot([mean_vps,mean_vps],[0,1],c='red',linewidth=1,alpha=1)

    ax4.plot([median_vsvs,median_vsvs],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax5.plot([median_xis,median_xis],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax6.plot([median_vps,median_vps],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)

    ax7.plot([mean_vsvs_b,mean_vsvs_b],[0,1],c='red',linewidth=1,alpha=1)
    ax8.plot([mean_xis_b,mean_xis_b],[0,1],c='red',linewidth=1,alpha=1)
    ax9.plot([mean_vps_b,mean_vps_b],[0,1],c='red',linewidth=1,alpha=1)

    ax7.plot([median_vsvs_b,median_vsvs_b],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax8.plot([median_xis_b,median_xis_b],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)
    ax9.plot([median_vps_b,median_vps_b],[0,1],c='red', linestyle='--',linewidth=1,alpha=1)

    ax1.annotate('(a)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax2.annotate('(b)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax3.annotate('(c)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))

    ax4.annotate('(d)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax5.annotate('(e)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax6.annotate('(f)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))

    ax7.annotate('(g)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax8.annotate('(h)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax9.annotate('(i)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))


    ax1.annotate('median={:.2f}'.format(median_vsvs_a),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax1.annotate('mean={:.2f}'.format(mean_vsvs_a),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax2.annotate('median={:.2f}'.format(median_xis_a),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax2.annotate('mean={:.2f}'.format(mean_xis_a),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax3.annotate('median={:.2f}'.format(median_vps_a),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax3.annotate('mean={:.2f}'.format(mean_vps_a),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax4.annotate('median={:.2f}'.format(median_vsvs),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax4.annotate('mean={:.2f}'.format(mean_vsvs),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax5.annotate('median={:.2f}'.format(median_xis),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax5.annotate('mean={:.2f}'.format(mean_xis),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax6.annotate('median={:.2f}'.format(median_vps),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax6.annotate('mean={:.2f}'.format(mean_vps),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax7.annotate('median={:.2f}'.format(median_vsvs_b),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax7.annotate('mean={:.2f}'.format(mean_vsvs_b),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax8.annotate('median={:.2f}'.format(median_xis_b),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax8.annotate('mean={:.2f}'.format(mean_xis_b),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

    ax9.annotate('median={:.2f}'.format(median_vps_b),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
    ax9.annotate('mean={:.2f}'.format(mean_vps_b),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=8,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')


    ax1.set_title('$V_{SV}$ above '+str(dep_val)+'km')
    ax1.set_xlabel('$V_{SV}$ (m/s)', fontsize=10)
    ax1.set_ylabel('Norm. Probability',fontsize=10)
    ax1.set_xlim([vref_min,vref_max])
    ax1.set_ylim([0.0,1.0])

    ax2.set_title('Rad Anis. above '+str(dep_val)+'km')
    ax2.set_xlabel('$Xi$',fontsize=10)
    ax2.set_xlim([xi_min,xi_max])

    ax3.set_title('$V_{PH}$ above '+str(dep_val)+'km')
    ax3.set_xlim([vp_min,vp_max])
    ax3.set_xlabel('$V_{PH}$ (m/s)',fontsize=10)

    ax4.set_title('$V_{SV}$ at '+str(dep_val)+'km')
    ax4.set_xlabel('$V_{SV}$ (m/s)', fontsize=10)
    ax4.set_ylabel('Norm. Probability',fontsize=10)
    ax4.set_xlim([vref_min,vref_max])
    ax4.set_ylim([0.0,1.0])

    ax5.set_title('Rad Anis. at '+str(dep_val)+'km')
    ax5.set_xlabel('$Xi$',fontsize=10)
    ax5.set_xlim([xi_min,xi_max])

    ax6.set_title('$V_{PH}$ at '+str(dep_val)+'km')
    ax6.set_xlim([vp_min,vp_max])
    ax6.set_xlabel('$V_{PH}$ (m/s)',fontsize=10)

    ax7.set_title('$V_{SV}$ below '+str(dep_val)+'km')
    ax7.set_xlabel('$V_{SV}$ (m/s)', fontsize=10)
    ax7.set_ylabel('Norm. Probability',fontsize=10)
    ax7.set_xlim([vref_min,vref_max])
    ax7.set_ylim([0.0,1.0])

    ax8.set_title('Rad Anis. below '+str(dep_val)+'km')
    ax8.set_xlabel('$Xi$',fontsize=10)
    ax8.set_xlim([xi_min,xi_max])

    ax9.set_title('$V_{PH}$ below '+str(dep_val)+'km')
    ax9.set_xlim([vp_min,vp_max])
    ax9.set_xlabel('$V_{PH}$ (m/s)',fontsize=10)

    plt.subplots_adjust(wspace=0.35,hspace=0.35)

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Hist_dist_'+str(dep_val)+'.png',dpi=200)
    plt.close()

    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))

    vsvs=np.linspace(vref_min,vref_max,disv+1)
    xis=np.linspace(xi_min,xi_max,disv+1)
    vps=np.linspace(vp_min,vp_max,disv+1)
    depths=np.linspace(0,prof,disd+1)

    ax0.invert_yaxis()
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('$V_{SV}$ (m/s)', fontsize=10)
    ax0.set_title('S-wave velocity')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel('$Xi$',fontsize=10)
    ax1.set_xlim([xi_min,xi_max])
    # ax1.set_xlim([0.5, 1.5])
    ax1.set_title('Radial Anisotropy')
    ax2.set_xlim([vp_min,vp_max])
    # ax2.set_xlim([-0.5,0.5])
    ax0.set_ylabel('Depth (km)',fontsize=10)
    # ax2.set_xlabel(r'VpVs*(1+X)',fontsize=10)
    # ax2.set_title('Vp/Vs deviation')
    ax2.set_xlabel('$V_{PH}$ (m/s)',fontsize=10)
    ax2.set_title('P-wave velocity')

    ax1.pcolormesh(xis,depths,xid,cmap='viridis')
    ax0.pcolormesh(vsvs,depths,vsvd,cmap='viridis')

    # true model overlaid on the posterior (only for synthetic tests)
    if os.path.isfile(directory+'/'+'true_model.out'):
        file=open(directory+'/'+'true_model.out','r')
        lines=file.readlines()
        file.close()

        true_depth=[]
        true_vsv=[]
        true_xi=[]
        true_vp=[]
        for line in lines[1:]:
            data=line.split()
            true_depth.append(float(data[0]))
            true_vsv.append(float(data[1]))
            try:
                true_xi.append(float(data[2]))
                true_vp.append(float(data[3]))
                # pass
            except:
                pass

        # True models in cyan.
        true,=ax0.plot(true_vsv,true_depth,c='white',linewidth=1)
        try:
            ax1.plot(true_xi,true_depth,c='white',linewidth=1,alpha=1,marker='o',markersize=2,mfc='k')
            ax2.plot(true_vp,true_depth,c='white',linewidth=1)

        except:
            pass

    ax2.pcolormesh(vps,depths,vpd,cmap='viridis')
    plt.setp(ax2.get_yticklabels(), visible=False)

    if os.path.isfile(directory+'/'+'Change_points.out'):
        # histogram of depths for layer boundaries
        file=open(directory+'/'+'Change_points.out','r')
        lines=file.readlines()
        file.close()

        change_depths=[]
        change_hist=[]
        for line in lines:
            data=line.split()
            change_depths.append(float(data[0]))
            change_hist.append(float(data[1]))

        ax5.plot(change_hist,change_depths)
    ax5.set_title('Change Points')
    ax5.set_xlabel('Frequency',fontsize=10)

    # average model overlaid on the posterior (only for synthetic tests)
    # and anisotropy probability

    depths=[]
    average_vs=[]
    average_xi=[]
    average_vp=[]
    average_probani=[]

    if os.path.isfile(directory+'/'+'Average.out'):

        file=open(directory+'/'+'Average.out','r')
        lines=file.readlines()
        file.close()

        for line in lines:
            data=line.split()
            depths.append(float(data[0]))
            average_vs.append(float(data[1]))
            average_xi.append(float(data[2]))
            average_vp.append(float(data[3]))
            average_probani.append(float(data[4]))
        
    # Average models in red.
    ave,=ax0.plot(average_vs,depths,c='r',linewidth=1)
    ax1.plot(average_xi,depths,c='r',linewidth=1)
    ax2.plot(average_vp,depths,c='r',linewidth=1)
    ax4.plot(average_probani,depths,c='k',linewidth=1)

    ax4.set_xlabel('Probability',fontsize=10)
    ax4.set_title('Anisotropy')
    ax4.set_xlim([0,100])


    # Try plotting LSQR CPinSeis results too...
    
    # model_1=get_inversion_model(filename='/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_fund/RM_craton.txt')
    # model_2=get_inversion_model(filename='/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_fund/RM_craton_meta.txt')
    # model_3=get_inversion_model(filename='/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_over/RM_craton.txt')
    # model_4=get_inversion_model(filename='/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_over/RM_craton_meta.txt')


    # ax0.plot(model_1['vsv'],model_1['depth'],c='magenta',linewidth=1)
    # ax1.plot(model_1['xi'],model_1['depth'],c='magenta',linewidth=1)
    # ax2.plot(model_1['vph'],model_1['depth'],c='magenta',linewidth=1)






    if os.path.isfile(directory+'/'+'true_model.out'):
        plt.legend([true, ave], ['true', 'ave'])
    else:
        plt.legend([ave], ['ave'])
    fig.suptitle('Posterior and Averages')

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Posterior.png',dpi=200)
    plt.close()


if Posterior2:


    if os.path.isfile(directory+'/'+'Proc_Posterior.out'):
        # Use the Output of Process_MCMC_output.py from transcale.
        file=open(directory+'/'+'Proc_Posterior.out','r')
    else:
        # Use with Original Posterior output from Fortran
        file=open(directory+'/'+'Posterior.out','r')

    lines=file.readlines()

    file.close()

    nlines = len(lines)
    data0=lines[0].split()
    prof=float(data0[0])
    disd=int(data0[1]) # Depth discretisation
    dmax=float(data0[2])

    burn_in=int(float(data0[3]))
    nsample=int(float(data0[4]))
    thinning=int(float(data0[5]))
    cores=int(float(data0[6]))


    [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[1].split()]
    disv=int(disv) # Velocity/aniso discretisation

    vsvd=np.zeros((disd,disv))
    xid=np.zeros((disd,disv))
    vpd=np.zeros((disd,disv))

    vsvs=np.linspace(vref_min,vref_max,disv)
    xis=np.linspace(xi_min,xi_max,disv)
    vps=np.linspace(vp_min,vp_max,disv)
    depths=np.linspace(0,prof,disd)

    i=0
    j=0

    for line in lines[2:]:
        [vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]
        j+=1
        if j==disv:
            j=0
            i+=1

    vsvd_means=[]; vsvd_medians=[]; vsvd_975=[]; vsvd_825=[]; vsvd_175=[]; vsvd_025=[]; vsvd_std=[]
    xid_means=[]; xid_medians=[]; xid_975=[]; xid_825=[]; xid_175=[]; xid_025=[]; xid_std=[]
    vpd_means=[]; vpd_medians=[]; vpd_975=[]; vpd_825=[]; vpd_175=[]; vpd_025=[]; vpd_std=[]

    x_vals=[]; y_vals=[]; vsvd_vals=[]; xid_vals=[]; vpd_vals=[]
    for i in range(np.shape(vsvd)[0]):
        for j in range(np.shape(vsvd)[1]):
            vsvd_vals.append(int(vsvd[i,j]))
            xid_vals.append(int(vsvd[i,j]))
            vpd_vals.append(int(vsvd[i,j]))
            y_vals.append(depths[i])
            x_vals.append(vsvs[j])

    hist_vsvd, bins_y, bins_x=np.histogram2d(y_vals, x_vals, weights=vsvd_vals, bins=(len(depths),len(vsvs)))
    hist_xid,  bins_y, bins_x=np.histogram2d(y_vals, x_vals, weights=xid_vals,  bins=(len(depths),len(xis)))
    hist_vpd,  bins_y, bins_x=np.histogram2d(y_vals, x_vals, weights=vpd_vals,  bins=(len(depths),len(vps)))



    for i in range(np.shape(hist_vsvd)[0]):
        print('layer: '+str(i))
        temp_vsvd=[]
        temp_xid =[]
        temp_vpd =[]
        for j in range(np.shape(hist_vsvd)[1]):
            val_vsvd=int(vsvd[i,j])
            for k in range(0,val_vsvd):
                temp_vsvd.append(vsvs[j])

            val_xid=int(xid[i,j])
            for k in range(0,val_xid):
                temp_xid.append(xis[j])

            val_vpd=int(vpd[i,j])
            for k in range(0,val_vpd):
                temp_vpd.append(vps[j])

        # vsvd_means.append(np.round_(np.mean(temp_vsvd),4))
        vsvd_medians.append(np.round_(np.median(temp_vsvd),4))
        vsvd_975.append(np.round_(np.quantile(temp_vsvd, 0.975),4))
        vsvd_825.append(np.round_(np.quantile(temp_vsvd, 0.825),4))
        vsvd_175.append(np.round_(np.quantile(temp_vsvd, 0.175),4))
        vsvd_025.append(np.round_(np.quantile(temp_vsvd, 0.025),4))
        vsvd_std.append(np.round_(np.std(temp_vsvd),4))

        # xid_means.append(np.round_(np.mean(temp_xid),4))
        xid_medians.append(np.round_(np.median(temp_xid),4))
        xid_975.append(np.round_(np.quantile(temp_xid, 0.975),4))
        xid_825.append(np.round_(np.quantile(temp_xid, 0.825),4))
        xid_175.append(np.round_(np.quantile(temp_xid, 0.175),4))
        xid_025.append(np.round_(np.quantile(temp_xid, 0.025),4))
        xid_std.append(np.round_(np.std(temp_xid),4))

        # vpd_means.append(np.round_(np.mean(temp_vpd),4))
        vpd_medians.append(np.round_(np.median(temp_vpd),4))
        vpd_975.append(np.round_(np.quantile(temp_vpd, 0.975),4))
        vpd_825.append(np.round_(np.quantile(temp_vpd, 0.825),4))
        vpd_175.append(np.round_(np.quantile(temp_vpd, 0.175),4))
        vpd_025.append(np.round_(np.quantile(temp_vpd, 0.025),4))
        vpd_std.append(np.round_(np.std(temp_vpd),4))

    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))

    # vsvs=np.linspace(vref_min,vref_max,disv+1)
    # xis=np.linspace(xi_min,xi_max,disv+1)
    # vps=np.linspace(vp_min,vp_max,disv+1)
    # depths=np.linspace(0,prof,disd+1)

    ax0.invert_yaxis()
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('$V_{SV}$ (m/s)', fontsize=10)
    ax0.set_title('S-wave velocity')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel('$Xi$',fontsize=10)
    ax1.set_xlim([xi_min,xi_max])
    ax1.set_title('Radial Anisotropy')
    ax2.set_xlim([vp_min,vp_max])
    ax0.set_ylabel('Depth (km)',fontsize=10)
    ax2.set_xlabel('$V_{PH}$ (m/s)',fontsize=10)
    ax2.set_title('P-wave velocity')

    # ax1.pcolormesh(xis,depths,xid,cmap='viridis')
    # ax0.pcolormesh(vsvs,depths,vsvd,cmap='viridis')

    ax0.fill_betweenx(depths,vsvd_025,vsvd_975,facecolor=[1.0, 0.0, 0.0], rasterized=False, alpha=0.15)
    ax0.fill_betweenx(depths,vsvd_175,vsvd_825,facecolor=[1.0, 0.0, 0.0], rasterized=False, alpha=0.25)

    ax1.fill_betweenx(depths,xid_025,xid_975,facecolor=[0.0, 0.0, 1.0], rasterized=False, alpha=0.15)
    ax1.fill_betweenx(depths,xid_175,xid_825,facecolor=[0.0, 0.0, 1.0], rasterized=False, alpha=0.25)

    ax2.fill_betweenx(depths,vpd_025,vpd_975,facecolor='darkgreen', rasterized=False, alpha=0.15)
    ax2.fill_betweenx(depths,vpd_175,vpd_825,facecolor='darkgreen', rasterized=False, alpha=0.25)

    # mean2,=ax0.plot(vsvd_means,depths,c='b',linewidth=1)
    # ax1.plot(xid_means,depths,c='b',linewidth=1)
    # ax2.plot(vpd_means,depths,c='b',linewidth=1)

    median2,=ax0.plot(vsvd_medians,depths,c='b',linewidth=2)
    ax1.plot(xid_medians,depths,c='b',linewidth=2)
    ax2.plot(vpd_medians,depths,c='b',linewidth=2)

    # true model overlaid on the posterior (only for synthetic tests)
    if os.path.isfile(directory+'/'+'true_model.out'):
        file=open(directory+'/'+'true_model.out','r')
        lines=file.readlines()
        file.close()

        true_depth=[]
        true_vsv=[]
        true_xi=[]
        true_vp=[]
        for line in lines[1:]:
            data=line.split()
            true_depth.append(float(data[0]))
            true_vsv.append(float(data[1]))
            try:
                true_xi.append(float(data[2]))
                true_vp.append(float(data[3]))
                # pass
            except:
                pass

        # True models in cyan.
        true,=ax0.plot(true_vsv,true_depth,c='k',linewidth=1)
        try:
            ax1.plot(true_xi,true_depth,c='k',linewidth=1) # ,alpha=1,marker='o',markersize=2,mfc='k'
            ax2.plot(true_vp,true_depth,c='k',linewidth=1)
        except:
            pass

    # ax2.pcolormesh(vps,depths,vpd,cmap='viridis')
    plt.setp(ax2.get_yticklabels(), visible=False)

    change_depths=[]
    change_hist=[]
    if os.path.isfile(directory+'/'+'Change_points.out'):

        # histogram of depths for layer boundaries
        file=open(directory+'/'+'Change_points.out','r')
        lines=file.readlines()
        file.close()


        for line in lines:
            data=line.split()
            change_depths.append(float(data[0]))
            change_hist.append(float(data[1]))

    ax5.plot(change_hist,change_depths)
    ax5.set_title('Change Points')
    ax5.set_xlabel('Frequency',fontsize=10)

    # average model overlaid on the posterior (only for synthetic tests)
    # and anisotropy probability
    depths=[]
    average_vs=[]
    average_xi=[]
    average_vp=[]
    average_probani=[]

    if os.path.isfile(directory+'/'+'Average.out'):
        file=open(directory+'/'+'Average.out','r')
        lines=file.readlines()
        file.close()

        for line in lines:
            data=line.split()
            depths.append(float(data[0]))
            average_vs.append(float(data[1]))
            average_xi.append(float(data[2]))
            average_vp.append(float(data[3]))
            average_probani.append(float(data[4]))
    
    # Average models in red.
    mean,=ax0.plot(average_vs,depths,c='r',linewidth=2)
    ax1.plot(average_xi,depths,c='r',linewidth=2)
    ax2.plot(average_vp,depths,c='r',linewidth=2)
    ax4.plot(average_probani,depths,c='k',linewidth=1)
    ax4.set_xlabel('Probability',fontsize=10)
    ax4.set_title('Anisotropy')

    ax4.set_xlim([0,100])
    # Make proxy artists to make the legend work
    q95,=ax0.fill(np.NaN, np.NaN, 'r', alpha=0.15)
    q65,=ax0.fill(np.NaN, np.NaN, 'r', alpha=0.25)
    if os.path.isfile(directory+'/'+'true_model.out'):
        plt.legend([true, mean, median2, q95, (q95, q65)], ['true', 'mean', 'median', '95% mods.', '65% mods.'],loc='lower right')
    else:
        plt.legend([mean, median2, q95, (q95, q65)], ['mean', 'median', '95% mods.', '65% mods.'],loc='lower right')

    ax0.annotate('(a)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax1.annotate('(b)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax2.annotate('(c)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax4.annotate('(d)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
    ax5.annotate('(e)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=10,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))

    fig.suptitle('Posterior and Averages')

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Posterior2.png',dpi=300)
    plt.close()


#################### histogram of rayleigh uncertainty parameter

if Sigmad:
    if os.path.isfile(directory+'/'+'Sigmad_R.out'):
        file=open(directory+'/'+'Sigmad_R.out','r')
        lines=file.readlines()
        file.close()

        data_init=lines[0].split()
        ad_r_min=float(data_init[0])
        ad_r_max=float(data_init[1])
        disa=int(float(data_init[2]))

        d=[]
        sigmad_R=[]
        for line in lines[1:]:
            data=line.split()
            d.append(float(data[0]))
            sigmad_R.append(float(data[1]))

        plt.figure('sigmad_R')
        plt.title('sigmad_R')
        plt.xlim([ad_r_min, ad_r_max])
        plt.xlabel('Rayleigh uncertainty parameter')
        plt.ylabel('Frequency')
        plt.plot(d,sigmad_R)
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Sigmad_R.png',dpi=200)
        # histogram of love uncertainty parameter
        file=open(directory+'/'+'Sigmad_L.out','r')
        lines=file.readlines()
        file.close()

        data_init=lines[0].split()
        ad_l_min=float(data_init[0])
        ad_l_max=float(data_init[1])
        disa=int(float(data_init[2]))

        d=[]
        sigmad_L=[]
        for line in lines[1:]:
            data=line.split()
            d.append(float(data[0]))
            sigmad_L.append(float(data[1]))

        plt.figure('sigmad_L')
        plt.title('sigmad_L')
        plt.xlim([ad_l_min, ad_l_max])
        plt.xlabel('Love uncertainty parameter')
        plt.ylabel('Frequency')
        plt.plot(d,sigmad_L)
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Sigmad_L.png',dpi=200)
        plt.close()

    if os.path.isfile(directory+'/'+'Sigmad_PsPp.out'):
        
        # histogram of PsPp uncertainty parameter
        file=open(directory+'/'+'Sigmad_PsPp.out','r')
        lines=file.readlines()
        file.close()

        data_init=lines[0].split()
        ad_PsPp_min=float(data_init[0])
        ad_PsPp_max=float(data_init[1])
        disa=int(float(data_init[2]))

        d=[]
        sigmad_PsPp=[]
        for line in lines[1:]:
            data=line.split()
            d.append(float(data[0]))
            sigmad_PsPp.append(float(data[1]))

        plt.figure('sigmad_PsPp')
        plt.title('sigmad_PsPp')
        plt.xlim([ad_PsPp_min, ad_PsPp_max])
        plt.xlabel('PsPp uncertainty parameter')
        plt.ylabel('Frequency')
        plt.plot(d,sigmad_PsPp)
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Sigmad_PsPp.png',dpi=200)
        plt.close()


#################################### CONVERGENCE #################################

if Convergence:
    if os.path.isfile(directory+'/'+'Convergence_misfit.out'):

        ############### misfit change over time, for one core and average over all cores #################
        file=open(directory+'/'+'Convergence_misfit.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        conv_R_one=[]
        conv_R=[]
        conv_L_one=[]
        conv_L=[]
        conv_PsPp_one=[]
        conv_PsPp=[]
        for line in lines[1:]:
            data=line.split()
            conv_R_one.append(float(data[0]))
            conv_R.append(float(data[1]))
            conv_L_one.append(float(data[2]))
            conv_L.append(float(data[3]))
            try:
                conv_PsPp_one.append(float(data[4]))
                conv_PsPp.append(float(data[5]))
            except:
                pass        
                


        plt.figure('convergence_misfit')
        plt.title('Misfit Convergence')

        plt.plot(conv_R_one[burn_in:],label='Rayleigh, one core')
        plt.plot(conv_R[burn_in:],label='Rayleigh, all cores')
        plt.plot(conv_L_one[burn_in:],label='Love, one core')
        plt.plot(conv_L[burn_in:],label='Love, all cores')
        try:
            plt.plot(conv_PsPp_one[burn_in:],label='PsPp, one core')
            plt.plot(conv_PsPp[burn_in:],label='PsPp, all cores')
        except:
            pass
        plt.xlim([0,nsample])
        plt.xlabel('Iteration number')
        plt.ylabel('Misfit')
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_misfit.png',dpi=200)
        plt.close()
        
        ############### number of layers over time, for one core and average over all cores ###############
        file=open(directory+'/'+'Convergence_nb_layers.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        conv_n_one=[]
        conv_n=[]
        for line in lines[1:]:
            data=line.split()
            conv_n_one.append(float(data[0]))
            conv_n.append(float(data[1]))

        plt.figure('convergence_nlayers')
        plt.title('Number of Layers Convergence')
        plt.plot(conv_n_one[burn_in:],label='nblayers, one core')
        plt.plot(conv_n[burn_in:],label='nblayers, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('Number of Layers')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_layers.png',dpi=200)
        plt.close()


        ################ rayleigh uncertainty parameter over time, for one core and average over all cores ###############
        file=open(directory+'/'+'Convergence_sigma_R.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        conv_sigmaR_one=[]
        conv_sigmaR=[]
        for line in lines[1:]:
            data=line.split()
            conv_sigmaR_one.append(float(data[0]))
            conv_sigmaR.append(float(data[1]))

        plt.figure('convergence_sigmaR')
        plt.title('sigmaR convergence')
        plt.plot(conv_sigmaR_one[burn_in:],label='sigmaR, one core')
        plt.plot(conv_sigmaR[burn_in:],label='sigmaR, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('sigmaR')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_sigma_R.png',dpi=200)
        plt.close()

        ################ love uncertainty parameter over time, for one core and average over all cores ###############
        file=open(directory+'/'+'Convergence_sigma_L.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        conv_sigmaL_one=[]
        conv_sigmaL=[]
        for line in lines[1:]:
            data=line.split()
            conv_sigmaL_one.append(float(data[0]))
            conv_sigmaL.append(float(data[1]))

        plt.figure('convergence_sigmaL')
        plt.title('sigmaL convergence')
        plt.plot(conv_sigmaL_one[burn_in:],label='sigmaL, one core')
        plt.plot(conv_sigmaL[burn_in:],label='sigmaL, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('sigmaL')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_sigma_L.png',dpi=200)
        plt.close()

    ################ acceptance rates for birth of isotropic and anisotropic layers over time, for one core and average over all cores ###############
    if os.path.isfile(directory+'/'+'Convergence_Birth.out'):
        file=open(directory+'/'+'Convergence_Birth.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        convB_one=[]
        convB=[]
        convBa_one=[]
        convBa=[]
        for line in lines[1:]:
            data=line.split()
            convB_one.append(float(data[0]))
            convB.append(float(data[1]))
            convBa_one.append(float(data[2]))
            convBa.append(float(data[3]))

        plt.figure('convergence_birth')
        plt.title('Birth Rate Convergence')
        plt.plot(convB_one[burn_in:],label='birth, one core')
        plt.plot(convB[burn_in:],label='birth, all cores')
        plt.plot(convBa_one[burn_in:],label='birth anisotropic, one core')
        plt.plot(convBa[burn_in:],label='birth anisotropic, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_Birth.png',dpi=200)
        plt.close()

    ################ acceptance rates for death of isotropic and anisotropic layers over time, for one core and average over all cores ###############
    if os.path.isfile(directory+'/'+'Convergence_Death.out'):
        file=open(directory+'/'+'Convergence_Death.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        convD_one=[]
        convD=[]
        convDa_one=[]
        convDa=[]
        for line in lines[1:]:
            data=line.split()
            convD_one.append(float(data[0]))
            convD.append(float(data[1]))
            convDa_one.append(float(data[2]))
            convDa.append(float(data[3]))

        plt.figure('convergence_death')
        plt.title('Death Rate Convergence')
        plt.plot(convD_one[burn_in:],label='death, one core')
        plt.plot(convD[burn_in:],label='death, all cores')
        plt.plot(convDa_one[burn_in:],label='death anisotropic, one core')
        plt.plot(convDa[burn_in:],label='death anisotropic, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_Death.png',dpi=200)
        plt.close()

    ################ ####################################################### ###############
    if os.path.isfile(directory+'/'+'Convergence_xi.out'):
        file=open(directory+'/'+'Convergence_xi.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        convP_one=[]
        convP=[]
        for line in lines[1:]:
            data=line.split()
            convP_one.append(float(data[0]))
            convP.append(float(data[1]))

        plt.figure('convergence_xi')
        plt.title('Anisotropy Change Rate Convergence')
        plt.plot(convP_one[burn_in:],label='xi change, one core')
        plt.plot(convP[burn_in:],label='xi change, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_xi.png',dpi=200)
        plt.close()

    ################ #######################################################s ###############
    if os.path.isfile(directory+'/'+'Convergence_vp.out'):

        file=open(directory+'/'+'Convergence_vp.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        convP_one=[]
        convP=[]
        for line in lines[1:]:
            data=line.split()
            convP_one.append(float(data[0]))
            convP.append(float(data[1]))

        plt.figure('convergence_vp')
        plt.title('Vp Change Rate Convergence')
        plt.plot(convP_one[burn_in:],label='vp change, one core')
        plt.plot(convP[burn_in:],label='vp change, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_vp.png',dpi=200)
        plt.close()


    
    ################ ####################################################### ###############
    if os.path.isfile(directory+'/'+'Convergence_vs.out'):
        file=open(directory+'/'+'Convergence_vs.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        convsv1_one=[]
        convsv1=[]
        for line in lines[1:]:
            data=line.split()
            convsv1_one.append(float(data[0]))
            convsv1.append(float(data[1]))

        plt.figure('convergence_vs')
        plt.title('Vsv Change Rate Convergence, upper half')
        plt.plot(convsv1_one[burn_in:],label='vsv change upper half, one core')
        plt.plot(convsv1[burn_in:],label='vsv change upper half, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_vs.png',dpi=200)
        plt.close()

    ############################################################################
    if os.path.isfile(directory+'/'+'Convergence_dp.out'):

        file=open(directory+'/'+'Convergence_dp.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        condp1_one=[]
        condp1=[]
        for line in lines[1:]:
            data=line.split()
            condp1_one.append(float(data[0]))
            condp1.append(float(data[1]))

        plt.figure('convergence_dp')
        plt.title('Depth Change Rate Convergence, upper half')
        plt.plot(condp1_one[burn_in:],label='depth change upper half, one core')
        plt.plot(condp1[burn_in:],label='depth change upper half, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('XXX')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_dp.png',dpi=200)
        plt.close()




    if os.path.isfile(directory+'/'+'Convergence_sigma_PsPp.out'):

        ################ love uncertainty parameter over time, for one core and average over all cores ###############
        file=open(directory+'/'+'Convergence_sigma_PsPp.out','r')
        lines=file.readlines()
        file.close()

        burn_in=int(lines[0].split()[0])
        nsample=int(lines[0].split()[1])

        conv_sigmaPsPp_one=[]
        conv_sigmaPsPp=[]
        for line in lines[1:]:
            data=line.split()
            conv_sigmaPsPp_one.append(float(data[0]))
            conv_sigmaPsPp.append(float(data[1]))

        plt.figure('convergence_sigmaPsPp')
        plt.title('sigmaPsPp convergence')
        plt.plot(conv_sigmaPsPp_one[burn_in:],label='sigmaPsPp, one core')
        plt.plot(conv_sigmaPsPp[burn_in:],label='sigmaPsPp, all cores')
        plt.xlabel('Iteration number')
        plt.ylabel('sigma_PsPp')
        plt.xlim([0,nsample])
        plt.legend()
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Convergence_sigma_PsPp.png',dpi=200)
        plt.close()


################################# Average dispersion curves ################################################
if Dispersion:
    if os.path.isfile(directory+'/'+'Dispersion_mean.out'):
        file=open(directory+'/'+'Dispersion_mean.out','r')
        lines=file.readlines()
        file.close()

        ndatad_R=int(lines[0].split()[0])
        ndatad_L=int(lines[0].split()[1])

        period_R=[]
        n_R=[]
        c_R=[]
        dc_R=[]
        for line in lines[1:ndatad_R+1]:
            data=line.split()
            period_R.append(float(data[0]))
            n_R.append(int(float(data[1])))
            c_R.append(float(data[2]))
            dc_R.append(float(data[3]))

        period_L=[]
        n_L=[]
        c_L=[]
        dc_L=[]
        for line in lines[ndatad_R+1:]:
            data=line.split()
            period_L.append(float(data[0]))
            n_L.append(int(float(data[1])))
            c_L.append(float(data[2]))
            dc_L.append(float(data[3]))

        # true dispersion curves (data)
        file=open(directory+'/'+'Dispersion_obs.out','r')
        lines=file.readlines()
        file.close()

        ndatad_R=int(lines[0].split()[0])
        ndatad_L=int(lines[0].split()[1])

        period_R_obs=[]
        n_R_obs=[]
        c_R_obs=[]
        dc_R_obs=[]
        for line in lines[1:ndatad_R+1]:
            data=line.split()
            period_R_obs.append(float(data[0]))
            n_R_obs.append(int(float(data[1])))
            c_R_obs.append(float(data[2]))
            dc_R_obs.append(float(data[3]))

        period_L_obs=[]
        n_L_obs=[]
        c_L_obs=[]
        dc_L_obs=[]
        for line in lines[ndatad_R+1:]:
            data=line.split()
            period_L_obs.append(float(data[0]))
            n_L_obs.append(int(float(data[1])))
            c_L_obs.append(float(data[2]))
            dc_L_obs.append(float(data[3]))

        Modes_R=np.unique(n_R)
        Modes_L=np.unique(n_L)

        plt.figure('dispersion')

        if len(Modes_R)>0:
            for R_mode in Modes_R:
                ave_R='ave_R_'+str(R_mode)
                obs_R='obs_R_'+str(R_mode)
                # print(ave_R)
                ind=np.where(n_R==R_mode)
                ave_R=plt.errorbar(np.array(period_R)[ind[0]],np.array(c_R)[ind[0]],yerr=np.array(dc_R)[ind[0]],marker='o',zorder=0,label='Rayleigh average',mfc='blue',mec='blue', c='blue')
                obs_R=plt.errorbar(np.array(period_R_obs)[ind[0]]-0.1,np.array(c_R_obs)[ind[0]],yerr=np.array(dc_R_obs)[ind[0]],marker='o',zorder=0,label='Rayleigh observed',mfc='limegreen',mec='limegreen', c='limegreen')

        if len(Modes_L)>0:
            for L_mode in Modes_L:
                ave_L='ave_L_'+str(L_mode)
                obs_L='obs_L_'+str(L_mode)
                # print(ave_L)
                ind=np.where(n_L==L_mode)
                ave_L=plt.errorbar(np.array(period_L)[ind[0]]+0.1,np.array(c_L)[ind[0]],yerr=np.array(dc_L)[ind[0]],marker='o',zorder=0,label='Love average',mfc='orange',mec='orange', c='orange')
                obs_L=plt.errorbar(np.array(period_L_obs)[ind[0]]+0.2,np.array(c_L_obs)[ind[0]],yerr=np.array(dc_L_obs)[ind[0]],marker='o',zorder=0,label='Love observed',mfc='red',mec='red', c='red')


        # plt.errorbar(np.array(period_R),c_R,yerr=dc_R,marker='o',zorder=0,label='Rayleigh average')
        # plt.errorbar(np.array(period_L)+0.1,c_L,yerr=dc_L,marker='o',zorder=0,label='Love average')
        # plt.errorbar(np.array(period_R_obs)-0.1,c_R_obs,yerr=dc_R_obs,marker='o',zorder=0,label='Rayleigh observed')
        # plt.errorbar(np.array(period_L_obs)+0.2,c_L_obs,yerr=dc_L_obs,marker='o',zorder=0,label='Love observed')
        # plt.legend()
        
        if len(Modes_R)>0 and len(Modes_L)>0:
            plt.legend([ave_R, obs_R, ave_L, obs_L], ['Rayleigh average', 'Rayleigh observed', 'Love average', 'Love observed'])
        elif len(Modes_R)>0 and len(Modes_L)==0:
            plt.legend([ave_R, obs_R], ['Rayleigh average', 'Rayleigh observed'])
        elif len(Modes_R)==0 and len(Modes_L)>0:
            plt.legend([ave_L, obs_L], ['Love average', 'Love observed'])

        
        # plt.xlim([40,360])
        plt.ylim([2.5,9.5])
        plt.xlabel('Period (s)')
        plt.ylabel('Phase Velocity (m/s)')
        plt.title('Compare Dispersion Curves')

        plt.xlim([np.min(period_R+period_L),np.max(period_R+period_L)])

        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Dispersion.png',dpi=200)
        plt.close()


################################# Average PsPp Fit ################################################
if PsPp_Fit:
    if os.path.isfile(directory+'/'+'PsPp_mean.out'):
        file=open(directory+'/'+'PsPp_mean.out','r')
        lines=file.readlines()
        file.close()

        nrays=int(lines[0].split()[0])

        rayp=[]
        d_PsPp=[]
        d_PsPpe=[]
        for line in lines[1:nrays+1]:
            data=line.split()
            rayp.append(float(data[0]))
            d_PsPp.append(float(data[1]))
            d_PsPpe.append(float(data[2]))

        rayp=np.array(rayp)
        d_PsPp=np.array(d_PsPp)
        d_PsPpe=np.array(d_PsPpe)

        mean=np.zeros((nrays,3))
        mean[:,0]=rayp[:]
        mean[:,1]=d_PsPp[:]
        mean[:,2]=d_PsPpe[:]
        sorted_mean = mean[np.argsort(mean[:, 0])]

        # true dispersion curves (data)
        file=open(directory+'/'+'PsPp_obs.out','r')
        lines=file.readlines()
        file.close()

        nrays=int(lines[0].split()[0])

        rayp_obs=[]
        d_obsPsPp=[]
        d_obsPsPpe=[]
        for line in lines[1:nrays+1]:
            data=line.split()
            rayp_obs.append(float(data[0]))
            d_obsPsPp.append(float(data[1]))
            d_obsPsPpe.append(float(data[2]))

        rayp_obs=np.array(rayp_obs)
        d_obsPsPp=np.array(d_obsPsPp)
        d_obsPsPpe=np.array(d_obsPsPpe)

        obs=np.zeros((nrays,3))
        obs[:,0]=rayp_obs[:]
        obs[:,1]=d_obsPsPp[:]
        obs[:,2]=d_obsPsPpe[:]
        sorted_obs = obs[np.argsort(obs[:, 0])]

        plt.figure('PsPp_Fit')

        ave_PsPp=plt.errorbar(sorted_mean[:,0],sorted_mean[:,1],yerr=sorted_mean[:,2],marker='o',markersize=3,zorder=0,label='PsPp average',mfc='blue',mec='blue', c='blue')
        obs_PsPp=plt.errorbar(sorted_obs[:,0] ,sorted_obs[:,1], yerr=sorted_obs[:,2], marker='o',markersize=3,zorder=0,label='PsPp observed',mfc='limegreen',mec='limegreen', c='limegreen')
        
        plt.legend([ave_PsPp, obs_PsPp], ['PsPp average', 'PsPp observed'])

        plt.xlabel('Ray Parameter (s/deg)')
        plt.ylabel('PsPp TT (s)')
        plt.title('Compare PsPp Fit')

        plt.xlim([np.min(rayp_obs),np.max(rayp_obs)])
        
        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'PsPp_Fit.png',dpi=200)
        plt.close()


if Corr_Hist:
    if os.path.isfile(directory+'/'+'Posterior_corr.json'):
        print('Posterior Correlation file present')
        
        with open(directory+'/'+'Posterior_corr.json') as data_file:
            corrs = json.load(data_file)

        nsample  = int(corrs['params_inversion']['nsample'])
        thinning = int(corrs['params_inversion']['thinning'])
        cores    = int(corrs['params_inversion']['cores'])
        # print(nsample, thinning, cores)
        # print(corrs['stack']['200']['50'].keys())

        if not str(av_int) in corrs['nostack']:
            print('Bad av_int supplied...')
            print('Options [2] :     Averaging interval, 50,100,200')
            print('Options [2] :    ',corrs['nostack'].keys())
            exit('exiting....')
        
        if not p1_u_dep in corrs['nostack'][str(av_int)]['av_u_deps'] or not p2_u_dep in corrs['nostack'][str(av_int)]['av_u_deps']:
            print('Bad upper layer depths supplied...')
            print('Options [5/7] :     Parameter 1/2 av_u_dep: 100, 200 ....')
            print('Options [2] :    ',corrs['nostack'][str(av_int)]['av_u_deps'])
            exit('exiting....')

        print('Plotting Corr hist for av_int: '+str(av_int)+', dh: '+str(dh)+', p1: '+str(p1)+' '+str(p1_u_dep)+'km, p2: '+str(p2)+' '+str(p2_u_dep)+'km')

        # Script.py <dir> <av_int> <dis> <p1> <av_u_dep1> <p2> <av_u_dep2>

        p1_u_dep_ind=np.where(np.array(corrs['nostack'][str(av_int)]['av_u_deps'])==float(p1_u_dep))[0][0]
        p2_u_dep_ind=np.where(np.array(corrs['nostack'][str(av_int)]['av_u_deps'])==float(p2_u_dep))[0][0]

        if p1 == 'vph':
            # min_1=np.round_(float(corrs['params_inversion']['vpref_min']),2)
            # max_1=np.round_(float(corrs['params_inversion']['vpref_max']),2)
            min_1=np.round_(float(corrs['nostack'][str(av_int)]['av_vph_min'][p1_u_dep_ind]),2)
            max_1=np.round_(float(corrs['nostack'][str(av_int)]['av_vph_max'][p1_u_dep_ind]),2)
            p1_l='$V_{PH}$'

        if p1 == 'vsv':
            # min_1=np.round_(float(corrs['params_inversion']['vsref_min']),2)
            # max_1=np.round_(float(corrs['params_inversion']['vsref_max']),2)
            min_1=np.round_(float(corrs['nostack'][str(av_int)]['av_vsv_min'][p1_u_dep_ind]),2)
            max_1=np.round_(float(corrs['nostack'][str(av_int)]['av_vsv_max'][p1_u_dep_ind]),2)

            p1_l='$V_{SV}$'

        if p2 == 'vsv':
            # min_2=np.round_(float(corrs['params_inversion']['vsref_min']),2)
            # max_2=np.round_(float(corrs['params_inversion']['vsref_max']),2)
            min_2=np.round_(float(corrs['nostack'][str(av_int)]['av_vsv_min'][p2_u_dep_ind]),2)
            max_2=np.round_(float(corrs['nostack'][str(av_int)]['av_vsv_max'][p2_u_dep_ind]),2)
            p2_l='$V_{SV}$'

        if p2 == 'xi':
            min_2=np.round_(float(corrs['params_inversion']['xi_min']),2)
            max_2=np.round_(float(corrs['params_inversion']['xi_max']),2)
            p2_l='$Xi$'

        p1_l_dep = p1_u_dep+av_int
        p2_l_dep = p2_u_dep+av_int

        corr_name='h_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        p1_p2 = np.array(corrs['stack'][str(av_int)][str(dh)][str(corr_name)])

        # Find max value.
        max_loc_p1_p2=unravel_index(p1_p2.argmax(), p1_p2.shape)
        # To remove the dominance of the isotropic layers.
        print(max_loc_p1_p2)
        # if p2 == 'xi':
        #     p1_p2[:,max_loc_p1_p2[1]]=p1_p2[:,max_loc_p1_p2[1]-1]

        vals = p1_p2 * (nsample/thinning) * cores

        # print(min_1, max_1, min_2, max_2, np.shape(p1_p2))

        # Define the x and y data 
        y = np.sum(p1_p2,axis=1)
        x = np.sum(p1_p2,axis=0)

        y_array = np.linspace(min_1+(((max_1-min_1)/dh)/2), max_1-(((max_1-min_1)/dh)/2), dh) # Center of the bins
        x_array = np.linspace(min_2+(((max_2-min_2)/dh)/2), max_2-(((max_2-min_2)/dh)/2), dh)

        # print(x_array[dh//2-1])
        # print(x_array)
        # print(y_array)

        # Generate the actual data to make stats easier
        vals1=[]
        vals2=[]

        for i in range(0,dh):
            for j in range(0,dh):
                val=int(vals[i,j])
                for k in range(0,val):
                    vals2.append(x_array[j]) # vsv or xi
                    vals1.append(y_array[i]) # vph or vsv

        # Calculate stats
        r,p=stats.pearsonr(vals1, vals2)

        y_mean=np.round_(np.mean(vals1),4)
        y_median=np.round_(np.median(vals1),4)
        y_975=np.round_(np.quantile(vals1, 0.975),4)
        y_825=np.round_(np.quantile(vals1, 0.825),4)
        y_175=np.round_(np.quantile(vals1, 0.175),4)
        y_025=np.round_(np.quantile(vals1, 0.025),4)
        y_std=np.round_(np.std(vals1),4)
        x_mean=np.round_(np.mean(vals2),4)
        x_median=np.round_(np.median(vals2),4)
        x_975=np.round_(np.quantile(vals2, 0.975),4)
        x_825=np.round_(np.quantile(vals2, 0.825),4)
        x_175=np.round_(np.quantile(vals2, 0.175),4)
        x_025=np.round_(np.quantile(vals2, 0.025),4)
        x_std=np.round_(np.std(vals2),4)

        print('x_mean: '+str(x_mean)+', x_std: '+str(x_std)+', x_median: '+str(x_median)+', 2.5%: '+str(x_025)+', 97.5%: '+str(x_975))
        print('y_mean: '+str(y_mean)+', y_std: '+str(y_std)+', y_median: '+str(y_median)+', 2.5%: '+str(y_025)+', 97.5%: '+str(y_975))
        print('Pearson r: '+str(r))


        # Normalize all max values to 1. Above all values sum to 1.0.
        p1_p2=p1_p2/np.max(p1_p2[:,:])
        x=x/np.max(x[:])
        y=y/np.max(y[:])

        # Set up your x and y labels
        ylabel = str(p1_l)+' (m/s) '+str(int(p1_u_dep))+' - '+str(int(p1_l_dep))+' km'
        if p2 == 'xi':
            xlabel = str(p2_l)+' '+str(int(p2_u_dep))+' - '+str(int(p2_l_dep))+' km' 
        else:
            xlabel = str(p2_l)+' (m/s) '+str(int(p2_u_dep))+' - '+str(int(p2_l_dep))+' km' 
        # Define the locations for the axes
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = left_h = left+width+0.05
        
        # Set up the geometry of the three plots
        rect_corrxy = [left, bottom, width, height] # dimensions of temp plot
        rect_histx = [left, bottom_h, width, 0.25] # dimensions of x-histogram
        rect_histy = [left_h, bottom, 0.25, height] # dimensions of y-histogram
        
        # Set up the size of the figure
        fig = plt.figure('Corr_hist', figsize=(10,9.5))
        # plt.title('Correlation')

        # Make the three plots
        axCorr = fig.add_axes(rect_corrxy) # temperature plot
        axHistx = fig.add_axes(rect_histx) # x histogram
        axHisty = fig.add_axes(rect_histy) # y histogram
        
        # Remove the inner axes numbers of the histograms
        nullfmt = NullFormatter()
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        # Set up default x and y limits
        xlims = [0, dh]
        ylims = [0, dh]
            
        # Find the min/max of the data
        xmin = min(xlims)
        xmax = max(xlims)
        ymin = min(ylims)
        ymax = max(ylims)
        # print(xmin,xmax,ymin,ymax)


        # option 1
        # axCorr.imshow(p1_p2, extent=[ymin,ymax,xmin,xmax], aspect='auto',
        #     interpolation='nearest', origin='lower')
        # Option 2
        # axCorr.contourf(p1_p2, extent=[ymin-1,ymax+1,xmin-1,xmax+1], origin='lower')
        # Option 3
        # x_vals=[]; y_vals=[]; z_vals=[]
        # for i in range(0,dh):
        #     for j in range(0,dh):
        #         y_vals.append(y_array[i])
        #         x_vals.append(x_array[j])
        #         z_vals.append(int(vals[i,j]))
        # hist, bins_y, bins_x=np.histogram2d(y_vals, x_vals, weights=z_vals, bins=(ymax,xmax))
        # try:
        #     probs=[0, 0.8, 0.90, 0.95, 0.98, 0.99, 0.995,0.998]
        #     q = stats.mstats.mquantiles(hist,prob=probs)
        #     m=axCorr.contourf(hist, extent=[ymin-1,ymax+1,xmin-1,xmax+1], origin='lower', levels=q, cmap='viridis')
        # except:
        #     try:
        #         probs=[0, 0.85, 0.90, 0.95, 0.98, 0.99, 0.995,0.998]
        #         q = stats.mstats.mquantiles(hist,prob=probs)
        #         m=axCorr.contourf(hist, extent=[ymin-1,ymax+1,xmin-1,xmax+1], origin='lower', levels=q, cmap='viridis')
        #     except:
        #         try:
        #             probs=[0, 0.90, 0.95, 0.98, 0.99, 0.995,0.998]
        #             q = stats.mstats.mquantiles(hist,prob=probs)
        #             m=axCorr.contourf(hist, extent=[ymin-1,ymax+1,xmin-1,xmax+1], origin='lower', levels=q, cmap='viridis')
        #         except:
        #             probs=[0, 0.95, 0.98, 0.99, 0.995,0.998]
        #             q = stats.mstats.mquantiles(hist,prob=probs)
        #             m=axCorr.contourf(hist, extent=[ymin-1,ymax+1,xmin-1,xmax+1], origin='lower', levels=q, cmap='viridis')

        # axins1 = inset_axes(axCorr,
        #             width="60%",  # width = 40% of parent_bbox width
        #             height="5%",  # height : 5%
        #             loc='lower right')
        # cbar=plt.colorbar(m, cax=axins1,orientation="horizontal", ticks=q)
        # cbar.ax.set_xticklabels(probs, fontsize=12, rotation=90, color='w')
        # cbar.ax.set_title(label='Probability', color='w')

        # axins1.xaxis.set_ticks_position("top")
        # axins1.xaxis.set_label_position("top")


        # Option 4

        xmin=np.min(vals2)
        xmax=np.max(vals2)
        ymax=np.max(vals1)
        ymin=np.min(vals1)
        print('Here.... #1')
        values = np.vstack([vals2,vals1])
        kernel = stats.gaussian_kde(values)
        print('Calculated KDE  #2')
        interval=50
        X, Y = np.mgrid[min_2:max_2:50j, min_1:max_1:50j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        print('Calculated meshgrid.... #3')
        Z = np.reshape((1/np.size(positions))*kernel.evaluate(positions), X.shape)
        max_Z=unravel_index(Z.argmax(), Z.shape)
        print('Re-shaped Z-array.... #4')
        m=axCorr.contourf(X,Y, Z, extent=[xmin, xmax, ymin, ymax], cmap='viridis',levels=np.linspace(0,np.max(Z),25))
        print('Contoured Z-array.... #5')
        print(X[max_Z],Z[max_Z])
        axCorr.plot(X[max_Z],Y[max_Z],c='r',marker='x',markersize=5,linewidth=1, alpha=1) # centered using +0.5

        axCorr.set_xlim([min_2, max_2])
        axCorr.set_ylim([min_1, max_1])

        xlabels=[ str(i) for i in np.linspace(min_2, max_2,3)]
        ylabels=[ str(i) for i in np.linspace(min_1, max_1,3)]


        axCorr.xaxis.set_ticks(np.linspace(min_2, max_2,3))
        axCorr.xaxis.set_ticklabels(xlabels, fontsize=12)
        axCorr.yaxis.set_ticks(np.linspace(min_1, max_1,3))
        axCorr.yaxis.set_ticklabels(ylabels, fontsize=12)


        axins1 = inset_axes(axCorr,
                    width="40%",  # width = 40% of parent_bbox width
                    height="5%",  # height : 5%
                    loc='lower right')
        cbar=plt.colorbar(m, cax=axins1, orientation="horizontal")
        cbar.ax.locator_params(nbins=5)
        # COLORBAR
        fg_color = 'white'
        # set colorbar tick color
        cbar.ax.xaxis.set_tick_params(color=fg_color)
        # set colorbar edgecolor 
        cbar.outline.set_edgecolor(fg_color)
        # set colorbar ticklabels
        cbar.ax.xaxis.set_tick_params(color=fg_color, labelcolor=fg_color, labelsize=14)
        cbar.ax.set_title(label='KDE - Probability', color=fg_color, fontsize=14)

        axins1.xaxis.set_ticks_position("top")
        axins1.xaxis.set_label_position("top")

        axCorr.plot([min_2, max_2], [min_1, max_1], '--g', linewidth=1) # 

        #Plot the axes labels
        axCorr.set_xlabel(xlabel,fontsize=16)
        axCorr.set_ylabel(ylabel,fontsize=16)

        # xlabels=[ str(i) for i in np.linspace(min_2, max_2,5)]
        # ylabels=[ str(i) for i in np.linspace(min_1, max_1,5)]
        # axCorr.xaxis.set_ticks(np.linspace(0, dh, 5))
        # axCorr.xaxis.set_ticklabels(xlabels, fontsize=12)
        # axCorr.yaxis.set_ticks(np.linspace(0, dh, 5))
        # axCorr.yaxis.set_ticklabels(ylabels, fontsize=12)
        # axCorr.set_xlim( 0, dh )
        # axCorr.set_ylim( 0, dh )
        # This works if it is 100 but not if 200, or 50?
        x_pos=np.linspace(0.5,len(x)-0.5, len(x))
        y_pos=np.linspace(0.5,len(y)-0.5, len(y))
        # # #Plot the histograms

        # Interpolate for a smooth curve and continuous plotting using fill-between
        int_arr=np.linspace(0,len(x), 10000)
        x_int=np.interp(int_arr,x_pos,x)
        y_int=np.interp(int_arr,y_pos,y)
        axHistx.plot(int_arr,x_int, 'b-')
        axHisty.plot(y_int,int_arr, 'r-')

        # #Set up the histogram limits
        axHistx.set_xlim( 0, dh )
        axHistx.set_ylim( 0, 1 )

        axHisty.set_ylim( 0, dh )
        axHisty.set_xlim( 0, 1 )

        axHistx.fill_between(int_arr,x_int,0, where=x_int >= 0,facecolor=[0.0, 0., 1.0], rasterized=False, alpha=0.1)
        axHisty.fill_between(y_int,int_arr,0, where=y_int >= 0,facecolor=[1.0, 0., 0.], rasterized=False, alpha=0.1)

        # Calculate the position of the statistics
        x_mean_pos=(x_mean-min_2)/(max_2-min_2)*dh
        x_std_pos=(x_std)/(max_2-min_2)*dh
        x_median_pos=(x_median-min_2)/(max_2-min_2)*dh
        x_025_pos=(x_025-min_2)/(max_2-min_2)*dh
        x_175_pos=(x_175-min_2)/(max_2-min_2)*dh
        x_825_pos=(x_825-min_2)/(max_2-min_2)*dh
        x_975_pos=(x_975-min_2)/(max_2-min_2)*dh


        y_mean_pos=(y_mean-min_1)/(max_1-min_1)*dh
        y_std_pos=(y_std)/(max_1-min_1)*dh
        y_median_pos=(y_median-min_1)/(max_1-min_1)*dh
        y_025_pos=(y_025-min_1)/(max_1-min_1)*dh
        y_175_pos=(y_175-min_1)/(max_1-min_1)*dh
        y_825_pos=(y_825-min_1)/(max_1-min_1)*dh
        y_975_pos=(y_975-min_1)/(max_1-min_1)*dh

        ############### Lets add some labels here!!!!

        axCorr.annotate('(a)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
        axHistx.annotate('(b)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))
        axHisty.annotate('(c)',(0, 1),xytext=(5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='left', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))

        axCorr.annotate('r={:.2f}'.format(r),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top', bbox=dict(facecolor='white',edgecolor='black', pad=2.0))

        axHistx.annotate('-- median={:.2f}'.format(x_median),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
        axHistx.annotate('- mean={:.2f}'.format(x_mean),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
        # axHistx.annotate('std={:.2f}'.format(x_std),(1, 1),xytext=(-5,-35),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

        axHisty.annotate('-- median={:.2f}'.format(y_median),(1, 1),xytext=(-5,-5),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
        axHisty.annotate('- mean={:.2f}'.format(y_mean),(1, 1),xytext=(-5,-20),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')
        # axHisty.annotate('std={:.2f}'.format(y_std),(1, 1),xytext=(-5,-35),xycoords='axes fraction',fontsize=12,textcoords='offset points', color='k', backgroundcolor='none',ha='right', va='top')

        # Plot mean and std on the plots.
        axHistx.plot([x_mean_pos,x_mean_pos],[0, 1],c='b',linewidth=1, alpha=1)
        axHistx.plot([x_median_pos,x_median_pos],[0, 1],c='b', linestyle='--', linewidth=1, alpha=1)

        # For 95% quantile
        a=x_025_pos<=int_arr
        b=x_975_pos>=int_arr
        c=a==b
        axHistx.fill_between(int_arr,x_int,0, where=c,facecolor=[0.0, 0.0, 1.0], rasterized=False, alpha=0.1)
        # For 65% quantile
        a=x_175_pos<=int_arr
        b=x_825_pos>=int_arr
        c=a==b
        axHistx.fill_between(int_arr,x_int,0, where=c,facecolor=[0.0, 0.0, 1.0], rasterized=False, alpha=0.1)

        axHisty.plot([0, 1],[y_mean_pos,y_mean_pos],c='r',linewidth=1, alpha=1)
        axHisty.plot([0, 1],[y_median_pos,y_median_pos],c='r', linestyle='--', linewidth=1, alpha=1)

        # For 95% quantile
        a=y_025_pos<=int_arr
        b=y_975_pos>=int_arr
        c=a==b
        axHisty.fill_betweenx(int_arr,0,y_int, where=c,facecolor=[1.0, 0.0, 0.0], rasterized=False, alpha=0.1)
        # For 65% quantile
        a=y_175_pos<=int_arr
        b=y_825_pos>=int_arr
        c=a==b
        axHisty.fill_betweenx(int_arr,0,y_int, where=c,facecolor=[1.0, 0.0, 0.0], rasterized=False, alpha=0.1)



        #################### Plot Max points and LSQR / inital model points
        # Max position: x, r
        axCorr.plot(max_loc_p1_p2[1]+0.5,max_loc_p1_p2[0]+0.5,c='r',marker='x',markersize=5,linewidth=1, alpha=1) # centered using +0.5
        
        Modes='fund'
        # Modes='over'

        # ########################### Add PREM to plot #########################################
        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/PREM_4_MINEOS/PREM_corr.json'):
        #     print('Adding PREM as Green circle') # o, g
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/PREM_4_MINEOS/PREM_corr.json') as data_file:
        #         PREM_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     PREM_p1_p2 = np.array(PREM_corr['stack'][str(av_int)][str(corr_name)])
        #     print(PREM_p1_p2[0])
        #     axCorr.plot(PREM_p1_p2[0][1],PREM_p1_p2[0][0],c='g',marker='o',markersize=5,linewidth=1, alpha=1)
        # ########################### Add AK135 to plot #########################################
        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135sph_corr.json'):
        #     print('Adding AK135 as Grey triangle') # t, g
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135sph_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='gray',marker='^',markersize=5,linewidth=1, alpha=1)

        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_ak135_corr.json'):
        #     print('Adding AK135 inversion as Grey inverted triangle') # t, g
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_ak135_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='gray',marker='v',markersize=5,linewidth=1, alpha=1)


        # ########################### Add Craton to plot #########################################
        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135_craton_corr.json'):
        #     print('Adding AK135_craton input as blue triangle') # t, b
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135_craton_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='b',marker='^',markersize=5,linewidth=1, alpha=1)

        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_craton_corr.json'):
        #     print('Adding RM_craton inversion as blue inverted triangle') # v, b
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_craton_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='b',marker='v',markersize=5,linewidth=1, alpha=1)

        # ########################### Add Craton to plot #########################################
        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135_craton_meta_corr.json'):
        #     print('Adding AK135_craton_meta input as red triangle') # t, br
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/AK135_craton_meta_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='r',marker='^',markersize=5,linewidth=1, alpha=1)

        # if os.path.isfile('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_craton_meta_corr.json'):
        #     print('Adding RM_craton_meta inversion as red inverted triangle') # v, r
        #     with open('/Users/alistair/Google_Drive/Lyon_Pdoc/Surf_rad_ani/AL_surf_inv_'+str(Modes)+'/RM_craton_meta_corr.json') as data_file:
        #         AK135_corr = json.load(data_file)
        #     corr_name='corr_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)
        #     AK135_p1_p2 = np.array(AK135_corr['stack'][str(av_int)][str(corr_name)])
        #     print(AK135_p1_p2[0])
        #     axCorr.plot(AK135_p1_p2[0][1],AK135_p1_p2[0][0],c='r',marker='v',markersize=5,linewidth=1, alpha=1)




        #Make the tickmarks pretty
        ticklabels = axHistx.get_yticklabels()
        for label in ticklabels:
            label.set_fontsize(12)
            label.set_family('sans-serif')
        
        #Make the tickmarks pretty
        ticklabels = axHisty.get_xticklabels()
        for label in ticklabels:
            label.set_fontsize(12)
            label.set_family('sans-serif')
        
        # #Cool trick that changes the number of tickmarks for the histogram axes
        axHisty.xaxis.set_major_locator(LinearLocator(3))
        axHistx.yaxis.set_major_locator(LinearLocator(3))
        axHisty.yaxis.set_major_locator(LinearLocator(5))
        axHistx.xaxis.set_major_locator(LinearLocator(5))


        axCorr.yaxis.set_major_locator(LinearLocator(5))
        axCorr.xaxis.set_major_locator(LinearLocator(5))

        # plt.show()
        plt.savefig(directory+'/PLOTS/'+str(fname_pre)+'Corr_av_'+str(av_int)+'_dh_'+str(dh)+'_'+str(p1)+'_'+str(p1_u_dep)+'_'+str(p2)+'_'+str(p2_u_dep)+'.png',dpi=200)
        plt.close()


