#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###################################### RUN IMPORTS ############################
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
from matplotlib import rcParams
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
matplotlib.rcParams['font.size'] = 10
import os, sys


print('Required :       Script.py test_directory')
print('Arguments:      ',sys.argv)
print('Options [1] :     RJ_MCMC_Tests/XYZ_test/OUT_TEST')
# python plot_test_results.py ./OUT_TEST2

num_args=len(sys.argv)
if num_args < 2:
    print('Number of arguments (' + str(num_args) +') too low... exit')
    exit('exiting....')
    
directory = str(sys.argv[1])
print('Plotting results for: ' +str(directory))


####################### SET PARAMS ############################################
# input directory, contains all the outputs from the MCMC code
# directory='OUT_TEST2'
if not os.path.exists(directory+'/PLOTS/'):
    os.mkdir(directory+'/PLOTS/')
maxnlay = 80

Layer_Hist      = True
Layers_Aniso_Tr = True
Posterior       = True
Sigmad          = True
Dispersion      = True
Convergence     = True

########## histogram of number of layers
if Layer_Hist:
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
    plt.savefig(directory+'/PLOTS/'+'Layer_hist.pdf')
    plt.close()


########## matrix containing number of layers and number of anisotropic layers
if Layers_Aniso_Tr:
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
    plt.savefig(directory+'/PLOTS/'+'Layers_Aniso_Tr.pdf')
    plt.close()


############ Posterior #######################
if Posterior:
    file=open(directory+'/'+'Posterior.out','r')
    depthplot=60

    lines=file.readlines()

    file.close()

    nlines = len(lines)
    data0=lines[0].split()
    prof=float(data0[0])
    disd=int(data0[1]) # Depth discretisation
    [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[1].split()]
    disv=int(disv) # Velocity/aniso discretisation

    vsvd=np.zeros((disd,disv))
    xid=np.zeros((disd,disv))
    vpd=np.zeros((disd,disv))

    vsvs=np.linspace(vref_min,vref_max,disv+1)
    xis=np.linspace(xi_min,xi_max,disv+1)
    vps=np.linspace(vp_min,vp_max,disv+1)
    depths=np.linspace(0,prof,disd+1)

    i=0
    j=0

    depth=[]

    for line in lines[2:]:
        [vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]

        j+=1
        if j==disv:
            j=0
            i+=1
    xid[:,disv//2-1]=0 # else the isotropic layers dominate the anisotropy density plot

    for i in range(np.shape(vsvd)[0]): # normalisation by depth, not strictly necessary
        s=np.amax(vsvd[i,:])
        vsvd[i,:]=vsvd[i,:]/s

    for i in range(np.shape(vpd)[0]): # normalisation by depth, not strictly necessary
        s=np.amax(vpd[i,:])
        vpd[i,:]=vpd[i,:]/s

    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))


    ax0.invert_yaxis()
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('Vsv (km/s)', fontsize=10)
    ax0.set_title('S-wave velocity')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel(r'xi',fontsize=10)
    ax1.set_xlim([xi_min,xi_max])
    ax1.set_title('Radial Anisotropy')
    ax2.set_xlabel(r'Vp (km/s)',fontsize=10)
    ax2.set_xlim([vp_min,vp_max])
    ax0.set_ylabel('Depth (km)',fontsize=10)
    ax2.set_xlabel(r'Vp, (km/s)',fontsize=10)
    ax2.set_title('P-wave velocity')
    ax1.pcolormesh(xis,depths,xid,cmap='viridis')
    ax0.pcolormesh(vsvs,depths,vsvd,cmap='viridis')

    # true model overlaid on the posterior (only for synthetic tests)
    file=open(directory+'/'+'true_model.out','r')
    lines=file.readlines()
    file.close()

    true_depth=[]
    true_vsv=[]
    true_xi=[]
    true_vpvs=[]
    for line in lines:
        data=line.split()
        true_depth.append(float(data[0]))
        true_vsv.append(float(data[1]))
        try:
            true_xi.append(float(data[2]))
            true_vpvs.append(float(data[3]))
        except:
            pass


    ax0.plot(true_vsv,true_depth,c='c',linewidth=3)
    try:
        ax1.plot(true_xi,true_depth,c='c',linewidth=3)
        ax2.plot(true_vpvs,true_depth,c='c',linewidth=3)
    except:
        pass

    ax2.pcolormesh(vps,depths,vpd,cmap='viridis')
    plt.setp(ax2.get_yticklabels(), visible=False)


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
    file=open(directory+'/'+'Average.out','r')
    lines=file.readlines()
    file.close()

    depths=[]
    average_vs=[]
    average_xi=[]
    average_vpvs=[]
    average_probani=[]
    for line in lines:
        data=line.split()
        depths.append(float(data[0]))
        average_vs.append(float(data[1]))
        average_xi.append(float(data[2]))
        average_vpvs.append(float(data[3]))
        average_probani.append(float(data[4]))

    ax0.plot(average_vs,depths,c='r',linewidth=3)
    ax1.plot(average_xi,depths,c='r',linewidth=3)
    ax2.plot(average_vpvs,depths,c='r',linewidth=3)
    ax4.plot(average_probani,depths,c='k',linewidth=3)
    ax4.set_xlabel('Probability',fontsize=10)
    ax4.set_title('Anisotropy')

    ax4.set_xlim([0,100])

    fig.suptitle('Posterior and Averages')

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Posterior.pdf')
    plt.close()


#################### histogram of rayleigh uncertainty parameter

if Sigmad:
    
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
    plt.savefig(directory+'/PLOTS/'+'Sigmad_R.pdf')
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
    plt.savefig(directory+'/PLOTS/'+'Sigmad_L.pdf')
    plt.close()

#################################### CONVERGENCE #################################

if Convergence:
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
    for line in lines[1:]:
        data=line.split()
        conv_R_one.append(float(data[0]))
        conv_R.append(float(data[1]))
        conv_L_one.append(float(data[2]))
        conv_L.append(float(data[3]))

    plt.figure('convergence_misfit')
    plt.title('Misfit Convergence')

    plt.plot(conv_R_one[burn_in:],label='Rayleigh, one core')
    plt.plot(conv_R[burn_in:],label='Rayleigh, all cores')
    plt.plot(conv_L_one[burn_in:],label='Love, one core')
    plt.plot(conv_L[burn_in:],label='Love, all cores')
    plt.xlim([0,nsample])
    plt.xlabel('Iteration number')
    plt.ylabel('Misfit')
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_misfit.pdf')
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_layers.pdf')
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
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_sigma_R.pdf')
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
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_sigma_L.pdf')
    plt.close()

    ################ acceptance rates for birth of isotropic and anisotropic layers over time, for one core and average over all cores ###############
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_Birth.pdf')
    plt.close()

    ################ acceptance rates for death of isotropic and anisotropic layers over time, for one core and average over all cores ###############
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_Death.pdf')
    plt.close()

    ################ acceptance rates for change in anisotropy over time, for one core and average over all cores ###############
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_xi.pdf')
    plt.close()

    ################ acceptance rates for change in vp/vsv over time, for one core and average over all cores ###############
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_vp.pdf')
    plt.close()

################################# Average dispersion curves ################################################
if Dispersion:
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

    for R_mode in Modes_R:
        ave_R='ave_R_'+str(R_mode)
        obs_R='obs_R_'+str(R_mode)
        # print(ave_R)
        ind=np.where(n_R==R_mode)
        ave_R=plt.errorbar(np.array(period_R)[ind[0]],np.array(c_R)[ind[0]],yerr=np.array(dc_R)[ind[0]],marker='o',zorder=0,label='Rayleigh average',mfc='blue',mec='blue', c='blue')
        obs_R=plt.errorbar(np.array(period_R_obs)[ind[0]]-0.1,np.array(c_R_obs)[ind[0]],yerr=np.array(dc_R_obs)[ind[0]],marker='o',zorder=0,label='Rayleigh observed',mfc='limegreen',mec='limegreen', c='limegreen')

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
    
    plt.legend([ave_R, obs_R, ave_L, obs_L], ['Rayleigh average', 'Rayleigh observed', 'Love average', 'Love observed'])
    plt.xlabel('Period (s)')
    plt.ylabel('Phase Velocity (km/s)')
    plt.title('Compare Dispersion Curves')

    plt.xlim([np.min(period_R+period_L),np.max(period_R+period_L)])

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Dispersion.pdf')
    # plt.close()
