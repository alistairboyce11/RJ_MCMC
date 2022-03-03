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
import matplotlib.patches as patches


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
PsPp_Fit        = True

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
    plt.savefig(directory+'/PLOTS/'+'Layer_hist.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Layers_Aniso_Tr.png',dpi=200)
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

    vsvs=np.linspace(vref_min,vref_max,disv)
    xis=np.linspace(xi_min,xi_max,disv)
    vps=np.linspace(vp_min,vp_max,disv)
    depths=np.linspace(0,prof,disd)

    i=0
    j=0

    # depth=[]

    for line in lines[2:]:
        [vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]

        j+=1
        if j==disv:
            j=0
            i+=1
    # xid[:,disv//2-1]=0 # else the isotropic layers dominate the anisotropy density plot
    xid[:,disv//2-1]=xid[:,disv//2] # else the isotropic layers dominate the anisotropy density plot

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
    ind=np.where(depths<=dep_val)[0]
    splice_vsvd_a=vsvd[ind,:]
    splice_xid_a=xid[ind,:]
    splice_vpd_a=vpd[ind,:]

    above_splice_vsvd=np.mean(splice_vsvd_a, axis=0)
    above_splice_xid=np.mean(splice_xid_a, axis=0)
    above_splice_vpd=np.mean(splice_vpd_a, axis=0)

    ################# Average distribution BELOW depth

    ind2=np.where(depths>=dep_val)[0]
    splice_vsvd_b=vsvd[ind2,:]
    splice_xid_b=xid[ind2,:]
    splice_vpd_b=vpd[ind2,:]

    below_splice_vsvd=np.mean(splice_vsvd_b, axis=0)
    below_splice_xid=np.mean(splice_xid_b, axis=0)
    below_splice_vpd=np.mean(splice_vpd_b, axis=0)

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

    std_vsvs=np.std(dist_vsvs,axis=0)
    mean_vsvs=np.mean(dist_vsvs,axis=0)
    std_vsvs_a=np.std(dist_vsvs_a,axis=0)
    mean_vsvs_a=np.mean(dist_vsvs_a,axis=0)
    std_vsvs_b=np.std(dist_vsvs_b,axis=0)
    mean_vsvs_b=np.mean(dist_vsvs_b,axis=0)

    # print(mean_vsvs_a,std_vsvs_a,mean_vsvs,std_vsvs,mean_vsvs_b,std_vsvs_b)

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

    std_xis=np.std(dist_xis,axis=0)
    mean_xis=np.mean(dist_xis,axis=0)
    std_xis_a=np.std(dist_xis_a,axis=0)
    mean_xis_a=np.mean(dist_xis_a,axis=0)
    std_xis_b=np.std(dist_xis_b,axis=0)
    mean_xis_b=np.mean(dist_xis_b,axis=0)

    # print(mean_xis_a,std_xis_a,mean_xis,std_xis,mean_xis_b,std_xis_b)

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

    std_vps=np.std(dist_vps,axis=0)
    mean_vps=np.mean(dist_vps,axis=0)
    std_vps_b=np.std(dist_vps_b,axis=0)
    mean_vps_b=np.mean(dist_vps_b,axis=0)
    std_vps_a=np.std(dist_vps_a,axis=0)
    mean_vps_a=np.mean(dist_vps_a,axis=0)

    # print(mean_vps_a,std_vps_a,mean_vps,std_vps,mean_vps_b,std_vps_b)


    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, figsize=(9, 9), sharey=True)

    rect = patches.Rectangle((mean_vsvs_a , 0), std_vsvs_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect2 = patches.Rectangle((mean_vsvs_a-std_vsvs_a , 0), std_vsvs_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax1.add_patch(rect)
    ax1.add_patch(rect2)


    rect3 = patches.Rectangle((mean_vsvs , 0), std_vsvs, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect4 = patches.Rectangle((mean_vsvs-std_vsvs , 0), std_vsvs, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax4.add_patch(rect3)
    ax4.add_patch(rect4)

    rect5 = patches.Rectangle((mean_vsvs_b , 0), std_vsvs_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect6 = patches.Rectangle((mean_vsvs_b-std_vsvs_b , 0), std_vsvs_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax7.add_patch(rect5)
    ax7.add_patch(rect6)




    rect = patches.Rectangle((mean_xis_a , 0), std_xis_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect2 = patches.Rectangle((mean_xis_a-std_xis_a , 0), std_xis_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax2.add_patch(rect)
    ax2.add_patch(rect2)


    rect3 = patches.Rectangle((mean_xis , 0), std_xis, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect4 = patches.Rectangle((mean_xis-std_xis , 0), std_xis, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax5.add_patch(rect3)
    ax5.add_patch(rect4)

    rect5 = patches.Rectangle((mean_xis_b , 0), std_xis_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect6 = patches.Rectangle((mean_xis_b-std_xis_b , 0), std_xis_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax8.add_patch(rect5)
    ax8.add_patch(rect6)

    rect = patches.Rectangle((mean_vps_a , 0), std_vps_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect2 = patches.Rectangle((mean_vps_a-std_vps_a , 0), std_vps_a, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax3.add_patch(rect)
    ax3.add_patch(rect2)


    rect3 = patches.Rectangle((mean_vps , 0), std_vps, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect4 = patches.Rectangle((mean_vps-std_vps , 0), std_vps, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax6.add_patch(rect3)
    ax6.add_patch(rect4)

    rect5 = patches.Rectangle((mean_vps_b , 0), std_vps_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    rect6 = patches.Rectangle((mean_vps_b-std_vps_b , 0), std_vps_b, 1, fill=True, fc='b', edgecolor=None, alpha=0.5, visible=True)
    ax9.add_patch(rect5)
    ax9.add_patch(rect6)

    ax1.plot([mean_vsvs_a,mean_vsvs_a],[0,1],c='red',linewidth=1,alpha=1)
    ax2.plot([mean_xis_a,mean_xis_a],[0,1],c='red',linewidth=1,alpha=1)
    ax3.plot([mean_vps_a,mean_vps_a],[0,1],c='red',linewidth=1,alpha=1)

    ax4.plot([mean_vsvs,mean_vsvs],[0,1],c='red',linewidth=1,alpha=1)
    ax5.plot([mean_xis,mean_xis],[0,1],c='red',linewidth=1,alpha=1)
    ax6.plot([mean_vps,mean_vps],[0,1],c='red',linewidth=1,alpha=1)

    ax7.plot([mean_vsvs_b,mean_vsvs_b],[0,1],c='red',linewidth=1,alpha=1)
    ax8.plot([mean_xis_b,mean_xis_b],[0,1],c='red',linewidth=1,alpha=1)
    ax9.plot([mean_vps_b,mean_vps_b],[0,1],c='red',linewidth=1,alpha=1)


    ax1.plot(vsvs,above_splice_vsvd,c='black',linewidth=1,alpha=1)
    ax2.plot(xis,above_splice_xid,c='black',linewidth=1,alpha=1)
    ax3.plot(vps,above_splice_vpd,c='black',linewidth=1,alpha=1)

    ax4.plot(vsvs,splice_vsvd,c='black',linewidth=1,alpha=1)
    ax5.plot(xis,splice_xid,c='black',linewidth=1,alpha=1)
    ax6.plot(vps,splice_vpd,c='black',linewidth=1,alpha=1)

    ax7.plot(vsvs,below_splice_vsvd,c='black',linewidth=1,alpha=1)
    ax8.plot(xis,below_splice_xid,c='black',linewidth=1,alpha=1)
    ax9.plot(vps,below_splice_vpd,c='black',linewidth=1,alpha=1)

    ax1.set_title('S-vel above '+str(dep_val)+'km')
    ax1.set_xlabel('Vsv (km/s)', fontsize=10)
    ax1.set_ylabel('Probability',fontsize=10)
    ax1.set_xlim([vref_min,vref_max])
    ax1.set_ylim([0.0,1.0])

    ax2.set_title('Rad Anis. above '+str(dep_val)+'km')
    ax2.set_xlabel(r'Xi',fontsize=10)
    ax2.set_xlim([xi_min,xi_max])

    ax3.set_title('Vp/Vs above '+str(dep_val)+'km')
    ax3.set_xlim([vp_min,vp_max])
    ax3.set_xlabel(r'VpVs*(1+X)',fontsize=10)

    ax4.set_title('S-vel at '+str(dep_val)+'km')
    ax4.set_xlabel('Vsv (km/s)', fontsize=10)
    ax4.set_ylabel('Probability',fontsize=10)
    ax4.set_xlim([vref_min,vref_max])
    ax4.set_ylim([0.0,1.0])

    ax5.set_title('Rad Anis. at '+str(dep_val)+'km')
    ax5.set_xlabel(r'Xi',fontsize=10)
    ax5.set_xlim([xi_min,xi_max])

    ax6.set_title('Vp/Vs at '+str(dep_val)+'km')
    ax6.set_xlim([vp_min,vp_max])
    ax6.set_xlabel(r'VpVs*(1+X)',fontsize=10)

    ax7.set_title('S-vel below '+str(dep_val)+'km')
    ax7.set_xlabel('Vsv (km/s)', fontsize=10)
    ax7.set_ylabel('Probability',fontsize=10)
    ax7.set_xlim([vref_min,vref_max])
    ax7.set_ylim([0.0,1.0])

    ax8.set_title('Rad Anis. below '+str(dep_val)+'km')
    ax8.set_xlabel(r'Xi',fontsize=10)
    ax8.set_xlim([xi_min,xi_max])

    ax9.set_title('Vp/Vs below '+str(dep_val)+'km')
    ax9.set_xlim([vp_min,vp_max])
    ax9.set_xlabel(r'VpVs*(1+X)',fontsize=10)




    plt.subplots_adjust(wspace=0.35,hspace=0.35)

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Hist_dist_'+str(dep_val)+'.png',dpi=200)
    plt.close()


    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))

    vsvs=np.linspace(vref_min,vref_max,disv+1)
    xis=np.linspace(xi_min,xi_max,disv+1)
    vps=np.linspace(vp_min,vp_max,disv+1)
    depths=np.linspace(0,prof,disd+1)


    ax0.invert_yaxis()
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('Vsv (km/s)', fontsize=10)
    ax0.set_title('S-wave velocity')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel(r'xi',fontsize=10)
    ax1.set_xlim([xi_min,xi_max])
    # ax1.set_xlim([0.5, 1.5])
    ax1.set_title('Radial Anisotropy')
    ax2.set_xlim([vp_min,vp_max])
    # ax2.set_xlim([-0.5,0.5])
    ax0.set_ylabel('Depth (km)',fontsize=10)
    ax2.set_xlabel(r'VpVs*(1+X)',fontsize=10)
    ax2.set_title('Vp/Vs deviation')
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
    for line in lines[1:]:
        data=line.split()
        true_depth.append(float(data[0]))
        true_vsv.append(float(data[1]))
        try:
            true_xi.append(float(data[2]))
            true_vpvs.append(float(data[3]))
            # pass
        except:
            pass

    # True models in cyan.
    true,=ax0.plot(true_vsv,true_depth,c='white',linewidth=1)
    try:
        ax1.plot(true_xi,true_depth,c='white',linewidth=1,alpha=1,marker='o',markersize=2,mfc='k')
        # ax1.plot(true_xi,true_depth,c='magenta',linewidth=2,alpha=1,marker='o',markersize=2,mfc='k')
        
        ax2.plot(true_vpvs,true_depth,c='white',linewidth=1)
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
    
    # Average models in red.
    ave,=ax0.plot(average_vs,depths,c='r',linewidth=1)
    ax1.plot(average_xi,depths,c='r',linewidth=1)
    ax2.plot(average_vpvs,depths,c='r',linewidth=1)
    ax4.plot(average_probani,depths,c='k',linewidth=1)
    ax4.set_xlabel('Probability',fontsize=10)
    ax4.set_title('Anisotropy')

    ax4.set_xlim([0,100])

    plt.legend([true, ave], ['true', 'ave'])

    fig.suptitle('Posterior and Averages')

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Posterior.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Sigmad_R.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Sigmad_L.png',dpi=200)
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
        plt.savefig(directory+'/PLOTS/'+'Sigmad_PsPp.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_misfit.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_layers.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_sigma_R.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_sigma_L.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_Birth.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_Death.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_xi.png',dpi=200)
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
    plt.savefig(directory+'/PLOTS/'+'Convergence_vp.png',dpi=200)
    plt.close()


    
    ################ acceptance rates for change in vsv (upper half) over time, for one core and average over all cores ###############
    file=open(directory+'/'+'Convergence_vs1.out','r')
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

    plt.figure('convergence_vs1')
    plt.title('Vsv Change Rate Convergence, upper half')
    plt.plot(convsv1_one[burn_in:],label='vsv change upper half, one core')
    plt.plot(convsv1[burn_in:],label='vsv change upper half, all cores')
    plt.xlabel('Iteration number')
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_vs1.png',dpi=200)
    plt.close()
    
    ################ acceptance rates for change in vsv (lower half) over time, for one core and average over all cores ###############
    file=open(directory+'/'+'Convergence_vs2.out','r')
    lines=file.readlines()
    file.close()

    burn_in=int(lines[0].split()[0])
    nsample=int(lines[0].split()[1])

    convsv2_one=[]
    convsv2=[]
    for line in lines[1:]:
        data=line.split()
        convsv2_one.append(float(data[0]))
        convsv2.append(float(data[1]))

    plt.figure('convergence_vs2')
    plt.title('Vsv Change Rate Convergence, lower half')
    plt.plot(convsv2_one[burn_in:],label='vsv change lower half, one core')
    plt.plot(convsv2[burn_in:],label='vsv change lower half, all cores')
    plt.xlabel('Iteration number')
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_vs2.png',dpi=200)
    plt.close()
    
    ################ acceptance rates for change in depth (upper half) over time, for one core and average over all cores ###############
    file=open(directory+'/'+'Convergence_dp1.out','r')
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

    plt.figure('convergence_dp1')
    plt.title('Depth Change Rate Convergence, upper half')
    plt.plot(condp1_one[burn_in:],label='depth change upper half, one core')
    plt.plot(condp1[burn_in:],label='depth change upper half, all cores')
    plt.xlabel('Iteration number')
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_dp1.png',dpi=200)
    plt.close()
    
    ################ acceptance rates for change in depth (lower half) over time, for one core and average over all cores ###############
    file=open(directory+'/'+'Convergence_dp2.out','r')
    lines=file.readlines()
    file.close()

    burn_in=int(lines[0].split()[0])
    nsample=int(lines[0].split()[1])

    condp2_one=[]
    condp2=[]
    for line in lines[1:]:
        data=line.split()
        condp2_one.append(float(data[0]))
        condp2.append(float(data[1]))

    plt.figure('convergence_dp2')
    plt.title('Depth Change Rate Convergence, lower half')
    plt.plot(condp2_one[burn_in:],label='depth change lower half, one core')
    plt.plot(condp2[burn_in:],label='depth change lower half, all cores')
    plt.xlabel('Iteration number')
    plt.ylabel('XXX')
    plt.xlim([0,nsample])
    plt.legend()
    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Convergence_dp2.png',dpi=200)
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
        plt.savefig(directory+'/PLOTS/'+'Convergence_sigma_PsPp.png',dpi=200)
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
    plt.ylabel('Phase Velocity (km/s)')
    plt.title('Compare Dispersion Curves')

    plt.xlim([np.min(period_R+period_L),np.max(period_R+period_L)])

    # plt.show()
    plt.savefig(directory+'/PLOTS/'+'Dispersion.png',dpi=200)
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
        plt.savefig(directory+'/PLOTS/'+'PsPp_Fit.png',dpi=200)
        plt.close()
