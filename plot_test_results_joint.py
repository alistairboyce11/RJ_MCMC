#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:06:13 2021

@author: dorian
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import gc
from matplotlib.backends.backend_pdf import PdfPages

directory='OUT_REAL_2_1_LONG'

results=True
preprocessing=True
posterior=True#True
dispersion=False
tradeoff=False#False
alphas=True
many_plots=False

print(directory+'/PLOTS/')
#print(os.path.exists(directory+'/PLOTS/'))
#print(not os.path.exists(directory+'/PLOTS/'))
if not os.path.exists(directory+'/PLOTS/'):
    print(not os.path.exists(directory+'/PLOTS/'))
    os.mkdir(directory+'/PLOTS/')

if not os.path.exists(directory+'/PLOTS/Posterior/'):
    os.mkdir(directory+'/PLOTS/Posterior/')

if not os.path.exists(directory+'/PLOTS/Tradeoff/'):
    os.mkdir(directory+'/PLOTS/Tradeoff/')
    
if not os.path.exists(directory+'/PLOTS/Dispersion/'):
    os.mkdir(directory+'/PLOTS/Dispersion/')

if not os.path.exists(directory+'/PLOTS/Alphas/'):
    os.mkdir(directory+'/PLOTS/Alphas/')

if not os.path.exists(directory+'/PLOTS/Nlayers/'):
    os.mkdir(directory+'/PLOTS/Nlayers/')

if not os.path.exists(directory+'/PLOTS/SigmaRs/'):
    os.mkdir(directory+'/PLOTS/SigmaRs/')

if not os.path.exists(directory+'/PLOTS/SigmaLs/'):
    os.mkdir(directory+'/PLOTS/SigmaLs/')

if not os.path.exists(directory+'/PLOTS/Individuals/'):
    os.mkdir(directory+'/PLOTS/Individuals/')

maxnlay = 40


##################
# Alpha Histogram 
##################

if results and alphas:

    file=open(directory+'/'+'alpha_hist.out','r')
    lines=file.readlines()
    file.close()
    
    numdis=int(lines[0])
    data=lines[1].split()
    num_logalpha=int(data[0])
    logalpha_min=float(data[1])
    logalpha_max=float(data[2])
    
    alpharange=np.linspace(logalpha_min,logalpha_max,num_logalpha)
    
    alphas=np.zeros((numdis,num_logalpha))
    
    i=0
    j=0
    for line in lines[2:]:
        alphas[i,j]=float(line)
        
        i+=1
        if i==numdis:
            i=0
            j+=1
    
    plt.figure('alphas')
    plt.title('alpha histograms')
    for i in range(numdis):
        plt.plot(alpharange,alphas[i,:],label=float(i+1))
    
    file=open(directory+'/'+'alpharef_hist.out','r')
    lines=file.readlines()
    file.close()
    
    data=lines[0].split()
    num_logalpha=int(data[0])
    logalpha_min=float(data[1])
    logalpha_max=float(data[2])
    
    alpharange=np.linspace(logalpha_min,logalpha_max,num_logalpha)
    
    alpharefs=np.zeros((num_logalpha))
    
    i=0
    for line in lines[1:]:
        alpharefs[i]=float(line)
        
        i+=1
    plt.plot(alpharange,alpharefs,label='ref')
    
    #print(np.argmax(alphas,keepdims=True))
    i=0
    while alphas[0,i]>0:
        i+=1
    print(i)
    print(alphas[0,i])
    print(alphas[0,i-1])
    print(np.argmax(alphas[:,i-1]))
        
    
    #plt.legend()
    plt.savefig(directory+'/PLOTS/'+'alpha_hists_final.pdf')
    
    if many_plots:
        for i in range(numdis):
            print(i)
            fig=plt.figure('alpha '+ str(i+1))
            plt.title('alpha histogram '+ str(i+1)+', too small: '+ str(alphas[i,0]/np.sum(alphas[i,:])*100)+' %')
            plt.plot(alpharange[1:],alphas[i,1:],label=float(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Alphas/alpha_hist_'+str(i)+'.pdf')
            
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/alpha_hist.pdf')
            
            plt.close(fig)

##########################################
# Number of Layers histogram
##########################################

if results:
    # Reference
    
    f=open(directory+'/'+'NB_layers.out','r')
    lines=f.readlines()
    f.close()
    
    malay=int(float(lines[0]))
    n=np.zeros((malay))
    for i,line in enumerate(lines[1:]):
        n[i]=float(line)
    
    plt.figure('nlayers_histogram reference')
    plt.title('nlayers histogram: Reference')
    plt.plot(n)
    
    plt.savefig(directory+'/PLOTS/'+'nb_layers_ref.pdf')
    
    # Widened reference
    
    f=open(directory+'/'+'NB_layers_wide.out','r')
    lines=f.readlines()
    f.close()
    
    malay=int(float(lines[0]))
    n=np.zeros((malay))
    for i,line in enumerate(lines[1:]):
        n[i]=float(line)
    
    plt.figure('nlayers_histogram widened reference')
    plt.title('nlayers histogram: widened reference')
    plt.plot(n)
    plt.savefig(directory+'/PLOTS/'+'nb_layers_wide.pdf')
    
    # Cluster members
    f=open(directory+'/'+'NB_layers_alt.out','r')
    lines=f.readlines()
    f.close()
    
    numdis=int(float(lines[0].split()[1]))
    malay=int(float(lines[0].split()[0]))
    n=np.zeros((malay,numdis))
    for i,line in enumerate(lines[1:]):
        for j,d in enumerate(line.split()):
            n[i,j]=float(d)
    
    plt.figure('nlayers_histogram cluster')
    plt.title('nlayers histogram: Cluster members')
    for i in range(numdis):
        
        plt.plot(n[:,i],label=str(i+1))
    plt.legend()
    plt.savefig(directory+'/PLOTS/'+'nb_layers_cluster.pdf')

################################################
# Tradeoffs
################################################

if results and tradeoff:
    # Tradeoff for reference
    
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
    
    plt.figure('TRA')
    plt.title('Tradeoff reference')
    plt.pcolor(l)
    plt.colorbar()
    
    plt.savefig(directory+'/PLOTS/'+'tradeoff_ref.pdf')
    
    # Tradeoff for widened reference
    
    f=open(directory+'/'+'Tradeoff_wide.out','r')
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
    
    plt.figure('TRA_wide')
    plt.title('Tradeoff widened reference')
    plt.pcolor(l)
    plt.colorbar()
    
    plt.savefig(directory+'/PLOTS/'+'tradeoff_wide.pdf')
    
    # Tradeoff for cluster members
    if many_plots:
        
        
    
        f=open(directory+'/'+'Tradeoff_alt.out','r')
        lines=f.readlines()
        f.close()
        
        nlay=int(float(lines[0].split()[0]))
        numdis=int(float(lines[0].split()[1]))
        l=np.zeros((nlay+1,nlay,numdis))
        
        i=0
        j=0
        for line in lines[1:]:
            #l[j,i,:]=float(line.split())#/(i+1)
            l[j,i,:]=np.array([float(n) for n in line.split()])
            #for n in range(numdis):
            #    l[j,i,n]=float(line_tmp[n])
            j+=1
            if j==nlay+1:
                i+=1
                j=0
        #with PdfPages(directory+'/PLOTS/tradeoffs.pdf') as pdf:
        for i in range(numdis):
            fig=plt.figure('TRA_wide '+ str(i+1))
            plt.title('Tradeoff cluster '+ str(i+1))
            plt.pcolor(l[:,:,i])
            plt.colorbar()
            
            #pdf.savefig()
            plt.savefig(directory+'/PLOTS/Tradeoff/'+'tradeoff_'+str(i)+'.pdf')
            
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/tradeoff.pdf')
            
            plt.close(fig)
            print(i)
            #del fig
        #del l
        #del lines
#gc.collect(generation=2)
#################################################
# Posterior
#################################################
if results and posterior:
    # Reference 
    
    file=open(directory+'/'+'Posterior.out','r')
    depthplot=60
    
    lines=file.readlines()
    
    file.close()
    
    nlines = len(lines)
    data0=lines[0].split()
    prof=float(data0[0])
    disd=int(data0[1])
    [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[1].split()]
    disv=int(disv)
    
    vsvd=np.zeros((disd,disv))
    xid=np.zeros((disd,disv))
    vpd=np.zeros((disd,disv))
    
    vsvs=np.linspace(vref_min,vref_max,disv+1)
    xis=np.linspace(xi_min,xi_max,disv+1)
    vps=np.linspace(vp_min,vp_max,disv+1)
    depths=np.linspace(0,prof,disd)
    #print(disd)
    #print(prof)
    
    i=0
    j=0
    
    
    for line in lines[2:]:
        data=line.split()
        vsvd[i,j]=float(data[0])
        xid[i,j]=float(data[1])
        vpd[i,j]=float(data[2])
        
        #[vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]
        
        j+=1
        if j==disv:
            j=0
            i+=1
    xid[:,disv//2-1]=0 
    
    #print(np.sum(vsvd[30*2,:]))
    #print(np.sum(xid[30*2,:]))
    
    
    # for i in range(np.shape(vsvd)[0]):
    #     s=np.amax(vsvd[i,:])
    #     vsvd[i,:]=vsvd[i,:]/s
    
    # for i in range(np.shape(vpd)[0]):
    #     s=np.amax(vpd[i,:])
    #     vpd[i,:]=vpd[i,:]/s
    
    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))
    
    
    ax0.invert_yaxis()
    #f=open(directory+'Model_REF.out')
    #line=f.readlines()[0]
    #vref=float(line.split()[1])
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('Vsv, km/s')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel(r'xi',fontsize=15)
    ax1.set_xlim([xi_min,xi_max])
    ax2.set_xlabel(r'Vp, km/s',fontsize=15)
    ax2.set_xlim([vp_min,vp_max])
    ax0.set_ylabel('Depth, km',fontsize=15)
    ax2.set_xlabel(r'Vp, km/s',fontsize=15)
    ax0.pcolormesh(vsvs,depths,vsvd[:-1,:],cmap='magma_r')
    ax1.pcolormesh(xis,depths,xid[:-1,:],cmap='magma_r')
    ax2.pcolormesh(vps,depths,vpd[:-1,:],cmap='magma_r')
    
    
    # file=open(directory+'/'+'true_model.out','r')
    # lines=file.readlines()
    # file.close()
    
    # true_depth=[]
    # true_vsv=[]
    # true_xi=[]
    # true_vpvs=[]
    # for line in lines:
    #     data=line.split()
    #     true_depth.append(float(data[0]))
    #     true_vsv.append(float(data[1]))
    #     try:
    #         true_xi.append(float(data[2]))
    #         true_vpvs.append(float(data[3]))
    #     except:
    #         pass
    
    
    # ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
    # try:
    #     ax1.plot(true_xi,true_depth,c='c',linewidth=5)
    #     ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
    # except: 
    #     pass
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    
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
    ax5.set_xlabel('change point histogram',fontsize=15)
    
    file=open(directory+'/'+'Average.out','r')
    lines=file.readlines()
    file.close()
    
    average_vs=np.zeros((disd))
    average_xi=np.zeros((disd))
    average_vpvs=np.zeros((disd))
    average_probani=np.zeros((disd))
    for i,line in enumerate(lines):
        data=line.split()
        average_vs[i]=float(data[1])
        average_xi[i]=float(data[2])
        average_vpvs[i]=float(data[3])
        average_probani[i]=float(data[4])
    
    ax5.plot(change_hist,change_depths)
    ax5.set_xlabel('change point histogram',fontsize=15)
    
    ax0.plot(average_vs,depths,c='r',linewidth=1)
    ax1.plot(average_xi,depths,c='r',linewidth=1)
    ax2.plot(average_vpvs,depths,c='r',linewidth=1)
    ax4.plot(average_probani,depths,c='k',linewidth=1)
    ax4.set_xlabel('anisotropy probability',fontsize=15)
    ax4.set_xlim([0,100])
    
    fig.suptitle('posterior and averages - reference')
    
    plt.savefig(directory+'/PLOTS/'+'posterior_ref.pdf')
    
    # widened Reference 
    
    file=open(directory+'/'+'Posterior_wide.out','r')
    depthplot=60
    
    lines=file.readlines()
    
    file.close()
    
    fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                        figsize=(12, 6))
    
    nlines = len(lines)
    data0=lines[0].split()
    prof=float(data0[0])
    disd=int(data0[1])
    [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[1].split()]
    disv=int(disv)
    
    vsvd=np.zeros((disd,disv))
    xid=np.zeros((disd,disv))
    vpd=np.zeros((disd,disv))
    
    vsvs=np.linspace(vref_min,vref_max,disv+1)
    xis=np.linspace(xi_min,xi_max,disv+1)
    vps=np.linspace(vp_min,vp_max,disv+1)
    depths=np.linspace(0,prof,disd)
    #print(disd)
    #print(prof)
    
    i=0
    j=0
    
    for line in lines[2:]:
        [vsvd[i,j],xid[i,j],vpd[i,j]]=[float(i) for i in line.split()]
        
        j+=1
        if j==disv:
            j=0
            i+=1
    xid[:,disv//2-1]=0 
    
    #print(np.sum(vsvd[30*2,:]))
    #print(np.sum(xid[30*2,:]))
    
    
    # for i in range(np.shape(vsvd)[0]):
    #     s=np.amax(vsvd[i,:])
    #     vsvd[i,:]=vsvd[i,:]/s
    
    # for i in range(np.shape(vpd)[0]):
    #     s=np.amax(vpd[i,:])
    #     vpd[i,:]=vpd[i,:]/s
    
    # file=open(directory+'/'+'true_model.out','r')
    # lines=file.readlines()
    # file.close()
    
    # true_depth=[]
    # true_vsv=[]
    # true_xi=[]
    # true_vpvs=[]
    # for line in lines:
    #     data=line.split()
    #     true_depth.append(float(data[0]))
    #     true_vsv.append(float(data[1]))
    #     try:
    #         true_xi.append(float(data[2]))
    #         true_vpvs.append(float(data[3]))
    #     except:
    #         pass
    
    
    # ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
    # try:
    #     ax1.plot(true_xi,true_depth,c='c',linewidth=5)
    #     ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
    # except: 
    #     pass 
    
    
    ax0.invert_yaxis()
    #f=open(directory+'Model_REF.out')
    #line=f.readlines()[0]
    #vref=float(line.split()[1])
    ax0.set_xlim([vref_min,vref_max])
    ax0.set_xlabel('Vsv, km/s')
    ax1.set_ylim([prof,0.])
    ax1.set_xlabel(r'xi',fontsize=15)
    ax1.set_xlim([xi_min,xi_max])
    ax2.set_xlabel(r'Vp, km/s',fontsize=15)
    ax2.set_xlim([vp_min,vp_max])
    ax0.set_ylabel('Depth, km',fontsize=15)
    ax2.set_xlabel(r'Vp, km/s',fontsize=15)
    ax0.pcolormesh(vsvs,depths,vsvd[:-1,:],cmap='magma_r')
    ax1.pcolormesh(xis,depths,xid[:-1,:],cmap='magma_r')
    ax2.pcolormesh(vps,depths,vpd[:-1,:],cmap='magma_r')
    
    
    # file=open(directory+'/'+'true_model.out','r')
    # lines=file.readlines()
    # file.close()
    
    # true_depth=[]
    # true_vsv=[]
    # true_xi=[]
    # true_vpvs=[]
    # for line in lines:
    #     data=line.split()
    #     true_depth.append(float(data[0]))
    #     true_vsv.append(float(data[1]))
    #     try:
    #         true_xi.append(float(data[2]))
    #         true_vpvs.append(float(data[3]))
    #     except:
    #         pass
    
    
    # ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
    # try:
    #     ax1.plot(true_xi,true_depth,c='c',linewidth=5)
    #     ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
    # except: 
    #     pass
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    
    file=open(directory+'/'+'Change_points_wide.out','r')
    lines=file.readlines()
    file.close()
    
    change_depths=[]
    change_hist=[]
    for line in lines:
        data=line.split()
        change_depths.append(float(data[0]))
        change_hist.append(float(data[1]))
    
    ax5.plot(change_hist,change_depths)
    ax5.set_xlabel('change point histogram',fontsize=15)
    
    file=open(directory+'/'+'Average_wide.out','r')
    lines=file.readlines()
    file.close()
    
    average_vs=np.zeros((disd))
    average_xi=np.zeros((disd))
    average_vpvs=np.zeros((disd))
    average_probani=np.zeros((disd))
    for i,line in enumerate(lines):
        data=line.split()
        average_vs[i]=float(data[1])
        average_xi[i]=float(data[2])
        average_vpvs[i]=float(data[3])
        average_probani[i]=float(data[4])
    
    ax5.plot(change_hist,change_depths)
    ax5.set_xlabel('change point histogram',fontsize=15)
    
    ax0.plot(average_vs,depths,c='r',linewidth=1)
    ax1.plot(average_xi,depths,c='r',linewidth=1)
    ax2.plot(average_vpvs,depths,c='r',linewidth=1)
    ax4.plot(average_probani,depths,c='k',linewidth=1)
    ax4.set_xlabel('anisotropy probability',fontsize=15)
    ax4.set_xlim([0,100])
    
    fig.suptitle('posterior and averages - widened reference')
    
    plt.savefig(directory+'/PLOTS/'+'posterior_wide.pdf')
    
    
    # cluster members 
    
    if many_plots:
        
        file=open(directory+'/'+'Posterior_alt.out','r')
        
        lines=file.readlines()
        
        file.close()
        
        nlines = len(lines)
        numdis=int(float(lines[0]))
        data0=lines[1].split()
        prof=float(data0[0])
        disd=int(float(data0[1]))
        [vref_min,vref_max,disv,width,xi_min,xi_max,vp_min,vp_max]=[float(i) for i in lines[2].split()]
        disv=int(disv)
        
        vsvd=np.zeros((disd,disv,numdis))
        xid=np.zeros((disd,disv,numdis))
        vpd=np.zeros((disd,disv,numdis))
        
        vsvs=np.linspace(vref_min,vref_max,disv+1)
        xis=np.linspace(xi_min,xi_max,disv+1)
        vps=np.linspace(vp_min,vp_max,disv+1)
        depths=np.linspace(0,prof,disd)
        #print(disd)
        #print(prof)
        
        i=0
        j=0
        k=0
        
        
        for line in lines[3:]:
            data=line.split()
            vsvd[i,j,k]=float(data[0])
            xid[i,j,k]=float(data[1])
            vpd[i,j,k]=float(data[2])
            
            k+=1
            if k==numdis:
                j+=1
                k=0
            if j==disv:
                j=0
                i+=1
                
        xid[:,disv//2-1,:]=0 
        
        #print(np.sum(vsvd[30*2,:]))
        #print(np.sum(xid[30*2,:]))
        
        # for k in range(numdis):
        #     for i in range(np.shape(vsvd)[0]):
        #         s=np.amax(vsvd[i,:,k])
        #         vsvd[i,:,k]=vsvd[i,:,k]/s
            
        #     for i in range(np.shape(vpd)[0]):
        #         s=np.amax(vpd[i,:,k])
        #         vpd[i,:,k]=vpd[i,:,k]/s
        
        
        
        
        file=open(directory+'/'+'Change_points_alt.out','r')
        lines=file.readlines()
        file.close()
        
        numdis=int(float(lines[0]))
        change_depths=np.zeros((disd))
        change_hist=np.zeros((disd,numdis))
        
        for i,line in enumerate(lines[1:]):
            data=line.split()
            change_depths[i]=float(data[0])
            for j in range(numdis):
                
                change_hist[i,j]=float(data[j+1])
        
        file=open(directory+'/'+'Average_alt.out','r')
        lines=file.readlines()
        file.close()
        
        numdis=int(float(lines[0]))
        
        average_vs=np.zeros((disd,numdis))
        average_xi=np.zeros((disd,numdis))
        average_vpvs=np.zeros((disd,numdis))
        average_probani=np.zeros((disd,numdis))
        
        i=0
        j=0
        k=0
        for line in lines[3:]:
            average_vs[i,k]=float(line.split()[1])
            average_xi[i,k]=float(line.split()[2])
            average_vpvs[i,k]=float(line.split()[3])
            average_probani[i,k]=float(line.split()[4])
            
            k+=1
            if k==numdis:
                i+=1
                k=0
        
        for k in range(numdis):
            
            if os.path.isfile(directory+'/'+'true_model_'+"{:02d}".format(k+1)+'.out'):
                file=open(directory+'/'+'true_model_'+"{:02d}".format(k+1)+'.out','r')
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
            #with PdfPages(directory+'/PLOTS/posteriors.pdf') as pdf:
            fig, (ax0, ax1, ax2, ax4, ax5) = plt.subplots(nrows=1, ncols=5, sharey=True,
                                                figsize=(12, 6),num='Posterior '+str(k))
            
            
            ax0.invert_yaxis()
            #f=open(directory+'Model_REF.out')
            #line=f.readlines()[0]
            #vref=float(line.split()[1])
            ax0.set_xlim([vref_min,vref_max])
            ax0.set_xlabel('Vsv, km/s')
            ax1.set_ylim([prof,0.])
            ax1.set_xlabel(r'xi',fontsize=15)
            ax1.set_xlim([xi_min,xi_max])
            ax2.set_xlabel(r'Vp, km/s',fontsize=15)
            ax2.set_xlim([vp_min,vp_max])
            ax0.set_ylabel('Depth, km',fontsize=15)
            ax2.set_xlabel(r'Vp, km/s',fontsize=15)
            ax0.pcolormesh(vsvs,depths,vsvd[:-1,:,k],cmap='magma_r')
            ax1.pcolormesh(xis,depths,xid[:-1,:,k],cmap='magma_r')
            ax2.pcolormesh(vps,depths,vpd[:-1,:,k],cmap='magma_r')
            
            
            
            
            if os.path.isfile(directory+'/'+'true_model_'+"{:02d}".format(k+1)+'.out'):
                ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
                try:
                    ax1.plot(true_xi,true_depth,c='c',linewidth=5)
                    ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
                except: 
                    pass
        
            ax5.plot(change_hist[:,k],change_depths)
            ax5.set_xlabel('change point histogram',fontsize=15)
            
            ax0.plot(average_vs[:,k],depths,c='r',linewidth=1)
            ax1.plot(average_xi[:,k],depths,c='r',linewidth=1)
            ax2.plot(average_vpvs[:,k],depths,c='r',linewidth=1)
            ax4.plot(average_probani[:,k],depths,c='k',linewidth=1)
            ax4.set_xlabel('anisotropy probability',fontsize=15)
            ax4.set_xlim([0,100])
            
            plt.setp(ax2.get_yticklabels(), visible=False)
            
            fig.suptitle('posterior and averages '+str(k+1))
            
            plt.savefig(directory+'/PLOTS/Posterior/'+'posterior_'+str(k)+'.pdf')#
            
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/posterior.pdf')
            
            plt.close(fig)
            print(k)
        
        #del vsvd
        #del xid
        #del vpd
        #del lines

#################################################
# Sigma R
#################################################

if results:
    # reference
    
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
    
    plt.figure('sigmad_R ref')
    plt.title('sigmad_R reference')
    plt.plot(d,sigmad_R)
    
    plt.savefig(directory+'/PLOTS/'+'sigma_R_ref.pdf')
    
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
    
    plt.figure('sigmad_L ref')
    plt.title('sigmad_L reference')
    plt.plot(d,sigmad_L)
    plt.savefig(directory+'/PLOTS/'+'sigma_L_ref.pdf')
    
    # widened reference
    
    file=open(directory+'/'+'Sigmad_R_wide.out','r')
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
    plt.title('sigmad_R reference widened')
    plt.plot(d,sigmad_R)
    
    plt.savefig(directory+'/PLOTS/'+'sigma_R_wide.pdf')
    
    file=open(directory+'/'+'Sigmad_L_wide.out','r')
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
    plt.title('sigmad_L reference widened')
    plt.plot(d,sigmad_L)
    
    plt.savefig(directory+'/PLOTS/'+'sigma_L_wide.pdf')
    
    # cluster members
    
    file=open(directory+'/'+'Sigmad_R_alt.out','r')
    lines=file.readlines()
    file.close()
    
    numdis=int(float(lines[0]))
    data_init=lines[1].split()
    ad_r_min=float(data_init[0])
    ad_r_max=float(data_init[1])
    disa=int(float(data_init[2]))
    
    dR=np.zeros((disa))
    sigmad_R=np.zeros((disa,numdis))
    for i,line in enumerate(lines[2:]):
        data=line.split()
        dR[i]=float(data[0])
        for j in range(numdis):
            sigmad_R[i,j]=float(data[j+1])
    
    
    file=open(directory+'/'+'Sigmad_L_alt.out','r')
    lines=file.readlines()
    file.close()
    
    numdis=int(float(lines[0]))
    data_init=lines[1].split()
    ad_l_min=float(data_init[0])
    ad_l_max=float(data_init[1])
    disa=int(float(data_init[2]))
    
    dL=np.zeros((disa))
    sigmad_L=np.zeros((disa,numdis))
    for i,line in enumerate(lines[2:]):
        data=line.split()
        dL[i]=float(data[0])
        for j in range(numdis):
            sigmad_L[i,j]=float(data[j+1])
    
    plt.figure('sigmad_R cluster')
    
    for i in range(numdis):
        
        plt.title('sigmad_R cluster')
        plt.plot(dR,sigmad_R[:,i],label=str(i+1))
        
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'sigma_R_cluster.pdf')
    
    plt.figure('sigmad_L cluster')
    
    for i in range(numdis):
        
        plt.title('sigmad_L cluster')
        plt.plot(dL,sigmad_L[:,i],label=str(i+1))
    
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'sigma_L_cluster.pdf')
    
    # Tradeoff for cluster members
    if many_plots:
        
        #with PdfPages(directory+'/PLOTS/tradeoffs.pdf') as pdf:
        for i in range(numdis):
            fig=plt.figure('sigmad_R '+ str(i+1))
            plt.title('sigmad_R cluster '+ str(i+1))
            plt.plot(dR,sigmad_R[:,i],label=str(i+1))
            
            #pdf.savefig()
            plt.savefig(directory+'/PLOTS/SigmaRs/'+'sigmaR_'+str(i)+'.pdf')
            
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/sigmaR.pdf')
            
            plt.close(fig)
            
            fig=plt.figure('sigmad_L '+ str(i+1))
            plt.title('sigmad_L cluster '+ str(i+1))
            plt.plot(dL,sigmad_L[:,i],label=str(i+1))
            
            #pdf.savefig()
            plt.savefig(directory+'/PLOTS/SigmaLs/'+'sigmaL_'+str(i)+'.pdf')
            
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/sigmaL.pdf')
            
            plt.close(fig)
            #del fig
        #del l
        #del lines

#########################################################
# Convergences
########################################################

if results:
    file=open(directory+'/'+'Convergence_misfit.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
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
    plt.title('misfit convergence')
    #plt.plot(conv_R_one[burn_in:],label='Rayleigh, one core')
    plt.plot(conv_R[burn_in:],label='Rayleigh, all cores')
    #plt.plot(conv_L_one[burn_in:],label='Love, one core')
    plt.plot(conv_L[burn_in:],label='Love, all cores')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_misfit.pdf')
    
    
    file=open(directory+'/'+'Convergence_nb_layers.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    conv_n_one=[]
    conv_n=[]
    for line in lines[1:]:
        data=line.split()
        conv_n_one.append(float(data[0]))
        conv_n.append(float(data[1]))
    
    plt.figure('convergence_nlayers')
    plt.title('number of layers convergence')
    #plt.plot(conv_n_one[burn_in:],label='nblayers, one core')
    plt.plot(conv_n[burn_in:],label='nblayers, all cores')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_np_layers.pdf')
    
    file=open(directory+'/'+'Convergence_sigma_R.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    conv_sigmaR_one=[]
    conv_sigmaR=[]
    for line in lines[1:]:
        data=line.split()
        conv_sigmaR_one.append(float(data[0]))
        conv_sigmaR.append(float(data[1]))
    
    plt.figure('convergence_sigmaR')
    plt.title('sigmaR convergence')
    #plt.plot(conv_sigmaR_one[burn_in:],label='sigmaR, one core')
    plt.plot(conv_sigmaR[burn_in:],label='sigmaR, all cores')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_sigma_R.pdf')
    
    file=open(directory+'/'+'Convergence_sigma_L.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    conv_sigmaL_one=[]
    conv_sigmaL=[]
    for line in lines[1:]:
        data=line.split()
        conv_sigmaL_one.append(float(data[0]))
        conv_sigmaL.append(float(data[1]))
    
    plt.figure('convergence_sigmaL')
    plt.title('sigmaL convergence')
    #plt.plot(conv_sigmaL_one[burn_in:],label='sigmaL, one core')
    plt.plot(conv_sigmaL[burn_in:],label='sigmaL, all cores')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_sigma_L.pdf')
    
    file=open(directory+'/'+'Convergence_Birth.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
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
    plt.title('birth rate convergence')
    #plt.plot(convB_one[burn_in:],label='birth, one core')
    plt.plot(convB[burn_in:],label='birth, all cores')
    #plt.plot(convBa_one[burn_in:],label='birth anisotropic, one core')
    plt.plot(convBa[burn_in:],label='birth anisotropic, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_birth.pdf')
    
    file=open(directory+'/'+'Convergence_Death.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
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
    plt.title('death rate convergence')
    #plt.plot(convD_one[burn_in:],label='death, one core')
    plt.plot(convD[burn_in:],label='death, all cores')
    #plt.plot(convDa_one[burn_in:],label='death anisotropic, one core')
    plt.plot(convDa[burn_in:],label='death anisotropic, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_death.pdf')
    
    file=open(directory+'/'+'Convergence_xi.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    convP_one=[]
    convP=[]
    for line in lines[1:]:
        data=line.split()
        convP_one.append(float(data[0]))
        convP.append(float(data[1]))
    
    plt.figure('convergence_xi')
    plt.title('anisotropy change rate convergence')
    #plt.plot(convP_one[burn_in:],label='xi change, one core')
    plt.plot(convP[burn_in:],label='xi change, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_xi.pdf')
    
    file=open(directory+'/'+'Convergence_vp.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    convP_one=[]
    convP=[]
    for line in lines[1:]:
        data=line.split()
        convP_one.append(float(data[0]))
        convP.append(float(data[1]))
    
    plt.figure('convergence_vp')
    plt.title('vp change rate convergence')
    #plt.plot(convP_one[burn_in:],label='vp change, one core')
    plt.plot(convP[burn_in:],label='vp change, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vp.pdf')
    
    file=open(directory+'/'+'Convergence_vs1.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    convs1_one=[]
    convs1=[]
    for line in lines[1:]:
        data=line.split()
        convs1_one.append(float(data[0]))
        convs1.append(float(data[1]))
    
    plt.figure('convergence_vs1')
    plt.title('vs1 change rate convergence')
    #plt.plot(convs1_one[burn_in:],label='vs1 change, one core')
    plt.plot(convs1[burn_in:],label='vs1 change, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vs1.pdf')
    
    file=open(directory+'/'+'Convergence_vs2.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    
    convs2_one=[]
    convs2=[]
    for line in lines[1:]:
        data=line.split()
        convs2_one.append(float(data[0]))
        convs2.append(float(data[1]))
    
    plt.figure('convergence_vs2')
    plt.title('vs2 change rate convergence')
    #plt.plot(convs2_one[burn_in:],label='vs2 change, one core')
    plt.plot(convs2[burn_in:],label='vs2 change, all cores')
    plt.ylim([0.,100.])
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vs2.pdf')

###################################################
# Convergence during preparation
###################################################

if preprocessing:
    file=open(directory+'/'+'Convergence_misfit_prep.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    n_w=int(lines[0].split()[1])
    burn_in_widening=int(lines[0].split()[2])
    nsample_widening=int(lines[0].split()[3])
    
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
    
    plt.figure('convergence_misfit_prep')
    plt.title('misfit convergence')
    #plt.plot(conv_R_one,label='Rayleigh, one core')
    plt.plot(conv_R,label='Rayleigh, all cores')
    #plt.plot(conv_L_one,label='Love, one core')
    plt.plot(conv_L,label='Love, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_misfit_prep.pdf')
    
    
    file=open(directory+'/'+'Convergence_nb_layers_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    conv_n_one=[]
    conv_n=[]
    for line in lines[1:]:
        data=line.split()
        conv_n_one.append(float(data[0]))
        conv_n.append(float(data[1]))
    
    plt.figure('convergence_nlayers_prep')
    plt.title('number of layers convergence')
    #plt.plot(conv_n_one,label='nblayers, one core')
    plt.plot(conv_n,label='nblayers, all cores',linewidth=0.1)
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_np_layers_prep.pdf')
    
    file=open(directory+'/'+'Convergence_sigma_R_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    conv_sigmaR_one=[]
    conv_sigmaR=[]
    for line in lines[1:]:
        data=line.split()
        conv_sigmaR_one.append(float(data[0]))
        conv_sigmaR.append(float(data[1]))
    
    plt.figure('convergence_sigmaR_prep')
    plt.title('sigmaR convergence')
    #plt.plot(conv_sigmaR_one,label='sigmaR, one core')
    plt.plot(conv_sigmaR,label='sigmaR, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_sigmaR_prep.pdf')
    
    file=open(directory+'/'+'Convergence_sigma_L_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    conv_sigmaL_one=[]
    conv_sigmaL=[]
    for line in lines[1:]:
        data=line.split()
        conv_sigmaL_one.append(float(data[0]))
        conv_sigmaL.append(float(data[1]))
    
    plt.figure('convergence_sigmaL_prep')
    plt.title('sigmaL convergence')
    #plt.plot(conv_sigmaL_one,label='sigmaL, one core')
    plt.plot(conv_sigmaL,label='sigmaL, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    
    plt.savefig(directory+'/PLOTS/'+'convergence_sigmaL_prep.pdf')
    
    file=open(directory+'/'+'Convergence_Birth_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
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
    
    plt.figure('convergence_birth_prep')
    plt.title('birth rate convergence')
    #plt.plot(convB_one,label='birth, one core')
    plt.plot(convB,label='birth, all cores')
    #plt.plot(convBa_one,label='birth anisotropic, one core')
    plt.plot(convBa,label='birth anisotropic, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_birth_prep.pdf')
    
    file=open(directory+'/'+'Convergence_Death_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
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
    
    plt.figure('convergence_death_prep')
    plt.title('death rate convergence')
    #plt.plot(convD_one,label='death, one core')
    plt.plot(convD,label='death, all cores')
    #plt.plot(convDa_one,label='death anisotropic, one core')
    plt.plot(convDa,label='death anisotropic, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_death_prep.pdf')
    
    file=open(directory+'/'+'Convergence_xi_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    convP_one=[]
    convP=[]
    for line in lines[1:]:
        data=line.split()
        convP_one.append(float(data[0]))
        convP.append(float(data[1]))
    
    plt.figure('convergence_xi_prep')
    plt.title('anisotropy change rate convergence')
    #plt.plot(convP_one,label='xi change, one core')
    plt.plot(convP,label='xi change, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_xi_prep.pdf')
    
    file=open(directory+'/'+'Convergence_vp_prep.out','r')
    lines=file.readlines()
    file.close()
    
    burn_in=int(lines[0].split()[0])
    n_w=int(lines[0].split()[1])
    burn_in_widening=int(lines[0].split()[2])
    nsample_widening=int(lines[0].split()[3])
    
    convP_one=[]
    convP=[]
    for line in lines[1:]:
        data=line.split()
        convP_one.append(float(data[0]))
        convP.append(float(data[1]))
    
    plt.figure('convergence_vp_prep')
    plt.title('vp change rate convergence')
    #plt.plot(convP_one,label='vp change, one core')
    plt.plot(convP,label='vp change, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vp_prep.pdf')
    
    file=open(directory+'/'+'Convergence_vs1_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    convs1_one=[]
    convs1=[]
    for line in lines[1:]:
        data=line.split()
        convs1_one.append(float(data[0]))
        convs1.append(float(data[1]))
    
    plt.figure('convergence_vs1_prep')
    plt.title('vs1 change rate convergence')
    #plt.plot(convs1_one,label='vs1 change, one core')
    plt.plot(convs1,label='vs1 change, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vs1_prep.pdf')
    
    file=open(directory+'/'+'Convergence_vs2_prep.out','r')
    lines=file.readlines()
    file.close()
    
    try:
        burn_in=int(lines[0].split()[0])
        n_w=int(lines[0].split()[1])
        burn_in_widening=int(lines[0].split()[2])
        nsample_widening=int(lines[0].split()[3])
    except:
        pass
    
    convs2_one=[]
    convs2=[]
    for line in lines[1:]:
        data=line.split()
        convs2_one.append(float(data[0]))
        convs2.append(float(data[1]))
    
    plt.figure('convergence_vs2_prep')
    plt.title('vs2 change rate convergence')
    #plt.plot(convs2_one,label='vs2 change, one core')
    plt.plot(convs2,label='vs2 change, all cores')
    plt.axvline(burn_in,c='k')
    for i in range(n_w):
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening),c='b')
        plt.axvline(burn_in+i*(burn_in_widening+nsample_widening)+burn_in_widening,c='r')
    plt.legend()
    plt.ylim([0.,100.])
    
    plt.savefig(directory+'/PLOTS/'+'convergence_vs2_prep.pdf')
    
###################################################
# Dispersion curves
###################################################

if results and dispersion:
    # reference
    
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
    period_R=np.array(period_R)
    n_R=np.array(n_R)
    c_R=np.array(c_R)
    dc_R=np.array(dc_R)
        
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
    period_L=np.array(period_L)
    n_L=np.array(n_L)
    c_L=np.array(c_L)
    dc_L=np.array(dc_L)
        
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
    period_R_obs=np.array(period_R_obs)
    n_R_obs=np.array(n_R_obs)
    c_R_obs=np.array(c_R_obs)
    dc_R_obs=np.array(dc_R_obs)
        
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
    period_L_obs=np.array(period_L_obs)
    n_L_obs=np.array(n_L_obs)
    c_L_obs=np.array(c_L_obs)
    dc_L_obs=np.array(dc_L_obs)
    
    plt.figure('dispersion')
    # plt.errorbar(period_R,c_R,yerr=dc_R,marker='o',zorder=0,label='Rayleigh average')
    # plt.errorbar(np.array(period_L)+0.1,c_L,yerr=dc_L,marker='o',zorder=0,label='Love average')
    # plt.errorbar(np.array(period_R_obs)-0.1,c_R_obs,yerr=dc_R_obs,marker='o',zorder=0,label='Rayleigh observed')
    # plt.errorbar(np.array(period_L_obs)+0.2,c_L_obs,yerr=dc_L_obs,marker='o',zorder=0,label='Love observed')
    for j in np.unique(n_R_obs):
        ints=np.where(n_R_obs==j)
        if j==0:
            plt.errorbar(period_R[ints[0]],c_R[ints[0]],yerr=dc_R[ints[0]],marker='o',zorder=0,label='Rayleigh average',c='C0')
            plt.errorbar(period_R_obs[ints[0]]-0.1,c_R_obs[ints[0]],yerr=dc_R_obs[ints[0]],marker='o',zorder=0,label='Rayleigh observed',c='C1')
        else:
            plt.errorbar(period_R[ints[0]],c_R[ints[0]],yerr=dc_R[ints[0]],marker='o',zorder=0,c='C0')
            plt.errorbar(period_R_obs[ints[0]]-0.1,c_R_obs[ints[0]],yerr=dc_R_obs[ints[0]],marker='o',zorder=0,c='C1')
    for j in np.unique(n_L_obs):
        ints=np.where(n_L_obs==j)
        if j==0:
            plt.errorbar(period_L[ints[0]]+0.1,c_L[ints[0]],yerr=dc_L[ints[0]],marker='o',zorder=0,label='Love average',c='C2')
            plt.errorbar(period_L_obs[ints[0]]+0.2,c_L_obs[ints[0]],yerr=dc_L_obs[ints[0]],marker='o',zorder=0,label='Love observed',c='C3')
        else:
            plt.errorbar(period_L[ints[0]]+0.1,c_L[ints[0]],yerr=dc_L[ints[0]],marker='o',zorder=0,c='C2')
            plt.errorbar(period_L_obs[ints[0]]+0.2,c_L_obs[ints[0]],yerr=dc_L_obs[ints[0]],marker='o',zorder=0,c='C3')
    plt.legend()
    plt.xlabel('period, s')
    plt.ylabel("phase velocity, km/s")
    plt.title('compared dispersion curves')
    #plt.show()
    
    plt.savefig(directory+'/PLOTS/'+'dispersion_ref.pdf')
    
    # widened reference
    
    file=open(directory+'/'+'Dispersion_mean_wide.out','r')
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
    period_R=np.array(period_R)
    n_R=np.array(n_R)
    c_R=np.array(c_R)
    dc_R=np.array(dc_R)
        
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
    period_L=np.array(period_L)
    n_L=np.array(n_L)
    c_L=np.array(c_L)
    dc_L=np.array(dc_L)
        
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
    period_R_obs=np.array(period_R_obs)
    n_R_obs=np.array(n_R_obs)
    c_R_obs=np.array(c_R_obs)
    dc_R_obs=np.array(dc_R_obs)
        
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
    period_L_obs=np.array(period_L_obs)
    n_L_obs=np.array(n_L_obs)
    c_L_obs=np.array(c_L_obs)
    dc_L_obs=np.array(dc_L_obs)
    
    plt.figure('dispersion widened')
    # plt.errorbar(period_R,c_R,yerr=dc_R,marker='o',zorder=0,label='Rayleigh average')
    # plt.errorbar(np.array(period_L)+0.1,c_L,yerr=dc_L,marker='o',zorder=0,label='Love average')
    # plt.errorbar(np.array(period_R_obs)-0.1,c_R_obs,yerr=dc_R_obs,marker='o',zorder=0,label='Rayleigh observed')
    # plt.errorbar(np.array(period_L_obs)+0.2,c_L_obs,yerr=dc_L_obs,marker='o',zorder=0,label='Love observed')
    for j in np.unique(n_R_obs):
        ints=np.where(n_R_obs==j)
        if j==0:
            plt.errorbar(period_R[ints[0]],c_R[ints[0]],yerr=dc_R[ints[0]],marker='o',zorder=0,label='Rayleigh average',c='C0')
            plt.errorbar(period_R_obs[ints[0]]-0.1,c_R_obs[ints[0]],yerr=dc_R_obs[ints[0]],marker='o',zorder=0,label='Rayleigh observed',c='C1')
        else:
            plt.errorbar(period_R[ints[0]],c_R[ints[0]],yerr=dc_R[ints[0]],marker='o',zorder=0,c='C0')
            plt.errorbar(period_R_obs[ints[0]]-0.1,c_R_obs[ints[0]],yerr=dc_R_obs[ints[0]],marker='o',zorder=0,c='C1')
    for j in np.unique(n_L_obs):
        ints=np.where(n_L_obs==j)
        if j==0:
            plt.errorbar(period_L[ints[0]]+0.1,c_L[ints[0]],yerr=dc_L[ints[0]],marker='o',zorder=0,label='Love average',c='C2')
            plt.errorbar(period_L_obs[ints[0]]+0.2,c_L_obs[ints[0]],yerr=dc_L_obs[ints[0]],marker='o',zorder=0,label='Love observed',c='C3')
        else:
            plt.errorbar(period_L[ints[0]]+0.1,c_L[ints[0]],yerr=dc_L[ints[0]],marker='o',zorder=0,c='C2')
            plt.errorbar(period_L_obs[ints[0]]+0.2,c_L_obs[ints[0]],yerr=dc_L_obs[ints[0]],marker='o',zorder=0,c='C3')
    
    plt.legend()
    plt.xlabel('period, s')
    plt.ylabel("phase velocity, km/s")
    plt.title('compared dispersion curves, widened data')
    #plt.show()
    
    plt.savefig(directory+'/PLOTS/'+'dispersion_wide.pdf')
    
    # cluster
    if many_plots:
        file=open(directory+'/'+'Dispersion_mean_alt.out','r')
        lines=file.readlines()
        file.close()
        
        numdis=int(float(lines[0]))
        ndatad_R=int(lines[1].split()[0])
        ndatad_L=int(lines[1].split()[1])
        
        period_R=np.zeros((ndatad_R,numdis))
        n_R=np.zeros((ndatad_R,numdis))
        c_R=np.zeros((ndatad_R,numdis))
        dc_R=np.zeros((ndatad_R,numdis))
        
            
        period_L=np.zeros((ndatad_L,numdis))
        n_L=np.zeros((ndatad_L,numdis))
        c_L=np.zeros((ndatad_L,numdis))
        dc_L=np.zeros((ndatad_L,numdis))
        
        k=2
        for i in range(numdis):
            j=0
            while (j<ndatad_R):
                line=lines[k]
                data=line.split()
                period_R[j,i]=float(data[0])
                n_R[j,i]=int(float(data[1]))
                c_R[j,i]=float(data[2])
                dc_R[j,i]=float(data[3])
                j+=1
                k+=1
            j=0
            while (j<ndatad_L):
                line=lines[k]
                data=line.split()
                period_L[j,i]=float(data[0])
                n_L[j,i]=int(float(data[1]))
                c_L[j,i]=float(data[2])
                dc_L[j,i]=float(data[3])
                j+=1
                k+=1
        
        file=open(directory+'/'+'Dispersion_obs_alt.out','r')
        lines=file.readlines()
        file.close()
        
        numdis=int(float(lines[0]))
        ndatad_R=int(lines[1].split()[0])
        ndatad_L=int(lines[1].split()[1])
        
        period_R_obs=np.zeros((ndatad_R,numdis))
        n_R_obs=np.zeros((ndatad_R,numdis))
        c_R_obs=np.zeros((ndatad_R,numdis))
        dc_R_obs=np.zeros((ndatad_R,numdis))
        
            
        period_L_obs=np.zeros((ndatad_L,numdis))
        n_L_obs=np.zeros((ndatad_L,numdis))
        c_L_obs=np.zeros((ndatad_L,numdis))
        dc_L_obs=np.zeros((ndatad_L,numdis))
        
        k=2
        for i in range(numdis):
            j=0
            while (j<ndatad_R):
                line=lines[k]
                data=line.split()
                period_R_obs[j,i]=float(data[0])
                n_R_obs[j,i]=int(float(data[1]))
                c_R_obs[j,i]=float(data[2])
                dc_R_obs[j,i]=float(data[3])
                j+=1
                k+=1
            j=0
            while (j<ndatad_L):
                line=lines[k]
                data=line.split()
                period_L_obs[j,i]=float(data[0])
                n_L_obs[j,i]=int(float(data[1]))
                c_L_obs[j,i]=float(data[2])
                dc_L_obs[j,i]=float(data[3])
                j+=1
                k+=1
        
        for i in range(numdis):
            #with PdfPages(directory+'/PLOTS/dispersions.pdf') as pdf:
            fig=plt.figure('dispersion '+str(i+1))
            for j in np.unique(n_R_obs):
                ints=np.where(n_R_obs==j)
                if j==0:
                    plt.errorbar(period_R[ints[0],i],c_R[ints[0],i],yerr=dc_R[ints[0],i],marker='o',zorder=0,label='Rayleigh average '+str(i+1),c='C0')
                    plt.errorbar(period_R_obs[ints[0],i]-0.1,c_R_obs[ints[0],i],yerr=dc_R_obs[ints[0],i],marker='o',zorder=0,label='Rayleigh observed '+str(i+1),c='C1')
                else:
                    plt.errorbar(period_R[ints[0],i],c_R[ints[0],i],yerr=dc_R[ints[0],i],marker='o',zorder=0,c='C0')
                    plt.errorbar(period_R_obs[ints[0],i]-0.1,c_R_obs[ints[0],i],yerr=dc_R_obs[ints[0],i],marker='o',zorder=0,c='C1')
            for j in np.unique(n_L_obs):
                ints=np.where(n_L_obs==j)
                if j==0:
                    plt.errorbar(period_L[ints[0],i]+0.1,c_L[ints[0],i],yerr=dc_L[ints[0],i],marker='o',zorder=0,label='Love average '+str(i+1),c='C2')
                
                    plt.errorbar(period_L_obs[ints[0],i]+0.2,c_L_obs[ints[0],i],yerr=dc_L_obs[ints[0],i],marker='o',zorder=0,label='Love observed '+str(i+1),c='C3')
                else:
                    plt.errorbar(period_L[ints[0],i]+0.1,c_L[ints[0],i],yerr=dc_L[ints[0],i],marker='o',zorder=0,c='C2')
                
                    plt.errorbar(period_L_obs[ints[0],i]+0.2,c_L_obs[ints[0],i],yerr=dc_L_obs[ints[0],i],marker='o',zorder=0,c='C3')
            plt.legend()
            plt.xlabel('period, s')
            plt.ylabel("phase velocity, km/s")
            plt.title('compared dispersion curves, model '+str(i+1))
            
            #pdf.savefig()
            plt.savefig(directory+'/PLOTS/Dispersion/'+'dispersion_'+str(i)+'.pdf')
        #plt.show()
            if not os.path.exists(directory+'/PLOTS/Individuals/'+str(i+1)):
                os.mkdir(directory+'/PLOTS/Individuals/'+str(i+1))
            plt.savefig(directory+'/PLOTS/'+'/Individuals/'+str(i+1)+'/dispersion.pdf')
            
            plt.close(fig)
        del lines

##################################################
# widening test results
##################################################

if preprocessing:
    file=open(directory+'/'+'mean_prop.out','r')
    lines=file.readlines()
    file.close()
    
    data=lines[0].split()
    widening_start=round(float(data[0]),3)
    widening_step=round(float(data[1]),3)
    n_w=int(float(data[2]))
    numdis=int(float(data[3]))
    
    means=np.zeros((n_w,numdis))
    widenings=np.arange(widening_start,widening_start+(n_w)*widening_step,step=widening_step)
    
    ind=0
    for i,line in enumerate(lines[1:]):
        #print(len(line.split()))
        #print(numdis)
        if ind==1:
            data=line.split()
            for j in range(numdis):

                means[(i-1)//4,j]=float(data[j])
        ind+=1
        if ind==4:
            ind=0
    
    # print(widenings)
    # print(widening_step)
    # print(n_w)
    # print(len(means[:,0]))
    plt.figure('widening results')
    plt.xlabel('widening')
    plt.ylabel('distance to widened reference')
    for i in range(numdis):
        plt.plot(widenings,means[:,i],label=str(i+1))
    
    #plt.legend()
    plt.title('distance change with widening')
    
    plt.savefig(directory+'/PLOTS/'+'mean_prop.pdf')
    
    file=open(directory+'/'+'alphahist.out','r')
    lines=file.readlines()
    file.close()
    
    data=lines[0].split()
    logalpha_min=float(data[0])
    logalpha_max=float(data[1])
    num_logalpha=int(float(data[2]))
    
    data=lines[1].split()
    widening_start=round(float(data[0]),3)
    widening_step=round(float(data[1]),3)
    n_w=int(float(data[2]))
    numdis=int(float(data[3]))
    widenings=np.arange(widening_start,widening_start+(n_w)*widening_step,step=widening_step)
    
    alphas=np.zeros((num_logalpha,n_w,numdis))
    alpha_range=np.linspace(logalpha_min,logalpha_max,num=num_logalpha)
    
    num_wide=0
    num_dis=0
    for i,line in enumerate(lines[2:]):
        for j,data in enumerate(line.split()):
            alphas[j,num_wide,num_dis]=float(data)
        num_dis+=1
        if num_dis==numdis:
            num_dis=0
            num_wide+=1
    
    
    file=open(directory+'/'+'mean_prop.out','r')
    lines=file.readlines()
    file.close()
    
    data=lines[0].split()
    widening_start=round(float(data[0]),3)
    widening_step=round(float(data[1]),3)
    n_w=int(float(data[2]))
    numdis=int(float(data[3]))
    
    means=np.zeros((n_w,numdis))
    
    num_wide=0
    num_dis=0
    ind=0
    for i,line in enumerate(lines[1:]):
        #print(len(line.split()))
        #print(numdis)
        if ind==1:
            data=line.split()
            for j in range(numdis):
                #print(j)
                #print(numdis)
                
                means[(i-1)//4,j]=float(data[j])
        ind+=1
        if ind==4:
            ind=0
    
    for i in range(n_w):
        
        plt.figure('alpha histograms '+str(widenings[i]))
        plt.xlabel('alpha')
        for j in range(numdis):
            plt.plot(alpha_range[1:],alphas[1:,i,j],label=str(j+1))
            #plt.axvline(means[i,j])
        
    
        #plt.legend()
        plt.title('alpha histograms '+str(widenings[i]))
        
        plt.savefig(directory+'/PLOTS/'+'alpha_hist_prep_'+str(i)+'.pdf')
        
    plt.figure('mean curves')
    plt.plot(widenings,np.mean(means,axis=1))
    plt.plot(widenings,np.median(means,axis=1))
    plt.plot(widenings,np.amin(means,axis=1))
    plt.savefig(directory+'/PLOTS/'+'mean_median_min_prep.pdf')
    
    means_upper_half=[]
    medians_upper_half=[]
    mins_upper_half=[]
    for i in range(n_w):
        means_upper_half.append(np.mean(means[i,np.where(means[i,:]>np.median(means[i,:]))]))
        medians_upper_half.append(np.median(means[i,np.where(means[i,:]>np.median(means[i,:]))]))
        mins_upper_half.append(np.amin(means[i,np.where(means[i,:]>np.median(means[i,:]))]))
    plt.figure('upper half mean curves')
    plt.plot(widenings,means_upper_half)
    plt.plot(widenings,medians_upper_half)
    plt.plot(widenings,mins_upper_half)
    plt.savefig(directory+'/PLOTS/'+'mean_median_min_prep_upper_half.pdf')
#plt.close('all')