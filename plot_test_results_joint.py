#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:06:13 2021

@author: dorian
"""
import matplotlib.pyplot as plt
import numpy as np

directory='OUT_TESTJOINT2'


##################
# Alpha Histogram 
##################

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

plt.legend()

##########################################
# Number of Layers histogram
##########################################

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

################################################
# Tradeoffs
################################################

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

# Tradeoff for cluster members

f=open(directory+'/'+'Tradeoff_alt.out','r')
lines=f.readlines()
f.close()

nlay=int(float(lines[0].split()[0]))
numdis=int(float(lines[0].split()[1]))
l=np.zeros((nlay+1,nlay,numdis))

i=0
j=0
for line in lines[1:]:
    l[j,i]=float(line.split()[0])#/(i+1)
    j+=1
    if j==nlay+1:
        i+=1
        j=0

for i in range(numdis):
    plt.figure('TRA_wide '+ str(i+1))
    plt.title('Tradeoff cluster '+ str(i+1))
    plt.pcolor(l[:,:,i])
    plt.colorbar()

#################################################
# Posterior
#################################################

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


ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
try:
    ax1.plot(true_xi,true_depth,c='c',linewidth=5)
    ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
except: 
    pass

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

ax0.plot(average_vs,depths,c='r',linewidth=5)
ax1.plot(average_xi,depths,c='r',linewidth=5)
ax2.plot(average_vpvs,depths,c='r',linewidth=5)
ax4.plot(average_probani,depths,c='k',linewidth=5)
ax4.set_xlabel('anisotropy probability',fontsize=15)
ax4.set_xlim([0,100])

fig.suptitle('posterior and averages - reference')

# widened Reference 

file=open(directory+'/'+'Posterior_wide.out','r')
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

depth=[]

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


ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
try:
    ax1.plot(true_xi,true_depth,c='c',linewidth=5)
    ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
except: 
    pass

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

ax0.plot(average_vs,depths,c='r',linewidth=5)
ax1.plot(average_xi,depths,c='r',linewidth=5)
ax2.plot(average_vpvs,depths,c='r',linewidth=5)
ax4.plot(average_probani,depths,c='k',linewidth=5)
ax4.set_xlabel('anisotropy probability',fontsize=15)
ax4.set_xlim([0,100])

fig.suptitle('posterior and averages - widened reference')


# cluster members 

file=open(directory+'/'+'Posterior_alt.out','r')
depthplot=60

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

depth=[]

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
    for j,data2 in enumerate(data[1:]):
        change_hist[i,j]=float(data2)

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

    file=open(directory+'/'+'true_model_'+str(k+1)+'.out','r')
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
    
    
    
    
    
    ax0.plot(true_vsv,true_depth,c='c',linewidth=5)
    try:
        ax1.plot(true_xi,true_depth,c='c',linewidth=5)
        ax2.plot(true_vpvs,true_depth,c='c',linewidth=5)
    except: 
        pass

    ax5.plot(change_hist[:,k],change_depths)
    ax5.set_xlabel('change point histogram',fontsize=15)
    
    ax0.plot(average_vs[:,k],depths,c='r',linewidth=5)
    ax1.plot(average_xi[:,k],depths,c='r',linewidth=5)
    ax2.plot(average_vpvs[:,k],depths,c='r',linewidth=5)
    ax4.plot(average_probani[:,k],depths,c='k',linewidth=5)
    ax4.set_xlabel('anisotropy probability',fontsize=15)
    ax4.set_xlim([0,100])
    
    plt.setp(ax2.get_yticklabels(), visible=False)
    
    fig.suptitle('posterior and averages '+str(k+1))

#################################################
# Sigma R
#################################################

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

plt.figure('sigmad_L cluster')

for i in range(numdis):
    
    plt.title('sigmad_L cluster')
    plt.plot(dL,sigmad_L[:,i],label=str(i+1))

plt.legend()

#########################################################
# Convergences
########################################################

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
plt.plot(conv_R_one[burn_in:],label='Rayleigh, one core')
plt.plot(conv_R[burn_in:],label='Rayleigh, all cores')
plt.plot(conv_L_one[burn_in:],label='Love, one core')
plt.plot(conv_L[burn_in:],label='Love, all cores')
plt.legend()


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
plt.plot(conv_n_one[burn_in:],label='nblayers, one core')
plt.plot(conv_n[burn_in:],label='nblayers, all cores')
plt.legend()

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
plt.plot(conv_sigmaR_one[burn_in:],label='sigmaR, one core')
plt.plot(conv_sigmaR[burn_in:],label='sigmaR, all cores')
plt.legend()

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
plt.plot(conv_sigmaL_one[burn_in:],label='sigmaL, one core')
plt.plot(conv_sigmaL[burn_in:],label='sigmaL, all cores')
plt.legend()

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
plt.plot(convB_one[burn_in:],label='birth, one core')
plt.plot(convB[burn_in:],label='birth, all cores')
plt.plot(convBa_one[burn_in:],label='birth anisotropic, one core')
plt.plot(convBa[burn_in:],label='birth anisotropic, all cores')
plt.legend()

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
plt.plot(convD_one[burn_in:],label='death, one core')
plt.plot(convD[burn_in:],label='death, all cores')
plt.plot(convDa_one[burn_in:],label='death anisotropic, one core')
plt.plot(convDa[burn_in:],label='death anisotropic, all cores')
plt.legend()

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
plt.plot(convP_one[burn_in:],label='xi change, one core')
plt.plot(convP[burn_in:],label='xi change, all cores')
plt.legend()

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
plt.plot(convP_one[burn_in:],label='vp change, one core')
plt.plot(convP[burn_in:],label='vp change, all cores')
plt.legend()

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
plt.plot(convs1_one[burn_in:],label='vs1 change, one core')
plt.plot(convs1[burn_in:],label='vs1 change, all cores')
plt.legend()

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
plt.plot(convs2_one[burn_in:],label='vs2 change, one core')
plt.plot(convs2[burn_in:],label='vs2 change, all cores')
plt.legend()

###################################################
# Dispersion curves
###################################################

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

plt.figure('dispersion')
plt.errorbar(period_R,c_R,yerr=dc_R,marker='o',zorder=0,label='Rayleigh average')
plt.errorbar(np.array(period_L)+0.1,c_L,yerr=dc_L,marker='o',zorder=0,label='Love average')
plt.errorbar(np.array(period_R_obs)-0.1,c_R_obs,yerr=dc_R_obs,marker='o',zorder=0,label='Rayleigh observed')
plt.errorbar(np.array(period_L_obs)+0.2,c_L_obs,yerr=dc_L_obs,marker='o',zorder=0,label='Love observed')
plt.legend()
plt.xlabel('period, s')
plt.ylabel("phase velocity, km/s")
plt.title('compared dispersion curves')
#plt.show()

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

plt.figure('dispersion widened')
plt.errorbar(period_R,c_R,yerr=dc_R,marker='o',zorder=0,label='Rayleigh average')
plt.errorbar(np.array(period_L)+0.1,c_L,yerr=dc_L,marker='o',zorder=0,label='Love average')
plt.errorbar(np.array(period_R_obs)-0.1,c_R_obs,yerr=dc_R_obs,marker='o',zorder=0,label='Rayleigh observed')
plt.errorbar(np.array(period_L_obs)+0.2,c_L_obs,yerr=dc_L_obs,marker='o',zorder=0,label='Love observed')
plt.legend()
plt.xlabel('period, s')
plt.ylabel("phase velocity, km/s")
plt.title('compared dispersion curves, widened data')
#plt.show()

# cluster

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

    plt.figure('dispersion '+str(i+1))
    plt.errorbar(period_R[:,i],c_R[:,i],yerr=dc_R[:,i],marker='o',zorder=0,label='Rayleigh average '+str(i+1))
    plt.errorbar(period_L[:,i]+0.1,c_L[:,i],yerr=dc_L[:,i],marker='o',zorder=0,label='Love average '+str(i+1))
    plt.errorbar(period_R_obs[:,i]-0.1,c_R_obs[:,i],yerr=dc_R_obs[:,i],marker='o',zorder=0,label='Rayleigh observed '+str(i+1))
    plt.errorbar(period_L_obs[:,i]+0.2,c_L_obs[:,i],yerr=dc_L_obs[:,i],marker='o',zorder=0,label='Love observed '+str(i+1))
    plt.legend()
    plt.xlabel('period, s')
    plt.ylabel("phase velocity, km/s")
    plt.title('compared dispersion curves, model '+str(i+1))
#plt.show()

##################################################
# widening test results
##################################################

file=open(directory+'/'+'mean_prop.out','r')
lines=file.readlines()
file.close()

data=lines[0].split()
widening_start=float(data[0])
widening_step=float(data[1])
n_w=int(float(data[2]))
numdis=int(float(data[3]))

means=np.zeros((n_w,numdis))
widenings=np.arange(widening_start,widening_start+(n_w)*widening_step,step=widening_step)

for i,line in enumerate(lines[1:]):
    for j,data in enumerate(line.split()):
        means[i,j]=float(data)

plt.figure('widening results')
plt.xlabel('widening')
plt.ylabel('distance to widened reference')
for i in range(numdis):
    plt.plot(widenings,means[:,i],label=str(i+1))

plt.legend()
plt.title('distance change with widening')

file=open(directory+'/'+'alphahist.out','r')
lines=file.readlines()
file.close()

data=lines[0].split()
logalpha_min=float(data[0])
logalpha_max=float(data[1])
num_logalpha=int(float(data[2]))

data=lines[1].split()
widening_start=float(data[0])
widening_step=float(data[1])
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
widening_start=float(data[0])
widening_step=float(data[1])
n_w=int(float(data[2]))
numdis=int(float(data[3]))

means=np.zeros((n_w,numdis))

num_wide=0
num_dis=0
for i,line in enumerate(lines[1:]):
    for j,data in enumerate(line.split()):
        means[i,j]=float(data)

for i in range(n_w):
    
    plt.figure('alpha histograms '+str(widenings[i]))
    plt.xlabel('alpha')
    for j in range(numdis):
        plt.plot(alpha_range,alphas[:,i,j],label=str(j+1))
        plt.axvline(means[i,j])

    plt.legend()
    plt.title('alpha histograms '+str(widenings[i]))