#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 14:06:13 2021

@author: dorian
"""
import matplotlib.pyplot as plt
import numpy as np
import os
import h5py
import sys
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

directory = 'OUT_ALL_DEEP_OCEAN_0.1'

results = True
posterior = True  # True
dispersion = True
tradeoff = False
alphas = True
sigma = True
many_plots = True
nlayers = True


def weighted_quantile(values, quantiles, sample_weight):
    datalen = np.shape(sample_weight)[0]
    cumsum = np.cumsum(sample_weight, axis=1)
    for line in range(datalen):
        cumsum[line, :] /= cumsum[line, -1]

    out = np.zeros((len(quantiles), datalen))

    for iquant, quant in enumerate(quantiles):
        ind_min = np.argmin(np.abs(cumsum - quant), axis=1)
        out[iquant, :] = values[ind_min]

    return out


if os.path.exists(directory + '/OUT/true_model.out'):
    print('true model available')
    true_model = True
else:
    true_model = False

file = directory + '/' + 'dispersion.h5'
f = h5py.File(file, 'r')
# numdis = int(np.array(f['lat'].get('numdis')))
# clusterlist = np.array(f['cluster_params'].get('clusterlist'))
# clusters={}
# for cl in clusterlist:
#     clusters[cl] = np.array(f['cluster_params']['cluster'].get(str(cl)))
batch = int(np.array(f['cluster_params'].get('batchsize')))
lat = np.array(f['cluster_params'].get('lat'))
lon = np.array(f['cluster_params'].get('lon'))
numdis = lat.size

points = []
lats_alistair = [63., -27., 65., -5., -29., 11., -17., 19., 31., 37.,17., 37., 45., 49., 53., 51.]
lons_alistair = [29., 27., -111., 33., 121., 41., -113., -155., 89., -111.,-133., -111., -97., -85., -71., -31.]
for j in range(len(lats_alistair)):
    print(lats_alistair[j], lons_alistair[j])
    print(lon)
    print(np.where(np.logical_and(lat == lats_alistair[j], lon == lons_alistair[j])))
    i = np.where(np.logical_and(lat == lats_alistair[j], lon == lons_alistair[j]))[0][0]
    print(i)
    points.append(i)

#points=np.arange(15250,15260).tolist()
#points=np.arange(numdis).tolist()
# points=list(np.where(np.logical_and(clusters==1,np.arange(numdis)<3000))[0])
# points = np.arange(numdis)
# points=np.arange(500)
num_files = np.amax(points) // batch + 1
# numdis=100
# points=np.arange(100)
# points=np.where(clusters==0)[0]
# print(points)

f.close()
# numdis=10

print(directory + '/PLOTS/')
# print(os.path.exists(directory+'/PLOTS/'))
# print(not os.path.exists(directory+'/PLOTS/'))
if not os.path.exists(directory + '/PLOTS/'):
    os.mkdir(directory + '/PLOTS/')

if not os.path.exists(directory + '/PLOTS/Posterior/'):
    os.mkdir(directory + '/PLOTS/Posterior/')

if not os.path.exists(directory + '/PLOTS/Tradeoff/'):
    os.mkdir(directory + '/PLOTS/Tradeoff/')

if not os.path.exists(directory + '/PLOTS/Dispersion/'):
    os.mkdir(directory + '/PLOTS/Dispersion/')

if not os.path.exists(directory + '/PLOTS/Alphas/'):
    os.mkdir(directory + '/PLOTS/Alphas/')

if not os.path.exists(directory + '/PLOTS/Nlayers/'):
    os.mkdir(directory + '/PLOTS/Nlayers/')

if not os.path.exists(directory + '/PLOTS/Nlayers/'):
    os.mkdir(directory + '/PLOTS/Nlayers/')

if not os.path.exists(directory + '/PLOTS/SigmaRs/'):
    os.mkdir(directory + '/PLOTS/SigmaRs/')

if not os.path.exists(directory + '/PLOTS/SigmaLs/'):
    os.mkdir(directory + '/PLOTS/SigmaLs/')

if not os.path.exists(directory + '/PLOTS/Individuals/'):
    print(directory + '/PLOTS/Individuals/')
    os.mkdir(directory + '/PLOTS/Individuals/')

maxnlay = 40

# widening=5.
# prepare=True
# maxpercent=0.

##################
# Alpha Histogram
##################

if results and alphas:

    print('plotting alphas')

    function = 'get_histograms'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)
            for j in points:
                i_file = j // batch
                if i_file == i:
                    points.remove(j)
    # f = h5py.File(file,'r')
    for i in range(num_files):
        if fs[i] is None:
            continue
        alpharange = np.array(fs[i][function]['nostack'].get('range_alpha'))
        break

    if many_plots:

        # numdis=np.shape(alphas)[0]

        for i in points:
            i_file = i // batch
            i_new = i % batch
            
            if fs[i_file]==None:
                continue
            
            alphas = np.array(fs[i_file][function]['stack'].get('alpha_hist'))
            expsum = np.array(fs[i_file][function]['stack'].get('exp_sum'))

            fig = plt.figure('alpha')
            j = np.argmin(np.abs(alpharange + 20))
            plt.title('alpha histogram ' + str(i + 1) + ', big enough: ' + str(
                np.sum(alphas[j:, i_new]) * expsum[i_new]))  # /np.sum(alphas[:,i])*100 +' %'
            plt.plot(alpharange[1:], alphas[1:, i_new] * expsum[i_new], label=float(i + 1))
            #plt.xlim([-5,5])
            #plt.ylim([0,20])
            plt.savefig(directory + '/PLOTS/' + 'Alphas/alpha_hist_' + str(i) + '.pdf')

            # for cl in clusters:
            #     if i in clusters[cl]:


            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(
            #             directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/alpha_hist.pdf')

            plt.close(fig)
    for i in range(num_files):
        try:
            fs[i].close()
        except:
            pass
#########################################
# Number of Layers histogram
##########################################

if results and nlayers:

    print('plotting histograms')

    function = 'get_histograms'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)
            for j in points:
                i_file = j // batch
                if i_file == i:
                    points.remove(j)

    for i in range(num_files):
        if fs[i] is None:
            continue
        # Reference

        range_nlay = np.array(fs[i][function]['nostack'].get('range_nlay'))

        # Cluster members

        nlay = np.array(fs[i][function]['stack'].get('nlay_hist'))

    # numdis=np.shape(nlay)[1]

    # plt.figure('nlayers_histogram all')
    # plt.title('nlayers histogram: Cluster members')
    # #plt.ylim([0,0.02])
    # for i in range(numdis):
    #     plt.plot(range_nlay,nlay[:,i],label=str(i+1))
    # #plt.legend()
    # plt.savefig(directory+'/PLOTS/'+'nb_layers_all.pdf')

    plt.figure('nlayers_histogram cluster')
    plt.title('nlayers histogram: Cluster centroid, ' + directory)
    plt.plot(range_nlay, nlay[:, -1], label='centroid')
    # plt.ylim([0,0.2])
    plt.savefig(directory + '/PLOTS/' + 'nb_layers_cluster.pdf')

    if many_plots:

        for i in points:
            i_file = i // batch
            i_new = i % batch
            if fs[i_file] is None:
                continue
            nlay = np.array(fs[i_file][function]['stack'].get('nlay_hist'))
            # print(i)
            fig = plt.figure('nlayers')
            plt.title('nlayers histogram ' + str(i + 1))
            plt.plot(range_nlay, nlay[:, i_new], label=float(i + 1))
            plt.savefig(directory + '/PLOTS/' + 'Nlayers/nlayers_hist_' + str(i) + '.pdf')

            # for cl in clusters:
            #     if i in clusters[cl]:
            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(
            #             directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/nlayers_hist.pdf')

            plt.close(fig)

    for i in range(num_files):
        if fs[i] is not None: fs[i].close()

################################################
# Tradeoffs
################################################

if results and tradeoff:
    # Tradeoff for reference

    print('plotting tradeoffs')

    function = 'get_tradeoff'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)
            for j in points:
                i_file = j // batch
                if i_file == i:
                    points.remove(j)

    # Tradeoff for cluster members

    for i in range(num_files):
        if fs[i] is None:
            continue
        tradeoff = np.array(fs[i][function]['stack'].get('tradeoff'))

    fig = plt.figure('TRA_wide')
    plt.title('Tradeoff cluster centroid, ' + directory)
    plt.pcolor(np.transpose(tradeoff[:, :, -1]))
    plt.colorbar()

    plt.savefig(directory + '/PLOTS/' + 'tradeoff_cluster.pdf')
    plt.close(fig)

    if many_plots:

        for i in points:
            i_file = i // batch
            i_new = i % batch

            if fs[i_file] is None:
                continue
            tradeoff = np.array(fs[i_file][function]['stack'].get('tradeoff'))

            fig = plt.figure('TRA_wide')
            plt.title('Tradeoff cluster ' + str(i + 1))
            plt.pcolor(np.transpose(tradeoff[:, :, i_new]))
            plt.colorbar()

            plt.savefig(directory + '/PLOTS/Tradeoff/' + 'tradeoff_' + str(i) + '.pdf')

            # for cl in clusters:
            #     if i in clusters[cl]:
            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/tradeoff.pdf')

            plt.close(fig)
            # print(i)
    for i in range(num_files):
        if fs[i] is not None: fs[i].close()

#################################################
# Posterior
#################################################
if results and posterior:
    # Reference

    print('plotting posterior projections')

    function = 'create_posterior'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)

    for i in range(num_files):
        if fs[i] is None:
            continue
        depths = np.array(fs[i][function]['nostack'].get('depths'))
        vsvs = np.array(fs[i][function]['nostack'].get('vels_vsv'))
        xis = np.array(fs[i][function]['nostack'].get('vels_xi'))
        vps = np.array(fs[i][function]['nostack'].get('vels_vp'))
        vref_min = np.amin(vsvs)
        vref_max = np.amax(vsvs)
        prof = 600  # np.amax(depths)
        index = np.argmin(np.abs(depths - prof))
        xi_min = np.amin(xis)
        xi_max = np.amax(xis)
        vp_min = np.amin(vps)
        vp_max = np.amax(vps)

        vsvd_all = np.array(fs[i][function]['stack'].get('vsv_all'))
        xid_all = np.array(fs[i][function]['stack'].get('xi_all'))
        vpd_all = np.array(fs[i][function]['stack'].get('vp_all'))

        fs[i].close()

    function = 'get_average'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)
            for j in points:
                i_file = j // batch
                if i_file == i:
                    points.remove(j)
    for i in range(num_files):
        if fs[i] is None:
            continue
        depths_average = np.array(fs[i][function]['nostack'].get('depths'))
        vsv_av = np.array(fs[i][function]['stack'].get('vsv'))
        xi_av = np.array(fs[i][function]['stack'].get('xi'))
        vp_av = np.array(fs[i][function]['stack'].get('vp'))
        probani_all = np.array(fs[i][function]['stack'].get('probani'))

        fs[i].close()

    fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=1, ncols=4, sharey=True,
                                             figsize=(10, 6))

    # xid_all[:,np.shape(xid_all)[1]//2,-1]=0.

    ax0.invert_yaxis()
    ax0.set_xlim([vref_min, vref_max])
    ax0.set_xlabel('Vsv, km/s', fontsize=15)
    ax1.set_ylim([prof, 0.])
    ax1.set_xlabel(r'Xi', fontsize=15)
    ax1.set_xlim([xi_min, xi_max])
    ax2.set_xlabel(r'Vp, km/s', fontsize=15)
    ax2.set_xlim([vp_min, vp_max])
    ax0.set_ylabel('Depth, km', fontsize=15)
    ax2.set_xlabel(r'Vp, km/s', fontsize=15)
    ax0.pcolormesh(vsvs, depths[:-1], vsvd_all[:-1, :, -1], cmap='magma_r')
    ax1.pcolormesh(xis, depths[:-1], xid_all[:-1, :, -1], cmap='magma_r')
    ax2.pcolormesh(vps, depths[:-1], vpd_all[:-1, :, -1], cmap='magma_r')

    plt.setp(ax2.get_yticklabels(), visible=False)

    ax0.plot(vsv_av[:, -1], depths_average, c='r', linewidth=1)
    ax1.plot(xi_av[:, -1], depths_average, c='r', linewidth=1)
    ax2.plot(vp_av[:, -1], depths_average, c='r', linewidth=1)
    ax4.plot(probani_all[:, -1] * 100, depths_average, c='k', linewidth=1)

    ax4.set_xlabel('anisotropy probability', fontsize=15)
    ax4.set_xlim([0, 100])

    if true_model:
        depth_true = []
        vsv_true = []
        xi_true = []
        vpv_true = []
        f = open(directory + '/OUT/true_model.out')
        for line in f:
            data = line.split()
            depth_true.append(float(data[0]))
            vsv_true.append(float(data[1]) / 1000)
            xi_true.append(float(data[2]))
            vpv_true.append(float(data[3]) / 1000)

        f.close()
        ax0.plot(vsv_true, depth_true, c='r', linewidth=0.3)
        ax1.plot(xi_true, depth_true, c='r', linewidth=0.3)
        ax2.plot(vpv_true, depth_true, c='r', linewidth=0.3)

    fig.suptitle('posterior and averages centroid, ' + directory)

    plt.savefig(directory + '/PLOTS/Posterior/' + 'posterior_cluster.pdf')  #
    plt.close(fig)

    fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=1, ncols=4, sharey=True,
                                             figsize=(10, 6))

    # xid_all[:,np.shape(xid_all)[1]//2,-1]=0.

    ax0.invert_yaxis()
    ax0.set_xlim([vref_min, vref_max])
    ax0.set_xlabel('Vsv, km/s', fontsize=15)
    ax1.set_ylim([prof, 0.])
    ax1.set_xlabel(r'Xi', fontsize=15)
    ax1.set_xlim([xi_min, xi_max])
    ax2.set_xlabel(r'Vp, km/s', fontsize=15)
    ax2.set_xlim([vp_min, vp_max])
    ax0.set_ylabel('Depth, km', fontsize=15)
    ax2.set_xlabel(r'Vp, km/s', fontsize=15)

    percentiles = weighted_quantile(vsvs,
                                    [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2, 0.50 - 0.38 / 2,
                                     0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], vsvd_all[:-1, :, -1])
    #sys.exit()
    vsv_95min = percentiles[0, :]
    vsv_95max = percentiles[1, :]
    vsv_68min = percentiles[2, :]
    vsv_68max = percentiles[3, :]
    vsv_38min = percentiles[4, :]
    vsv_38max = percentiles[5, :]
    vsv_median = percentiles[6, :]
    vsv_10min = percentiles[7, :]
    vsv_10max = percentiles[8, :]

    ax0.fill_betweenx(depths[:-1], vsv_95min, vsv_95max, color='r', alpha=0.2,label='95%')
    ax0.fill_betweenx(depths[:-1], vsv_68min, vsv_68max, color='r', alpha=0.2,label='68%')
    ax0.fill_betweenx(depths[:-1], vsv_38min, vsv_38max, color='r', alpha=0.2,label='38%')
    ax0.fill_betweenx(depths[:-1], vsv_10min, vsv_10max, color='r', alpha=0.2,label='10%')
    ax0.plot(vsv_median, depths[:-1], c='r', linewidth=1,label='median')

    percentiles = weighted_quantile(xis,
                                    [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2, 0.50 - 0.38 / 2,
                                     0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], xid_all[:-1, :, -1])
    # sys.exit()
    xi_95min = percentiles[0, :]
    xi_95max = percentiles[1, :]
    xi_68min = percentiles[2, :]
    xi_68max = percentiles[3, :]
    xi_38min = percentiles[4, :]
    xi_38max = percentiles[5, :]
    xi_median = percentiles[6, :]
    xi_10min = percentiles[7, :]
    xi_10max = percentiles[8, :]

    ax1.fill_betweenx(depths[:-1], xi_95min, xi_95max, color='r', alpha=0.2)
    ax1.fill_betweenx(depths[:-1], xi_68min, xi_68max, color='r', alpha=0.2)
    ax1.fill_betweenx(depths[:-1], xi_38min, xi_38max, color='r', alpha=0.2)
    ax1.fill_betweenx(depths[:-1], xi_10min, xi_10max, color='r', alpha=0.2)
    ax1.plot(xi_median, depths[:-1], c='r', linewidth=1)

    percentiles = weighted_quantile(vps,
                                    [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2, 0.50 - 0.38 / 2,
                                     0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], vpd_all[:-1, :, -1])
    # sys.exit()
    vp_95min = percentiles[0, :]
    vp_95max = percentiles[1, :]
    vp_68min = percentiles[2, :]
    vp_68max = percentiles[3, :]
    vp_38min = percentiles[4, :]
    vp_38max = percentiles[5, :]
    vp_median = percentiles[6, :]
    vp_10min = percentiles[7, :]
    vp_10max = percentiles[8, :]

    ax2.fill_betweenx(depths[:-1], vp_95min, vp_95max, color='r', alpha=0.2)
    ax2.fill_betweenx(depths[:-1], vp_68min, vp_68max, color='r', alpha=0.2)
    ax2.fill_betweenx(depths[:-1], vp_38min, vp_38max, color='r', alpha=0.2)
    ax2.fill_betweenx(depths[:-1], vp_10min, vp_10max, color='r', alpha=0.2)
    ax2.plot(vp_median, depths[:-1], c='r', linewidth=1)

    plt.setp(ax2.get_yticklabels(), visible=False)

    ax0.plot(vsv_av[:, -1], depths_average, c='b', linewidth=1, label='mean')
    ax1.plot(xi_av[:, -1], depths_average, c='b', linewidth=1)
    ax2.plot(vp_av[:, -1], depths_average, c='b', linewidth=1)
    ax4.plot(probani_all[:, -1] * 100, depths_average, c='k', linewidth=1)

    ax4.set_xlabel('anisotropy probability', fontsize=15)
    ax4.set_xlim([0, 100])
    
    
    legend_elements = [Patch(facecolor='r', edgecolor='r', alpha=0.2,
                         label='95%'),Patch(facecolor='r', edgecolor='r', alpha=0.2*0.8+0.2,
                         label='68%'),Patch(facecolor='r', edgecolor='r', alpha=(0.2*0.8+0.2)*0.8+0.2,
                         label='38%'),Patch(facecolor='r', edgecolor='r', alpha=((0.2*0.8+0.2)*0.8+0.2)*0.8+0.2,
                         label='10%'),Line2D([0], [0], color='r', lw=1, label='median'),
                         Line2D([0], [0], color='b', lw=1, label='mean')]
    ax2.legend(handles=legend_elements)

    if true_model:
        depth_true = []
        vsv_true = []
        xi_true = []
        vpv_true = []
        f = open(directory + '/OUT/true_model.out')
        for line in f:
            data = line.split()
            depth_true.append(float(data[0]))
            vsv_true.append(float(data[1]) / 1000)
            xi_true.append(float(data[2]))
            vpv_true.append(float(data[3]) / 1000)

        f.close()
        ax0.plot(vsv_true, depth_true, c='r', linewidth=0.3)
        ax1.plot(xi_true, depth_true, c='r', linewidth=0.3)
        ax2.plot(vpv_true, depth_true, c='r', linewidth=0.3)

    fig.suptitle('posterior and averages centroid, ' + directory)

    plt.savefig(directory + '/PLOTS/Posterior/' + 'posterior_cluster_quantiles.pdf')  #
    plt.close(fig)

    if many_plots:

        function = 'create_posterior'
        fs = []
        for i in range(num_files):
            file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'

            if os.path.isfile(file):
                fs.append(h5py.File(file, 'r'))
            else:
                fs.append(None)
                for j in points:
                    i_file = j // batch
                    if i_file == i:
                        points.remove(j)
        function = 'get_average'
        fs2 = []
        for i in range(num_files):
            file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
            if os.path.isfile(file):
                fs2.append(h5py.File(file, 'r'))
            else:
                fs2.append(None)
                for j in points:
                    i_file = j // batch
                    if i_file == i:
                        points.remove(j)

        for i in points:
            i_file = i // batch
            i_new = i % batch
            if fs[i_file] is None:
                continue

            vsvd_all = np.array(fs[i_file]['create_posterior']['stack'].get('vsv_all'))
            xid_all = np.array(fs[i_file]['create_posterior']['stack'].get('xi_all'))
            vpd_all = np.array(fs[i_file]['create_posterior']['stack'].get('vp_all'))

            vsv_av = np.array(fs2[i_file]['get_average']['stack'].get('vsv'))
            xi_av = np.array(fs2[i_file]['get_average']['stack'].get('xi'))
            vp_av = np.array(fs2[i_file]['get_average']['stack'].get('vp'))
            probani_all = np.array(fs2[i_file]['get_average']['stack'].get('probani'))

            # print(i)
            fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=1, ncols=4, sharey=True,
                                                     figsize=(10, 6))

            ax0.invert_yaxis()
            ax0.set_xlim([vref_min, vref_max])
            ax0.set_xlabel('Vsv, km/s')
            ax1.set_ylim([prof, 0.])
            ax1.set_xlabel(r'Xi', fontsize=15)
            ax1.set_xlim([xi_min, xi_max])
            ax2.set_xlabel(r'Vp, km/s', fontsize=15)
            ax2.set_xlim([vp_min, vp_max])
            ax0.set_ylabel('Depth, km', fontsize=15)
            ax2.set_xlabel(r'Vp, km/s', fontsize=15)
            ax0.pcolormesh(vsvs, depths[:-1], vsvd_all[:-1, :, i_new], cmap='magma_r')
            ax1.pcolormesh(xis, depths[:-1], xid_all[:-1, :, i_new], cmap='magma_r')
            ax2.pcolormesh(vps, depths[:-1], vpd_all[:-1, :, i_new], cmap='magma_r')

            plt.setp(ax2.get_yticklabels(), visible=False)

            ax0.plot(vsv_av[:, i_new], depths_average, c='r', linewidth=1)
            ax1.plot(xi_av[:, i_new], depths_average, c='r', linewidth=1)
            ax2.plot(vp_av[:, i_new], depths_average, c='r', linewidth=1)
            ax4.plot(probani_all[:, i_new] * 100, depths_average, c='k', linewidth=1)

            ax4.set_xlabel('anisotropy probability', fontsize=15)
            ax4.set_xlim([0, 100])

            fig.suptitle('posterior and averages ' + str(i + 1) + ', lat: ' + str(lat[i]) + ', lon: ' + str(lon[i]))

            plt.savefig(directory + '/PLOTS/Posterior/' + 'posterior_' + str(i) + '.pdf')  #

            # for cl in clusters:
            #     if i in clusters[cl]:
            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/posterior.pdf')

            plt.close(fig)

            fig, (ax0, ax1, ax2, ax4) = plt.subplots(nrows=1, ncols=4, sharey=True,
                                                     figsize=(10, 6))

            # xid_all[:,np.shape(xid_all)[1]//2,-1]=0.

            ax0.invert_yaxis()
            ax0.set_xlim([vref_min, vref_max])
            ax0.set_xlabel('Vsv, km/s', fontsize=15)
            ax1.set_ylim([prof, 0.])
            ax1.set_xlabel(r'Xi', fontsize=15)
            ax1.set_xlim([xi_min, xi_max])
            ax2.set_xlabel(r'Vp, km/s', fontsize=15)
            ax2.set_xlim([vp_min, vp_max])
            ax0.set_ylabel('Depth, km', fontsize=15)
            ax2.set_xlabel(r'Vp, km/s', fontsize=15)

            percentiles = weighted_quantile(vsvs,
                                            [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2,
                                             0.50 - 0.38 / 2,
                                             0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], vsvd_all[:-1, :, i_new])
            # sys.exit()
            vsv_95min = percentiles[0, :]
            vsv_95max = percentiles[1, :]
            vsv_68min = percentiles[2, :]
            vsv_68max = percentiles[3, :]
            vsv_38min = percentiles[4, :]
            vsv_38max = percentiles[5, :]
            vsv_median = percentiles[6, :]
            vsv_10min = percentiles[7, :]
            vsv_10max = percentiles[8, :]

            ax0.fill_betweenx(depths[:-1], vsv_95min, vsv_95max, color='r', alpha=0.2)
            ax0.fill_betweenx(depths[:-1], vsv_68min, vsv_68max, color='r', alpha=0.2)
            ax0.fill_betweenx(depths[:-1], vsv_38min, vsv_38max, color='r', alpha=0.2)
            ax0.fill_betweenx(depths[:-1], vsv_10min, vsv_10max, color='r', alpha=0.2)
            ax0.plot(vsv_median, depths[:-1], c='r', linewidth=1)

            percentiles = weighted_quantile(xis,
                                            [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2,
                                             0.50 - 0.38 / 2,
                                             0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], xid_all[:-1, :, i_new])
            # sys.exit()
            xi_95min = percentiles[0, :]
            xi_95max = percentiles[1, :]
            xi_68min = percentiles[2, :]
            xi_68max = percentiles[3, :]
            xi_38min = percentiles[4, :]
            xi_38max = percentiles[5, :]
            xi_median = percentiles[6, :]
            xi_10min = percentiles[7, :]
            xi_10max = percentiles[8, :]

            ax1.fill_betweenx(depths[:-1], xi_95min, xi_95max, color='r', alpha=0.2)
            ax1.fill_betweenx(depths[:-1], xi_68min, xi_68max, color='r', alpha=0.2)
            ax1.fill_betweenx(depths[:-1], xi_38min, xi_38max, color='r', alpha=0.2)
            ax1.fill_betweenx(depths[:-1], xi_10min, xi_10max, color='r', alpha=0.2)
            ax1.plot(xi_median, depths[:-1], c='r', linewidth=1)

            percentiles = weighted_quantile(vps,
                                            [0.5 - 0.95 / 2, 0.50 + 0.95 / 2, 0.50 - 0.68 / 2, 0.50 + 0.68 / 2,
                                             0.50 - 0.38 / 2,
                                             0.50 + 0.38 / 2, 0.50, 0.5-0.05, 0.5 + 0.05], vpd_all[:-1, :, i_new])
            # sys.exit()
            vp_95min = percentiles[0, :]
            vp_95max = percentiles[1, :]
            vp_68min = percentiles[2, :]
            vp_68max = percentiles[3, :]
            vp_38min = percentiles[4, :]
            vp_38max = percentiles[5, :]
            vp_median = percentiles[6, :]
            vp_10min = percentiles[7, :]
            vp_10max = percentiles[8, :]

            ax2.fill_betweenx(depths[:-1], vp_95min, vp_95max, color='r', alpha=0.2)
            ax2.fill_betweenx(depths[:-1], vp_68min, vp_68max, color='r', alpha=0.2)
            ax2.fill_betweenx(depths[:-1], vp_38min, vp_38max, color='r', alpha=0.2)
            ax2.fill_betweenx(depths[:-1], vp_10min, vp_10max, color='r', alpha=0.2)
            ax2.plot(vp_median, depths[:-1], c='r', linewidth=1)

            plt.setp(ax2.get_yticklabels(), visible=False)

            ax0.plot(vsv_av[:, i_new], depths_average, c='b', linewidth=1)
            ax1.plot(xi_av[:, i_new], depths_average, c='b', linewidth=1)
            ax2.plot(vp_av[:, i_new], depths_average, c='b', linewidth=1)
            ax4.plot(probani_all[:, i_new] * 100, depths_average, c='k', linewidth=1)

            ax4.set_xlabel('anisotropy probability', fontsize=15)
            ax4.set_xlim([0, 100])
            legend_elements = [Patch(facecolor='r', edgecolor='r', alpha=0.2,
                         label='95%'),Patch(facecolor='r', edgecolor='r', alpha=0.2*0.8+0.2,
                         label='68%'),Patch(facecolor='r', edgecolor='r', alpha=(0.2*0.8+0.2)*0.8+0.2,
                         label='38%'),Patch(facecolor='r', edgecolor='r', alpha=((0.2*0.8+0.2)*0.8+0.2)*0.8+0.2,
                         label='10%'),Line2D([0], [0], color='r', lw=1, label='median'),
                         Line2D([0], [0], color='b', lw=1, label='mean')]
            ax2.legend(handles=legend_elements)

            fig.suptitle('posterior and averages ' + str(i + 1) + ', lat: ' + str(lat[i]) + ', lon: ' + str(lon[i]))

            plt.savefig(directory + '/PLOTS/Posterior/' + 'posterior_' + str(i) + '_quantiles.pdf')  #

            # for cl in clusters:
            #     if i in clusters[cl]:

            #         plt.savefig(directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/posterior_quantiles.pdf')

            plt.close(fig)

        for i in range(num_files):
            if fs[i] is not None: fs[i].close()

        for i in range(num_files):
            if fs2[i] is not None: fs2[i].close()

#################################################
# Sigma R
#################################################

if results and sigma:
    # reference

    print('plotting uncertainty parameters')

    function = 'get_histograms'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)
            for j in points:
                i_file = j // batch
                if i_file == i:
                    points.remove(j)
    for i in range(num_files):
        if fs[i] is None:
            continue
        sigmaRrange = np.array(fs[i][function]['nostack'].get('range_sigmaR'))
        sigmaLrange = np.array(fs[i][function]['nostack'].get('range_sigmaL'))

        sigmaR_all = np.array(fs[i][function]['stack'].get('sigmaR_hist'))
        sigmaL_all = np.array(fs[i][function]['stack'].get('sigmaL_hist'))

    plt.figure('sigmad_R cluster')
    plt.title('sigmad_R cluster, ' + directory)
    plt.plot(sigmaRrange, sigmaR_all[:, -1], label='centroid')
    plt.legend()
    plt.savefig(directory + '/PLOTS/' + 'sigma_R_cluster.pdf')
    plt.close()

    plt.figure('sigmad_L cluster')
    plt.title('sigmad_L cluster, ' + directory)
    plt.plot(sigmaLrange, sigmaL_all[:, -1], label='centroid')
    plt.legend()
    plt.savefig(directory + '/PLOTS/' + 'sigma_L_cluster.pdf')
    plt.close()

    # Tradeoff for cluster members
    if many_plots:

        # with PdfPages(directory+'/PLOTS/tradeoffs.pdf') as pdf:
        for i in points:
            i_file = i // batch
            i_new = i % batch
            if fs[i_file] is None:
                continue

            sigmaR_all = np.array(fs[i_file][function]['stack'].get('sigmaR_hist'))
            sigmaL_all = np.array(fs[i_file][function]['stack'].get('sigmaL_hist'))

            # print(i)
            fig = plt.figure('sigmad_R')
            plt.title('sigmad_R cluster ' + str(i + 1))
            plt.plot(sigmaRrange, sigmaR_all[:, i_new], label=str(i + 1))

            # pdf.savefig()
            plt.savefig(directory + '/PLOTS/SigmaRs/' + 'sigmaR_' + str(i) + '.pdf')

            # for cl in clusters:
            #     if i in clusters[cl]:

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(
            #             directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/sigmaR.pdf')

            plt.close(fig)

            fig = plt.figure('sigmad_L')
            plt.title('sigmad_L cluster ' + str(i + 1) + ', lat: ' + str(lat[i]) + ', lon: ' + str(lon[i]))
            plt.plot(sigmaLrange, sigmaL_all[:, i_new], label=str(i + 1))

            # pdf.savefig()
            plt.savefig(directory + '/PLOTS/SigmaLs/' + 'sigmaL_' + str(i) + '.pdf')

            # for cl in clusters:
            #     if i in clusters[cl]:

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(
            #             directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/sigmaL.pdf')

            plt.close(fig)
            # del fig
        # del l
        # del lines
    for i in range(num_files):
        if fs[i] is not None: fs[i].close()

###################################################
# Dispersion curves
###################################################

if results and dispersion:

    print('plotting dispersion curves')

    function = 'get_dispersion_mean'
    fs = []
    for i in range(num_files):
        file = directory + '/' + 'processing_' + '_' + function + '_outputs_' + str(i) + '.h5'
        if os.path.isfile(file):
            fs.append(h5py.File(file, 'r'))
        else:
            fs.append(None)

    # f.close()
    # reference
    file = directory + '/' + 'dispersion.h5'
    f = h5py.File(file, 'r')

    # numdis=int(f['cluster_params'].get('numdis'))
    period_R_obs = np.array(f['dispersion_params']['R'].get('periods'))
    period_R = period_R_obs
    n_R_obs = np.array(f['dispersion_params']['R'].get('modes'))
    n_R = n_R_obs
    # c_R_obs=np.array(f['reference']['R'].get('dispersion'))
    # dc_R_obs=np.array(f['reference']['R'].get('error'))

    period_L_obs = np.array(f['dispersion_params']['L'].get('periods'))
    period_L = period_L_obs
    n_L_obs = np.array(f['dispersion_params']['L'].get('modes'))
    n_L = n_L_obs

    c_R_obs = np.array(f['cluster']['R'].get('dispersion'))
    dc_R_obs = np.array(f['cluster']['R'].get('error'))

    c_L_obs = np.array(f['cluster']['L'].get('dispersion'))
    dc_L_obs = np.array(f['cluster']['L'].get('error'))

    f.close()

    ndatad_R = len(period_R_obs)
    ndatad_L = len(period_L_obs)

    for i in range(num_files):
        if fs[i] is None:
            continue

        c_R = np.array(fs[i][function]['stack'].get('dispersion_R'))
        dc_R = np.array(fs[i][function]['stack'].get('dispersion_R_sq'))

        c_L = np.array(fs[i][function]['stack'].get('dispersion_L'))
        dc_L = np.array(fs[i][function]['stack'].get('dispersion_L_sq'))

    if many_plots:

        # numdis=np.shape(c_L)[1]

        for i in points:
            i_file = i // batch
            i_new = i % batch
            if fs[i_file] is None:
                continue

            c_R = np.array(fs[i_file][function]['stack'].get('dispersion_R'))
            dc_R = np.array(fs[i_file][function]['stack'].get('dispersion_R_sq'))

            c_L = np.array(fs[i_file][function]['stack'].get('dispersion_L'))
            dc_L = np.array(fs[i_file][function]['stack'].get('dispersion_L_sq'))

            # print(i)
            fig = plt.figure('dispersion')
            for j in np.unique(n_R_obs):
                ints = np.where(n_R_obs == j)
                if j == 0:
                    plt.errorbar(period_R[ints[0]], c_R[ints[0], i_new], yerr=dc_R[ints[0], i_new], marker='o',
                                 zorder=0, label='Rayleigh average ' + str(i + 1), c='C0', elinewidth=0.1,
                                 linewidth=0.1, markersize=0.3)
                    plt.errorbar(period_R_obs[ints[0]] - 0.1, c_R_obs[ints[0], i], yerr=dc_R_obs[ints[0], i],
                                 marker='o', zorder=0, label='Rayleigh observed ' + str(i + 1), c='C1',
                                 elinewidth=0.1,
                                 linewidth=0.1, markersize=0.3)
                else:
                    plt.errorbar(period_R[ints[0]], c_R[ints[0], i_new], yerr=dc_R[ints[0], i_new], marker='o',
                                 zorder=0, c='C0', elinewidth=0.1, linewidth=0.1, markersize=0.3)
                    plt.errorbar(period_R_obs[ints[0]] - 0.1, c_R_obs[ints[0], i], yerr=dc_R_obs[ints[0], i],
                                 marker='o', zorder=0, c='C1', elinewidth=0.1, linewidth=0.1, markersize=0.3)
            for j in np.unique(n_L_obs):
                ints = np.where(n_L_obs == j)
                if j == 0:
                    plt.errorbar(period_L[ints[0]] + 0.1, c_L[ints[0], i_new], yerr=dc_L[ints[0], i_new],
                                 marker='o',
                                 zorder=0, label='Love average ' + str(i + 1), c='C2', elinewidth=0.1,
                                 linewidth=0.1,
                                 markersize=0.3)

                    plt.errorbar(period_L_obs[ints[0]] + 0.2, c_L_obs[ints[0], i], yerr=dc_L_obs[ints[0], i],
                                 marker='o', zorder=0, label='Love observed ' + str(i + 1), c='C3', elinewidth=0.1,
                                 linewidth=0.1, markersize=0.3)
                else:
                    plt.errorbar(period_L[ints[0]] + 0.1, c_L[ints[0], i_new], yerr=dc_L[ints[0], i_new],
                                 marker='o',
                                 zorder=0, c='C2', elinewidth=0.1, linewidth=0.1, markersize=0.3)

                    plt.errorbar(period_L_obs[ints[0]] + 0.2, c_L_obs[ints[0], i], yerr=dc_L_obs[ints[0], i],
                                 marker='o', zorder=0, c='C3', elinewidth=0.1, linewidth=0.1, markersize=0.3)
            plt.legend()
            plt.xlabel('period, s')
            plt.ylabel("phase velocity, km/s")
            plt.title('compared dispersion curves, model ' + str(i + 1))

            # pdf.savefig()
            plt.savefig(directory + '/PLOTS/Dispersion/' + 'dispersion_' + str(i) + '.pdf')
            # plt.show()

            # for cl in clusters:
            #     if i in clusters[cl]:
            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl))

            #         if not os.path.exists(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1)):
            #             os.mkdir(directory + '/PLOTS/Individuals/' + str(cl) + '/' + str(i + 1))

            #         plt.savefig(
            #             directory + '/PLOTS/' + 'Individuals/' + str(cl) + '/' + str(i + 1) + '/dispersion.pdf')

            plt.close(fig)
    for i in range(num_files):
        if fs[i] is not None: fs[i].close()
