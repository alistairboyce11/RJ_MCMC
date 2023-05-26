#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:43:33 2022

@author: dorian
"""
import numpy as np
import os
import glob
import time
import h5py
import sys
from mpi4py import MPI
import postprocess_util
import postprocess_functions
import subprocess


def process_alphamax(comm, maxpercent, alphafile, dispersion_ref, dispersion_all, outdir_ref, points, widenings):
    size = comm.Get_size()
    rank = comm.Get_rank()

    if maxpercent == 0.:
        return False

    alpha_max = []

    #if rank == 0:
    # get which points are in which file
    whoiswhere = {}
    for file in glob.glob(outdir_ref + '/' + 'processing_' + '_' + 'get_histograms' + '_outputs_*.h5'):
        filenum = int(file.split('.')[-2].split('_')[-1])
        f = h5py.File(file, 'r')
        whoiswhere[file] = (np.array(f['get_histograms']['nostack'].get('points')), filenum)
        f.close()

    # get the file for each point
    file_loc = []
    for point in points:
        for file in whoiswhere:
            if point in whoiswhere[file][0]:
                file_loc.append((file, whoiswhere[file][1]))

    # if I don't have points for all files, calculate normally
    if len(file_loc) < len(points):
        return False

    # If I have all the files I need
    # get total amount of models
    clusters = params_dispersion['clusterlist']
    num_model_all = 0
    for cluster in clusters:
        for widening in widenings:
            num_model_all += params_dispersion['num_models'][cluster][widening]

    # amount of models capped
    num_alphas = num_model_all * maxpercent / 100

    # get alphamax histogram for each point
    alpha_max = []
    for i, point in enumerate(points):

        filenum = file_loc[i][1]
        f = h5py.File(file_loc[i][0], 'r')
        ind = np.where(np.array(f['get_histograms']['nostack'].get('points')) == point)[0][0]
        exp_sum = np.array(f['get_histograms']['stack'].get('exp_sum'))[ind]

        range_alpha = np.flip(np.array(f['get_histograms']['nostack'].get('range_alpha')))
        indmax=len(range_alpha)
        alphaline = np.flip(np.array(f['get_histograms']['stack'].get('alpha_hist'))[:, ind])*exp_sum
        f.close()

        ind_alphamax = np.argmin(np.abs(np.cumsum(alphaline) - num_alphas))

        if ind_alphamax == indmax:
            print('capping too strong for it to work')
            #return False

        # shift caused by capping
        alphamax_shift = range_alpha[ind_alphamax]

        # get alphamax without capping for every point
        if not os.path.isfile(outdir_ref + '/alphamax_' + str(filenum) + '.txt'):
            return False
        f = open(outdir_ref + '/alphamax_' + str(filenum) + '.txt')
        alphamax_file = f.readline().split()
        f.close()

        # sum the two
        alpha_max.append(float(alphamax_file[ind]) + alphamax_shift)
    

    if rank == 0:
        file_alpha = open(alphafile, 'w')
        file_alpha.write(' '.join(np.array(alpha_max).astype('str')) + '\n')
        for cluster in num_models_dict:
            for widening in num_models_dict[cluster]:
                file_alpha.write(
                    ' '.join([str(cluster), str(widening), str(num_models_dict[cluster][widening])]) + '\n')
        file_alpha.write(str(len(alpha_max)))
        file_alpha.close()
    alpha_max = comm.bcast(alpha_max, root=0)
    comm.Barrier()
    dispersion_ref['alpha_max'] = 0.
    dispersion_all['alpha_max'] = alpha_max
    return True


def get_alpha_max(comm, directories, widenings_in, clusters_in, params_dispersion, dispersion_ref, dispersion_all,
                  alphafile='OUT_KL/Processing/alpha_max.txt', outdir_ref='', maxpercent=0., thin_target=None):
    '''
    gets the maximum of alpha (ratios)
    Code runs in parallel, computes alpha_max for each file independently
    if no capping, takes the maximum of the results for each file
    if capping, takes the average of the results for each file

    Parameters
    ----------
    directory: str
        directory name where all the data is.
    params_dispersion : dict
        contains the metadata of the dispersion curves.
    dispersion_ref : dict
        contains the reference dispersion curve.
    dispersion_all : dict
        contains all dispersion curves.
    maxpercent : float, optional
        percentile at which we are capping the alphas. if 0, no capping. The default is 0..
    widening : float, optional
        which part of the preprocessing are we looking at? The default is 1..

    Returns
    -------
    None.
    Modifies dispersion_all and dispersion_ref in place, adding 'alpha_max' keys to them with the corresponding values

    '''

    numdis = dispersion_all['numdis']

    num_models_dict = params_dispersion['num_models']

    if os.path.isfile(alphafile):
        file = open(alphafile, 'r')
        dispersion_ref['alpha_max'] = 0.
        dispersion_all['alpha_max'] = np.array(file.readline().split()).astype('float')
        file.close()
        return

    files_all = []

    clusters_in = params_dispersion['clusterlist']

    # initiate stuff
    for cluster_cur in clusters_in:
        for directory in params_dispersion['directories'][cluster_cur]:
            for widening_alpha in widenings_in:
                # list files
                files_all_tmp = glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening_alpha)
                files_all_tmp.sort()
                files_all.extend(files_all_tmp)

    # create array for max alpha, potentially more than one if maxpercent is bigger than zero
    num_models_all = 0
    for cl in clusters_in:
        for w in widenings_in:
            num_models_all += num_models_dict[cl][w]
    if maxpercent == 0:
        num_to_store = 1
    else:
        num_to_store = round(num_models_all * maxpercent / 100.)

        print('num to store', num_to_store)

    numdis = dispersion_all['numdis']

    alpha_max = np.ones((num_to_store, numdis)) * float('-inf')

    # split it over the cores
    size_alpha = comm.Get_size()
    rank_alpha = comm.Get_rank()
    rank2 = rank_alpha

    if rank_alpha == 0:
        print(files_all)

    whodoeswhat = {}
    for i in range(len(files_all)):
        whodoeswhat[files_all[i]] = i % size_alpha

    for cluster in clusters_in:
        for widening_alpha in widenings_in:
            files_all_tmp = []
            for directory in params_dispersion['directories'][cluster]:
                files_all_tmp.extend(glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening_alpha))
            files_all_tmp.sort()
            files = []
            for i in range(len(files_all_tmp)):
                if i % size_alpha == rank2:
                    files.append(files_all_tmp[i])

            # if more cores than files: go to the next iteration of the loop
            # but keep core 0 in, it's the core that gathers everything
            computing = []
            if len(files) == 0:
                rank2 -= len(files_all_tmp)
                if rank_alpha != 0:
                    continue
            for f in files_all_tmp:
                try:
                    computing.append(whodoeswhat[f])
                except:
                    pass

            # alpha_ref_max_tmp=np.ones((num_to_store))*float('-inf')
            alpha_max_tmp = np.ones((num_to_store, numdis)) * float('-inf')
            if rank_alpha == 0:
                print(computing)

            # print(files)
            num_models_cur = num_models_dict[cluster][widening_alpha]
            # loop over the files
            for i in range(len(files)):
                file = files[i]

                file_probsum=outdir_ref+'/'+file
                # file_probsum = None

                # process one file
                print('computing alphamax for file', file, 'rank', rank_alpha, 'cluster', cluster, 'widening',widening_alpha)
                alpha_max_prop, num_models_prop = postprocess_util.process_file_alphamax_all(file, 
                                                                                            widenings_in, cluster,
                                                                                            widening_alpha,
                                                                                            params_dispersion,
                                                                                            dispersion_ref,
                                                                                            dispersion_all, maxpercent,
                                                                                            num_to_store,
                                                                                            num_models_cur, file_probsum=file_probsum, thin_target=thin_target)

                # merge results on one core
                alpha_max_tmp = np.append(alpha_max_tmp, alpha_max_prop, axis=0)
                alpha_max_tmp.sort(axis=0)
                alpha_max_tmp = alpha_max_tmp[-num_to_store:, :]
            computing[:] = [x for x in computing if x != 0]
            # merge results between cores on core 0
            if rank_alpha in np.unique(computing):
                print('sending', rank_alpha)
                comm.Ssend([alpha_max_tmp, MPI.FLOAT], dest=0, tag=1 + int(widening_alpha) * 100 + int(cluster) * 1000)
                print('sent', rank_alpha)
            if rank_alpha == 0:
                # recieve and merge data from all other cores
                print('recieving data for cluster '+str(cluster)+', widening '+str(widening_alpha))
                for i in np.unique(computing):
                    print('recieving', i)
                    alpha_max_prop = np.zeros_like(alpha_max_tmp)
                    comm.Recv([alpha_max_prop, MPI.FLOAT], i, 1 + int(widening_alpha) * 100 + int(cluster) * 1000)
                    alpha_max_tmp = np.append(alpha_max_tmp, alpha_max_prop, axis=0)
                    alpha_max_tmp.sort(axis=0)
                    alpha_max_tmp = alpha_max_tmp[-num_to_store:, :]

                    print('recieved', i)

                # take max over all clusters/widenings (incremental)
                alpha_max = np.append(alpha_max, alpha_max_tmp, axis=0)
                alpha_max.sort(axis=0)
                alpha_max = alpha_max[-num_to_store:, :]

                print('computed alphamax for cluster', cluster, 'and widening', widening_alpha)
            # special treatement for core 0, else
            # shift to take into account we have more cores than files (potentially)
            if rank_alpha != 0 or len(files) != 0:
                rank2 += size_alpha * len(files) - len(files_all_tmp)
    if num_models_all == 0:
        return

    if rank_alpha == 0:
        # handle alphamax on core 0

        # get true alpha_max: min of alpha_max
        print('getting true alpha_max')
        alpha_max.sort(axis=0)
        alpha_max = alpha_max[-num_to_store, :]
        # print alpha_max to file
        print('writing')
        file_alpha = open(alphafile, 'w')
        file_alpha.write(' '.join(np.array(alpha_max).astype('str')) + '\n')
        for cluster in num_models_dict:
            for widening_alpha in num_models_dict[cluster]:
                file_alpha.write(
                    ' '.join([str(cluster), str(widening_alpha), str(num_models_dict[cluster][widening_alpha])]) + '\n')
        file_alpha.write(str(len(alpha_max)))
        file_alpha.close()
        print('written')

        # broadcast to all other cores
        print('broadcasting')
    alpha_max = comm.bcast(alpha_max, root=0)
    comm.Barrier()
    # put alpha_max in dispersion_all and dispersion_ref for its cluster and widening
    dispersion_ref['alpha_max'] = 0.
    dispersion_all['alpha_max'] = alpha_max

    return


def apply_stuff(comm, directories, widenings_in, clusters_in, functions, params_dispersion, dispersion_ref,
                dispersion_all, model_ref, outdir_ref='', thin_target=None):
    '''
    Takes a list of functions, reads all models, applies each function to all of the models and stacks the results

    Parameters
    ----------
    directory: str
        directory name where all the data is.
    functions : list
        list of functions.
    params_inversion : dict

    params_dispersion : dict

    dispersion_ref : dict

    dispersion_all : dict

    model_ref : dict

    prepare : bool, optional

    widening : float, optional


    Returns
    -------
    outputs_all : dict
        contains the outputs of all the functions.

    The functions must take the following inputs:
    model: dict
        contains the information on the current model. Has keys
        npt_true: number of points for mineos
        npt: number of layer added in the inversion
        npt_ani: number of anisotropic layers added
        Ad_R, Ad_L: noise parameters
        depth,vsv,xi,vp: parameters of the model. All numpy arrays, vp contains the deviation of vpv/vsv to the reference
        like_w: log_likelihood of the model

    dispersion_one: dict
        contains the information on the current dispersion curve
        similar to dispersion_ref, has 3 additionnal keys:
        Ad_R and Ad_L: noise parameters
        like_w: log-likelihood

    params_dispersion: dict

    params_inversion: dict

    dispersion_ref: dict

    dispersion_all: dict

    alpha_ref: float

    alpha: numpy array

    first: bool
    whether it is the first model to be processes on the file. May provide a minor speed boost.

    The function must give the following output:
    output: dict
        contains the output of the function. Must have 2 keys, 'stack' and 'nostack' that each have sub-dicts
        values in output['stack'] will be added for key and each model
        values in output['nostack'] will be kept as in the first model
    '''

    # initialise outputs_all
    outputs_all = {}
    params_inversion = postprocess_util.get_metadata(directories, widenings_in)
    for function_process in functions:
        # print(params_inversion)
        outputs_all[function_process.__name__] = function_process('', {}, model_ref, {'cluster': 0}, params_dispersion,
                                                                  params_inversion, dispersion_all, np.zeros((1)), outputs={},
                                                                  init=True)

    # initialise weights
    numdis = dispersion_all['numdis']
    alpha_sum = np.zeros((numdis + 1))
    num_models_process = 0
    exp_sum = np.zeros_like(alpha_sum)

    files_all = []
    clusters_all = []
    widenings_all = []

    clusters = params_dispersion['clusterlist']

    for cluster_cur in clusters:
        for directory in params_dispersion['directories'][cluster_cur]:
            for widening_process in widenings_in:
                # list files
                files_all.extend(glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening_process))
                numfiles = len(glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening_process))
                clusters_all += numfiles * [cluster_cur]
                widenings_all += numfiles * [widening_process]

    # split it over the cores
    size = comm.Get_size()
    rank = comm.Get_rank()
    files = []
    clusters = []
    widenings = []
    for i_file_process in range(len(files_all)):
        if i_file_process % size == rank:
            files.append(files_all[i_file_process])
            clusters.append(clusters_all[i_file_process])
            widenings.append(widenings_all[i_file_process])
    print(files)

    # loop over files
    for i_file_process in range(len(files)):

        file = files[i_file_process]
        cluster_process = clusters[i_file_process]
        widening_process = widenings[i_file_process]
        print(file, rank)
        if params_dispersion['num_models'][cluster_process][widening_process] == 0:
            print('no models')
            continue

        file_probsum = outdir_ref + '/' + file
        # file_probsum = None

        (alpha_sum_prop, num_models_prop, exp_sum_prop) = postprocess_util.process_one_file_all(file, cluster_process,
                                                                                                widenings_in, functions,
                                                                                                params_dispersion,
                                                                                                dispersion_ref,
                                                                                                dispersion_all,
                                                                                                model_ref, outputs_all, file_probsum=file_probsum, thin_target=thin_target)
        print('done', rank, file)

        alpha_sum += alpha_sum_prop
        num_models_process += num_models_prop
        exp_sum += exp_sum_prop
        # for key in outputs_all['get_dispersion_mean']['stack']:
        #    print(key,outputs_all['get_dispersion_mean']['stack'][key])
    # send and gather alpha, alphasum, num_models
    if rank != 0:
        comm.Ssend([alpha_sum, MPI.FLOAT], dest=0, tag=101)
        comm.Ssend([np.array([num_models_process]), MPI.INT], dest=0, tag=103)
        comm.Ssend([exp_sum, MPI.FLOAT], dest=0, tag=104)

    elif rank == 0:
        for i_file_process in range(size)[1:]:
            alpha_sum_prop = np.zeros_like(alpha_sum)
            comm.Recv([alpha_sum_prop, MPI.FLOAT], i_file_process, 101)
            alpha_sum += alpha_sum_prop

            num_models_prop = np.zeros_like(num_models_process)
            comm.Recv([num_models_prop, MPI.INT], i_file_process, 103)
            num_models_process += num_models_prop

            exp_sum_prop = np.zeros_like(exp_sum)
            comm.Recv([exp_sum_prop, MPI.FLOAT], i_file_process, 104)
            exp_sum += exp_sum_prop

    # send and gather function outputs
    for function_process in outputs_all:
        for key_process in outputs_all[function_process]['stack']:

            # make every float an array for sending
            if type(outputs_all[function_process]['stack'][key_process]) == float:
                outputs_all[function_process]['stack'][key_process] = np.array([outputs_all[function_process]['stack'][key_process]])

            # send to core 0
            if rank != 0:
                # print(function,key)
                comm.Send([outputs_all[function_process]['stack'][key_process], MPI.FLOAT], dest=0, tag=3)

            # recieve and sum up
            if rank == 0:
                for i_file_process in range(size)[1:]:
                    outputs_all_prop = np.zeros_like(outputs_all[function_process]['stack'][key_process])
                    comm.Recv([outputs_all_prop, MPI.FLOAT], i_file_process, 3)
                    outputs_all[function_process]['stack'][key_process] += outputs_all_prop
                # normalise by sum of weights

                # print('before norm',key,outputs_all['get_dispersion_mean']['stack'][key])
                outputs_all[function_process]['stack'][key_process][..., :] /= exp_sum
                # print('exp_sum',exp_sum)
                # print('after norm',key,outputs_all['get_dispersion_mean']['stack'][key])

    # calculate standart deviations
    if rank == 0:
        for function_process in outputs_all:

            if 'meansum_sq' in outputs_all[function_process].keys():
                for i_file_process in range(len(outputs_all[function_process]['meansum_sq'])):

                    # appropriate keys
                    key_mean_sq = outputs_all[function_process]['meansum_sq'][i_file_process]
                    key_mean_shift = outputs_all[function_process]['meansum_sq_shift'][i_file_process]

                    # mean of squares, shifted
                    mean_sq = outputs_all[function_process]['stack'][key_mean_sq]

                    # square of means, shifted
                    mean_shifted = np.square(outputs_all[function_process]['stack'][key_mean_shift])
                    try:
                        outputs_all[function_process]['stack'][key_mean_sq] = np.sqrt(mean_sq - mean_shifted)
                    except:
                        pass

                    del outputs_all[function_process]['stack'][key_mean_shift]
            outputs_all[function_process]['stack']['exp_sum'] = exp_sum
            outputs_all[function_process]['stack']['num_models'] = num_models_process

    return outputs_all  # careful, it is only on core 0!


############################################################################
#        Example
############################################################################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(size, rank)

# directories=['OUT_POINTS_1228_RUN/OUT']
# clusters=[0]
directories = ['OUT_CLUSTER0_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER1_NEW_CLUSTERS_RUN/OUT',
              'OUT_CLUSTER2_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER3_NEW_CLUSTERS_RUN/OUT',
              'OUT_CLUSTER4_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER5_NEW_CLUSTERS_RUN/OUT',
              'OUT_CLUSTER6_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER7_NEW_CLUSTERS_RUN/OUT',
              'OUT_CLUSTER8_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER0_NEW_CLUSTERS_RUN2/OUT',
              'OUT_CLUSTER1_NEW_CLUSTERS_RUN2/OUT', 'OUT_CLUSTER2_NEW_CLUSTERS_RUN2/OUT',
              'OUT_CLUSTER3_NEW_CLUSTERS_RUN2/OUT', 'OUT_CLUSTER4_NEW_CLUSTERS_RUN2/OUT',
              'OUT_CLUSTER5_NEW_CLUSTERS_RUN2/OUT', 'OUT_CLUSTER6_NEW_CLUSTERS_RUN2/OUT',
              'OUT_CLUSTER7_NEW_CLUSTERS_RUN2/OUT', 'OUT_CLUSTER8_NEW_CLUSTERS_RUN2/OUT']
#              'OUT_CLUSTER0_DEEP/OUT', 'OUT_CLUSTER1_DEEP/OUT',
#              'OUT_CLUSTER2_DEEP/OUT', 'OUT_CLUSTER3_DEEP/OUT', 
#              'OUT_CLUSTER4_DEEP/OUT']
clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8]#, 10, 11, 12, 13, 14]

directories = ['OUT_CLUSTER0_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER1_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER2_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER3_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER4_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER5_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER6_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER7_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER8_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER0_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER1_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER2_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER3_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER4_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER5_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER6_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER7_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER8_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER0_DEEP/OUT', 'OUT_CLUSTER1_DEEP/OUT',
              'OUT_CLUSTER2_DEEP/OUT', 'OUT_CLUSTER3_DEEP/OUT', 
              'OUT_CLUSTER4_DEEP/OUT']
clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14]

directories = ['OUT_CLUSTER0_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER1_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER2_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER3_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER4_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER5_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER6_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER7_NEW_CLUSTERS_RUN/OUT_250',
              'OUT_CLUSTER8_NEW_CLUSTERS_RUN/OUT_250', 'OUT_CLUSTER0_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER1_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER2_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER3_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER4_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER5_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER6_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER7_NEW_CLUSTERS_RUN2/OUT_250', 'OUT_CLUSTER8_NEW_CLUSTERS_RUN2/OUT_250',
              'OUT_CLUSTER0_DEEP/OUT_250', 'OUT_CLUSTER1_DEEP/OUT_250',
              'OUT_CLUSTER2_DEEP/OUT_250', 'OUT_CLUSTER3_DEEP/OUT_250', 
              'OUT_CLUSTER4_DEEP/OUT_250',
              'OUT_CLUSTER24_OCEANS/OUT','OUT_CLUSTER25_OCEANS/OUT','OUT_CLUSTER26_OCEANS/OUT']
clusters = [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 24, 25, 26]

#directories = ['OUT_CLUSTER0_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER1_NEW_CLUSTERS_RUN/OUT',
#              'OUT_CLUSTER0_NEW_CLUSTERS_RUN2/OUT',
#              'OUT_CLUSTER1_NEW_CLUSTERS_RUN2/OUT', ]
#clusters = [0, 1, 0, 1]

# directories = ['OUT_CLUSTER8_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER8_NEW_CLUSTERS_RUN2/OUT']
# clusters = [8,8]

# directories = ['OUT_CLUSTER2_NEW_CLUSTERS_RUN/OUT', 'OUT_CLUSTER2_NEW_CLUSTERS_RUN2/OUT']
# clusters = [2,2]

# directories=['OUT_POINTS_1228_RUN/OUT']
# clusters=[1228]

#directories=['OUT_POINTS_15254_RUN/OUT']
#clusters=[15254]

# directories=['OUT_CLUSTER0_NARROW/OUT','OUT_CLUSTER1_NARROW/OUT','OUT_CLUSTER2_NARROW/OUT','OUT_CLUSTER6_NARROW/OUT','OUT_CLUSTER7_NARROW/OUT','OUT_CLUSTER8_NARROW/OUT','OUT_CLUSTER9_NARROW/OUT','OUT_CLUSTER0_NARROW_HIGHT/OUT','OUT_CLUSTER1_NARROW_HIGHT/OUT','OUT_CLUSTER2_NARROW_HIGHT/OUT','OUT_CLUSTER6_NARROW_HIGHT/OUT','OUT_CLUSTER7_NARROW_HIGHT/OUT','OUT_CLUSTER8_NARROW_HIGHT/OUT','OUT_CLUSTER9_NARROW_HIGHT/OUT']
# directories.extend(['OUT_CLUSTER'+str(i)+'_NARROW/OUT' for i in range(10,30)])
# clusters=[0,1,2,6,7,8,9,0,1,2,6,7,8,9]
# clusters.extend(range(10,30))

maxpercent = 10.
suffix = '_10.'
batch = 500
outdir_ref = "OUT_ALL_DEEP_OCEAN"# 'OUT_POINTS_15254'#'
#outdir_ref = 'OUT_POINTS_15254'
outdir = outdir_ref
widenings = [1., 2., 3., 4.]
points = 'all'  # [1228,1273,1487,5459,1839,9215,3989,11172,11123,11186]
thin_target=None

outdir = outdir_ref + suffix

if not os.path.exists(outdir_ref) and rank == 0:
    os.mkdir(outdir_ref)

if not os.path.exists(outdir) and rank == 0:
    os.mkdir(outdir)

#functions = [postprocess_functions.get_histograms]
functions = [postprocess_functions.create_posterior,postprocess_functions.get_average,postprocess_functions.get_histograms,postprocess_functions.get_dispersion_mean]#,postprocess_functions.get_tradeoff]


if rank == 0:
    print('start')
params_dispersion, dispersion_ref, dispersion_all_ref = postprocess_util.get_dispersion(directories, clusters,
                                                                                        points=points)
alphafiles_present = glob.glob(outdir + '/alphamax_*.txt')
print('counting models')
num_models_done = False
if rank == 0:
    if os.path.isfile(outdir_ref + '/dispersion.h5'):
        f = h5py.File(outdir_ref + '/dispersion.h5', 'r')
        num_models = {}
        if not 'num_models' in f:
            f.close()
        else:
            for key in f['num_models']:
                num_models[int(key)] = {}
                for key2 in f['num_models'][key]:
                    print(key, key2, f['num_models'][key][key2][...])
                    num_models[int(key)][float(key2)] = f['num_models'][key][key2][...]
            params_dispersion['num_models'] = num_models
            num_models_done = True
            f.close()
num_models_done = comm.bcast(num_models_done, root=0)
if num_models_done:
    params_dispersion = comm.bcast(params_dispersion, root=0)
else:
    postprocess_util.count_models(directories, widenings, params_dispersion, comm, thin_target=thin_target)
# else:
#    f=open(alphafiles_present[0])
#    numdis_dict={}
#    for line in f.readlines()[1:]:
#        data=line.split()
#        cl=int(data[0])
#        w=float(data[1])
#        num=int(data[2])
#        if not cl in numdis_dict:
#            numdis_dict[cl]={}
#        numdis_dict[cl][w]=num
#    f.close()
#    params_dispersion['num_models']=numdis_dict


if rank == 0:
    print('get dispersion')
model_ref = postprocess_util.get_model_ref(filename='Modified_PREM_GLOBAL.in')
if rank == 0:
    print('get model ref')

# save dispersion to file
if rank == 0:
    print('saving dispersion')
    filename = 'dispersion.h5'
    f = h5py.File(outdir_ref + '/' + filename, 'w')
    grp = f.create_group('cluster_params')
    grp.create_dataset('batchsize', data=batch)
    grp.create_dataset('numdis', data=dispersion_all_ref['numdis'])
    grp.create_dataset('lat', data=dispersion_all_ref['lat'])
    grp.create_dataset('lon', data=dispersion_all_ref['lon'])
    grp.create_dataset('clusterlist', data=np.array(params_dispersion['clusterlist']))
    grp2 = grp.create_group('cluster')
    for cl in dispersion_all_ref['cluster']:
        grp2.create_dataset(str(cl), data=dispersion_all_ref['cluster'][cl])
    grp = f.create_group('dispersion_params')
    grp2 = grp.create_group('L')
    grp2.create_dataset('periods', data=params_dispersion['L']['periods'])
    grp2.create_dataset('modes', data=params_dispersion['L']['modes'])
    grp2 = grp.create_group('R')
    grp2.create_dataset('periods', data=params_dispersion['R']['periods'])
    grp2.create_dataset('modes', data=params_dispersion['R']['modes'])
    grp = f.create_group('reference')
    for cl in dispersion_all_ref['cluster']:
        grp3 = grp.create_group(str(cl))
        grp2 = grp3.create_group('L')
        grp2.create_dataset('dispersion', data=dispersion_ref[cl]['L']['dispersion'])
        grp2.create_dataset('error', data=dispersion_ref[cl]['L']['error'])
        grp2 = grp3.create_group('R')
        grp2.create_dataset('dispersion', data=dispersion_ref[cl]['R']['dispersion'])
        grp2.create_dataset('error', data=dispersion_ref[cl]['R']['error'])
    grp = f.create_group('cluster')
    grp2 = grp.create_group('L')
    grp2.create_dataset('dispersion', data=dispersion_all_ref['L']['dispersion'])
    grp2.create_dataset('error', data=dispersion_all_ref['L']['error'])
    grp2 = grp.create_group('R')
    grp2.create_dataset('dispersion', data=dispersion_all_ref['R']['dispersion'])
    grp2.create_dataset('error', data=dispersion_all_ref['R']['error'])
    grp2 = f.create_group('num_models')
    for cluster in np.unique(clusters):
        grp3 = grp2.create_group(str(cluster))
        for widening in widenings:
            grp3.create_dataset(str(widening), data=params_dispersion['num_models'][cluster][widening])
    f.close()

    print('saving dispersion')
    filename = 'dispersion.h5'
    f = h5py.File(outdir + '/' + filename, 'w')
    grp = f.create_group('cluster_params')
    grp.create_dataset('batchsize', data=batch)
    grp.create_dataset('numdis', data=dispersion_all_ref['numdis'])
    grp.create_dataset('lat', data=dispersion_all_ref['lat'])
    grp.create_dataset('lon', data=dispersion_all_ref['lon'])
    grp.create_dataset('clusterlist', data=np.array(params_dispersion['clusterlist']))
    grp2 = grp.create_group('cluster')
    for cl in dispersion_all_ref['cluster']:
        grp2.create_dataset(str(cl), data=dispersion_all_ref['cluster'][cl])
    grp = f.create_group('dispersion_params')
    grp2 = grp.create_group('L')
    grp2.create_dataset('periods', data=params_dispersion['L']['periods'])
    grp2.create_dataset('modes', data=params_dispersion['L']['modes'])
    grp2 = grp.create_group('R')
    grp2.create_dataset('periods', data=params_dispersion['R']['periods'])
    grp2.create_dataset('modes', data=params_dispersion['R']['modes'])
    grp = f.create_group('reference')
    for cl in dispersion_all_ref['cluster']:
        grp3 = grp.create_group(str(cl))
        grp2 = grp3.create_group('L')
        grp2.create_dataset('dispersion', data=dispersion_ref[cl]['L']['dispersion'])
        grp2.create_dataset('error', data=dispersion_ref[cl]['L']['error'])
        grp2 = grp3.create_group('R')
        grp2.create_dataset('dispersion', data=dispersion_ref[cl]['R']['dispersion'])
        grp2.create_dataset('error', data=dispersion_ref[cl]['R']['error'])
    grp = f.create_group('cluster')
    grp2 = grp.create_group('L')
    grp2.create_dataset('dispersion', data=dispersion_all_ref['L']['dispersion'])
    grp2.create_dataset('error', data=dispersion_all_ref['L']['error'])
    grp2 = grp.create_group('R')
    grp2.create_dataset('dispersion', data=dispersion_all_ref['R']['dispersion'])
    grp2.create_dataset('error', data=dispersion_all_ref['R']['error'])
    grp2 = f.create_group('num_models')
    for cluster in np.unique(clusters):
        grp3 = grp2.create_group(str(cluster))
        for widening in widenings:
            grp3.create_dataset(str(widening), data=params_dispersion['num_models'][cluster][widening])
    f.close()
    print('dispersion saved')

comm.Barrier()
num_models_dict = params_dispersion['num_models']

i = 0
# for i in range(100,200):
# for i in [100]:
#for i in range(0,10):
while i * batch < dispersion_all_ref['numdis'] + 1:
    # while i==0:
    if points == "all":
        points2 = []
        points2 = range(i * batch, min([(i + 1) * batch, dispersion_all_ref['numdis']]))
        file_suffix = str(i)
    else:
        points2 = points
        file_suffix = 'test'

    params_dispersion, dispersion_ref, dispersion_all = postprocess_util.get_dispersion(directories, clusters,
                                                                                        points=points2)
    params_dispersion['num_models'] = num_models_dict
    output = {}

    if rank == 0:
        print('get alphamax')
    alphafile = outdir + '/alphamax_' + file_suffix + '.txt'

    found_alphamax = process_alphamax(comm, maxpercent, alphafile, dispersion_ref, dispersion_all, outdir_ref, points2, widenings)

    found_alphamax = comm.bcast(found_alphamax, root=0)

    if not found_alphamax:
        print('calculating alphamax')
        get_alpha_max(comm, directories, widenings, clusters, params_dispersion, dispersion_ref, dispersion_all,
                      alphafile, outdir_ref=outdir_ref, maxpercent=maxpercent, thin_target=thin_target)

    if rank == 0:
        print('apply functions')
    output = apply_stuff(comm, directories, widenings, clusters, functions, params_dispersion, dispersion_ref,
                         dispersion_all, model_ref, outdir_ref=outdir_ref, thin_target=thin_target)  #

    if rank == 0:
        print('printing result to files')

        ############################################################
        # print to files
        ############################################################
        for function in output:
            filename = 'processing_' + '_' + function + '_outputs_' + file_suffix + '.h5'
            f = h5py.File(outdir + '/' + filename, 'w')
            # f.create_dataset('burn-in',data=params_inversion['burn-in'])
            # f.create_dataset('nsample',data=params_inversion['nsample'])

            grp = f.create_group(function)
            grp_stack = grp.create_group('stack')
            grp_nostack = grp.create_group('nostack')
            grp_nostack.create_dataset('points', data=points2)
            for key in output[function]['stack']:
                grp_stack.create_dataset(key, data=output[function]['stack'][key])
                output[function]['stack'][key] = []
            for key in output[function]['nostack']:
                grp_nostack.create_dataset(key, data=output[function]['nostack'][key])
                output[function]['nostack'][key] = []
            f.close()
            output[function] = {}
    i += 1
