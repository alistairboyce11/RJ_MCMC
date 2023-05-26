#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:18:47 2022

@author: dorian
"""
import numpy as np
import glob
import time
import warnings
import concurrent.futures
from mpi4py import MPI
import os

warnings.filterwarnings("error")


class Probsum1(Exception):
    pass


def get_dispersion(directories, clusters, points=[], filename='dispersion_all.in'):
    '''

    Parameters
    ----------
    directory : str
        directory containing all the results of the inversion.
    filename : str, optional
        file used for the dispersion curve for MCMc inversions. The default is 'dispersion_all.in'.

    Returns
    -------
    params_dispersion : dict
        dictionary containing the metadata of the dispersion curves.
    dispersion_ref : dict
        dictionary containing the dispersion curve of the cluster centre.
    dispersion_all : dict
        dictionary containing the individual dispersion curves of the cluster members
    and additional information such as cluster they are members of, latitude and longitude of individual points.

    All dictionaries are, when relevant, divided into Love (L) and Rayleigh (R) data

    '''

    # initialising dicts
    dispersion_ref = {}
    dispersion_ref['R'] = {}
    dispersion_ref['L'] = {}
    dispersion_all = {}
    dispersion_all['R'] = {}
    dispersion_all['L'] = {}
    params_dispersion = {}
    params_dispersion['R'] = {}
    params_dispersion['L'] = {}
    params_dispersion['clusterlist'] = []
    params_dispersion['directories'] = {}
    dispersion_all['cluster'] = {}


    for ind_cluster in range(len(clusters)):

        if clusters[ind_cluster] in dispersion_all['cluster']:
            params_dispersion['directories'][clusters[ind_cluster]].append(directories[ind_cluster])
            continue

        params_dispersion['clusterlist'].append(clusters[ind_cluster])
        params_dispersion['directories'][clusters[ind_cluster]] = [directories[ind_cluster]]

        dispersion_ref[clusters[ind_cluster]]={}
        dispersion_ref[clusters[ind_cluster]]['R']={}
        dispersion_ref[clusters[ind_cluster]]['L']={}

        # different references for different clusters
        file = open(directories[ind_cluster] + '/' + filename, 'r')

        numdis_tmp = int(file.readline())
        if type(points) == str and points == 'all':
            numdis = numdis_tmp
        elif len(points) == 0:
            numdis = 0
        else:
            numdis = len(points)
        dispersion_all['numdis'] = numdis

        dispersion_all['lat'] = np.array(file.readline().split()).astype('float')
        dispersion_all['lon'] = np.array(file.readline().split()).astype('float')
        clusterlist = np.array(file.readline().split()).astype('float')
        clusterlist = clusterlist.astype(np.int64)
        dispersion_all['cluster'][clusters[ind_cluster]] = np.where(clusterlist == clusters[ind_cluster])[0]

        # Read rayleigh data
        ndatad_R = int(file.readline())
        params_dispersion['R']['ndatad'] = ndatad_R
        mode_rayl = int(file.readline())

        mode_R = np.zeros((ndatad_R))
        period_R = np.zeros((ndatad_R))
        dispersion_R = np.zeros((ndatad_R, numdis))
        error_R = np.zeros((ndatad_R, numdis))
        dispersion_R_ref = np.zeros((ndatad_R))
        error_R_ref = np.zeros((ndatad_R))
        j = 0
        for i in range(mode_rayl):
            file.readline()
            ind_k = int(file.readline()) + 1
            for k in range(ind_k):
                data = file.readline().split()
                mode_R[j] = int(data[0])
                period_R[j] = float(data[1])
                dispersion_R_ref[j] = float(data[2])
                error_R_ref[j] = float(data[3])

                for l in range(numdis_tmp):
                    line = file.readline()
                    # print(line)
                    data = line.split()
                    # print(data)
                    if type(points) == str and points == 'all':
                        dispersion_R[j, l] = float(data[0])
                        error_R[j, l] = float(data[1])
                    else:
                        if l not in points:
                            pass
                        else:
                            l_ind = list(points).index(l)
                            dispersion_R[j, l_ind] = float(data[0])
                            error_R[j, l_ind] = float(data[1])
                j += 1

        # fill reference dispersion curves: one line per cluster

        dispersion_ref[clusters[ind_cluster]]['R']['dispersion'] = dispersion_R_ref
        dispersion_ref[clusters[ind_cluster]]['R']['error'] = error_R_ref
        dispersion_ref[clusters[ind_cluster]]['R']['error_sum'] = np.sum(np.log(error_R_ref))

        # individual data, identical every time
        dispersion_all['R']['dispersion'] = dispersion_R
        dispersion_all['R']['error'] = error_R
        dispersion_all['R']['error_sum'] = np.sum(np.log(error_R), axis=0)
        params_dispersion['R']['periods'] = period_R
        params_dispersion['R']['modes'] = mode_R

        # read love data
        ndatad_L = int(file.readline())
        params_dispersion['L']['ndatad'] = ndatad_L
        mode_love = int(file.readline())

        mode_L = np.zeros((ndatad_L))
        period_L = np.zeros((ndatad_L))
        dispersion_L = np.zeros((ndatad_L, numdis))
        error_L = np.zeros((ndatad_L, numdis))
        dispersion_L_ref = np.zeros((ndatad_L))
        error_L_ref = np.zeros((ndatad_L))
        j = 0
        for i in range(mode_love):
            file.readline()
            ind_k = int(file.readline()) + 1
            for k in range(ind_k):
                data = file.readline().split()
                mode_L[j] = int(data[0])
                period_L[j] = float(data[1])
                dispersion_L_ref[j] = float(data[2])
                error_L_ref[j] = float(data[3])

                for l in range(numdis_tmp):
                    data = file.readline().split()
                    if type(points) == str and points == 'all':
                        dispersion_L[j, l] = float(data[0])
                        error_L[j, l] = float(data[1])
                    else:
                        if l not in points:
                            continue
                        else:
                            l_ind = list(points).index(l)
                            dispersion_L[j, l_ind] = float(data[0])
                            error_L[j, l_ind] = float(data[1])

                j += 1

        dispersion_ref[clusters[ind_cluster]]['L']['dispersion'] = dispersion_L_ref
        dispersion_ref[clusters[ind_cluster]]['L']['error'] = error_L_ref
        dispersion_ref[clusters[ind_cluster]]['L']['error_sum'] = np.sum(np.log(error_L_ref))
        
        dispersion_all['L']['dispersion'] = dispersion_L
        dispersion_all['L']['error'] = error_L
        dispersion_all['L']['error_sum'] = np.sum(np.log(error_L), axis=0)
        params_dispersion['L']['periods'] = period_L
        params_dispersion['L']['modes'] = mode_L

        file.close()

    params_dispersion['clusterlist'] = np.array(params_dispersion['clusterlist'])
    return params_dispersion, dispersion_ref, dispersion_all


def get_metadata(directories, widenings):
    '''

    gets the metadata/parameters of the inversion
    For this, it reads the first of the output files and extracts the metadata from it
    metadata is:
        burn-in: length of the burn-in
        nsample: number of samples
        widening: widening, or temperature, of the inversion.
            For a regular inversion, this can be set to one or even better skipped altogether
        d_min,d_max: depth boundaries of the inversion
        width_vsv: boundaries of the vsv prior (+- deviation from PREM)
        xi_min,xi_max: boundaries of the prior for radial anisotropy
        vpvs: reference vpv/vsv ratio
        vpvs_min,vpvs_max: maximum deviations from the reference vpv/vsv ratio
        Ad_R_min,Ad_R_max: boundaries of the prior for Rayleigh noise parameter
        Ad_L_min,Ad_L_max: boundaries of the prior for Love noise parameter

    Parameters
    ----------
    directory : str
        see get_dispersion.
    prepare : bool, optional
        see get_alpha_max. The default is False.
    widening : float, optional
        see get_alpha_max. The default is 1.0.

    Returns
    -------
    params_inversion : dict
        dict containing all the metadata of the inversion.

    '''

    params_inversion = {}

    files_all = []
    for i in range(len(directories)):
        directory = directories[i]
        for j in range(len(widenings)):
            files_all.extend(glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widenings[j]))

    for file in files_all:

        f = open(file, 'r')
        lines = []
        for i in range(9):
            lines.append(f.readline())
        # print(file)

        params_inversion = read_header(lines)
        f.close()

        if not params_inversion:
            continue
        else:
            # print('found data! ')
            return params_inversion
    return


def process_file_alphamax_all(file, widenings, cluster, widening, params_dispersion, dispersion_ref,
                              dispersion_all, maxpercent, num_to_store, num_models_cur=np.exp(1.), file_probsum=None,
                              thin_target=None):
    '''processes one file to get alphamax
    input: input_dict
    input_dict contains:
        file: filename of the file to process


    Parameters
    ----------
    input_dict : dict
        Contains file (file name), params_inversion, params_dispersion, dispersion_ref, dispersion_all and maxpercent

    Returns
    -------
    alpha_ref_max : float
        max of alpha_ref, eventually capped
    alpha_max : numpy array, length numdis
        max of alpha, eventually capped

    '''

    # reading the file
    f = open(file, 'r')
    lines = f.readlines()

    if file_probsum is None:
        pass
    elif os.path.isfile(file_probsum):
        compute_probsum = False
        f_probsum = open(file_probsum, 'r')
        lines_probsum = f_probsum.readlines()
        f_probsum.close()
        i_line_probsum = 0
        line = lines_probsum[i_line_probsum]
        num_models_prop = {}
        while len(line.split()) == 3:
            data = line.split()
            cluster_tmp = int(data[0])
            widening_tmp = float(data[1])
            num_models_tmp = int(data[2])
            if cluster_tmp in num_models_prop:
                num_models_prop[cluster_tmp][widening_tmp] = num_models_tmp
            else:
                num_models_prop[cluster_tmp] = {}
                num_models_prop[cluster_tmp][widening_tmp] = num_models_tmp
            i_line_probsum += 1
            line = lines_probsum[i_line_probsum]

        compute_probsum = params_dispersion['num_models'] != num_models_prop
    else:
        compute_probsum = True
        dir = file_probsum.rpartition('/')[0]
        os.makedirs(dir, exist_ok=True)

    num_model = 0

    numdis = dispersion_all['numdis']
    alpha_max = np.ones((num_to_store, numdis)) * float('-inf')

    params_inversion = read_header(lines[:9])

    # thinning to get all files to same 'effective thinning'
    if thin_target is None:
        thin_true = 1
    else:
        thin_true = int(round(thin_target / params_inversion['thin']))
    print('thin_true',thin_true)

    if not params_inversion:
        f.close()
        return

    widening = params_inversion['widening']

    # loop on all models
    valid_model = True
    cnt = 0
    i_model = 0
    # lines=f.readlines()
    probsums_out = []
    maxlines = len(lines)
    while valid_model:
        if 9 + 11 * (i_model) > maxlines:
            valid_model = False
            continue
        lines_model = lines[9 + 11 * i_model:9 + 11 * (i_model + 1)]
        i_model += 1
        # thinning to get all files to same 'effective thinning'
        if (i_model // params_inversion['num_proc']) % thin_true != 0:
            continue

        valid_model, dispersion_one, model = read_model(lines_model, widening, cluster)
        if not valid_model:
            continue
        elif len(dispersion_one.keys()) == 0:
            continue
        elif len(dispersion_one['R']['dispersion']) == 0:
            continue

        # elif model['npt_true']<8:
        #    continue
        cnt += 1

        if file_probsum is None:
            probsum_prop = None
        elif compute_probsum:
            probsum_prop = None
        else:
            probsum_prop = float(lines_probsum[i_line_probsum])
            i_line_probsum += 1

        # if num_model==999:
        #    valid_model=False

        # calculate alpha
        alpha, probsum = get_alpha_all(dispersion_one, widenings, params_dispersion, dispersion_ref,
                                       dispersion_all,
                                       num_models_cur, probsum_prop)
        probsums_out.append(probsum)
        if np.any(np.amin(alpha_max, axis=0) < alpha):
            min_indices = np.argmin(alpha_max, axis=0)
            alpha_max[min_indices, np.arange(numdis)] = np.maximum(alpha, np.amin(alpha_max, axis=0))
        num_model += 1

        # if (num_model%10000)==0:
        #     print(num_model)

    f.close()

    if file_probsum is None:
        pass
    elif compute_probsum:
        dir = file_probsum.rpartition('/')[0]
        os.makedirs(dir, exist_ok=True)
        f_probsum = open(file_probsum, 'w')
        for cluster_probsum in params_dispersion['num_models']:
            for widening_probsum in params_dispersion['num_models'][cluster_probsum]:
                f_probsum.write(str(cluster_probsum) + ' ' + str(widening_probsum) + ' ' + str(
                    params_dispersion['num_models'][cluster_probsum][widening_probsum]) + '\n')
        for probsum in probsums_out:
            f_probsum.write(str(probsum) + '\n')
        f_probsum.close()

    return alpha_max, num_model


def get_alpha_all(dispersion_one, widenings, params_dispersion, dispersion_ref, dispersion_all,
                  num_models_cur=np.exp(1), probsum=None):
    '''
    Gets alpha for one dispersion curve

    Parameters
    ----------
    dispersion_one : dict
        dispersion curve for one model
    params_dispersion : dict

    dispersion_ref : dict

    dispersion_all : dict


    Returns
    -------
    alpha_ref : float
        alpha for the reference dispersion curve.
    alpha : numpy array, length numdis
        alpha for all cluster members.

    '''
    clusters = params_dispersion['clusterlist']

    ndatad_R = params_dispersion['R']['ndatad']
    ndatad_L = params_dispersion['L']['ndatad']
    dispersion_R = dispersion_all['R']['dispersion']
    dispersion_L = dispersion_all['L']['dispersion']
    dispersion_R_one = dispersion_one['R']['dispersion']
    dispersion_L_one = dispersion_one['L']['dispersion']
    Ad_R = dispersion_one['R']['Ad']
    Ad_L = dispersion_one['L']['Ad']
    error_R = dispersion_all['R']['error']
    error_L = dispersion_all['L']['error']
    # error_R_sum=dispersion_all['R']['error_sum']
    # error_L_sum=dispersion_all['L']['error_sum']
    numdis = dispersion_all['numdis']

    if np.any(2 * Ad_R ** 2 * np.square(np.multiply(error_R, dispersion_R / 100.)) == 0):
        print('error division by zero', 'Ad_R', Ad_R, 'min error_R', np.amin(error_R), 'min dispersion_R',
              np.amin(dispersion_R))

    if np.any(2 * Ad_L ** 2 * np.square(np.multiply(error_L, dispersion_L / 100.)) == 0):
        print('error division by zero', 'Ad_L', Ad_L, 'min error_L', np.amin(error_L), 'min dispersion_L',
              np.amin(dispersion_L))
    like_alt = -np.sum(np.divide(
        np.square(dispersion_R - np.transpose(np.tile(dispersion_R_one, (numdis, 1)))),
        2 * Ad_R ** 2 *
        np.square(
            np.multiply(
                error_R, dispersion_R / 100.)))
        , axis=0)
    like_alt -= np.sum(np.divide(
        np.square(dispersion_L - np.transpose(np.tile(dispersion_L_one, (numdis, 1)))),
        2 * Ad_L ** 2 *
        np.square(
            np.multiply(
                error_L, dispersion_L / 100.)))
        , axis=0)
    like_alt -= ndatad_R * np.log(Ad_R)
    like_alt -= ndatad_L * np.log(Ad_L)

    if probsum is None:
        probsum = 0.
        shift = 0
        not_processed = True
        while not_processed:
            try:
                numprob = 0.
                probsum = 0.
                for cluster in clusters:
                    for widening in widenings:
                        if params_dispersion['num_models'][cluster][widening] == 0:
                            continue
                        dispersion_R_ref = dispersion_ref[cluster]['R']['dispersion']
                        dispersion_L_ref = dispersion_ref[cluster]['L']['dispersion']

                        error_R_ref = dispersion_ref[cluster]['R']['error']
                        error_L_ref = dispersion_ref[cluster]['L']['error']

                        error_R_ref_sum=dispersion_ref[cluster]['R']['error_sum']
                        error_L_ref_sum=dispersion_ref[cluster]['L']['error_sum']

                        if np.shape(dispersion_R_ref) == (0,):
                            print('here')
                        if np.shape(dispersion_R_one) == (0,):
                            print('here2')
                            print(dispersion_R_one)

                        like_one = -np.sum(np.divide(np.square(dispersion_R_ref - dispersion_R_one),
                                                     2 * Ad_R ** 2 * np.square(
                                                         np.multiply(error_R_ref / 100., dispersion_R_ref))), axis=0)
                        like_one -= np.sum(np.divide(np.square(dispersion_L_ref - dispersion_L_one),
                                                     2 * Ad_L ** 2 * np.square(
                                                         np.multiply(error_L_ref / 100., dispersion_L_ref))), axis=0)
                        like_one -= error_R_ref_sum + 150  # static shift to avoid exponential explosion
                        like_one -= error_L_ref_sum + 150  # static shift to avoid exponential explosion
                        like_one -= ndatad_R * np.log(Ad_R)
                        like_one -= ndatad_L * np.log(Ad_L)
                        like_one /= widening

                        like_one += shift

                        probsum += np.exp(like_one)

                        numprob += 1
                        not_processed = False
                if probsum == 0.:
                    raise Probsum1
            except RuntimeWarning:
                not_processed = True
                probsum = 0.
                shift -= 100
                print('too big probsum')
            except Probsum1:
                not_processed = True
                probsum = 0.
                shift += 1000
                print('too small probsum')
        probsum /= numprob

        probsum = np.log(probsum) - shift

    alpha = like_alt - probsum
    alpha -= np.log(num_models_cur)

    return alpha, probsum


def read_header(lines_header):
    params_inversion = {}

    d = lines_header.pop(
        0)  # f.readline() # contains rank of processor, number of file for the processor and number of models stored at most. Not used or tested currently
    params_inversion['num_proc'] = 112
    params_inversion['everyall'] = int(float(d.split()[1]))
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['burn-in'] = float(data[0])
    params_inversion['widening'] = float(data[1])
    params_inversion['thin'] = float(data[2])
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['d_min'] = float(data[0])
    params_inversion['d_max'] = float(data[1])
    params_inversion['width_vsv'] = float(lines_header.pop(0))  # f.readline())
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['xi_min'] = float(data[0])
    params_inversion['xi_max'] = float(data[1])
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['vp_min'] = float(data[0])
    params_inversion['vp_max'] = float(data[1])
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['Ad_R_min'] = float(data[0])
    params_inversion['Ad_R_max'] = float(data[1])
    data = lines_header.pop(0).split()  # f.readline().split()
    params_inversion['Ad_L_min'] = float(data[0])
    params_inversion['Ad_L_max'] = float(data[1])
    data = lines_header.pop(0).split()  # f.readline().split()
    if len(data) == 0:
        return
    params_inversion['milay'] = int(data[0])
    params_inversion['malay'] = int(data[1])

    return params_inversion


def read_model(lines_model, widening, cluster):
    dispersion_one = {}
    dispersion_one['R'] = {}
    dispersion_one['L'] = {}

    model = {}
    if len(lines_model) != 11:
        return False, {}, {}
    line = lines_model.pop(0)  # f.readline()
    # if not line:
    #    return False,{},{}

    data = line.split()
    npt_true = int(data[0])
    if npt_true == 0:
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        line = lines_model.pop(0)  # f.readline()
        return True, {}, {}
    # try:
    model['npt_true'] = npt_true
    model['npt'] = int(data[1])
    model['npt_ani'] = int(data[2])
    data = lines_model.pop(0).split()  # f.readline().split()
    dispersion_one['R']['Ad'] = float(data[0])
    dispersion_one['L']['Ad'] = float(data[1])
    model['Ad_R'] = dispersion_one['R']['Ad']
    model['Ad_L'] = dispersion_one['L']['Ad']

    # new output
    d = np.array(lines_model.pop(0).split()).astype('float')  # f.readline().split()).astype('float')
    vsv = np.array(lines_model.pop(0).split()).astype('float')  # np.array(f.readline().split()).astype('float')
    xi = np.array(lines_model.pop(0).split()).astype('float')  # np.array(f.readline().split()).astype('float')
    vp = np.array(lines_model.pop(0).split()).astype('float')  # np.array(f.readline().split()).astype('float')

    model['depth'] = d
    model['vsv'] = vsv
    model['xi'] = xi
    model['vp'] = vp

    dispersion_one['like_w'] = float(lines_model.pop(0))  # f.readline())
    model['like_w'] = dispersion_one['like_w']
    dispersion_one['cluster'] = cluster

    len_R = int(lines_model.pop(0))  # f.readline())
    dispersion_R_one = np.array(lines_model.pop(0).split()).astype(
        'float')  # np.array(f.readline().split()).astype('float')
    dispersion_one['R']['dispersion'] = dispersion_R_one

    len_L = int(lines_model.pop(0))  # f.readline())
    dispersion_L_one = np.array(lines_model.pop(0).split()).astype(
        'float')  # np.array(f.readline().split()).astype('float')

    dispersion_one['L']['dispersion'] = dispersion_L_one
    # except:
    #    return True, {}, {}

    if len(dispersion_L_one) != len_L or len(dispersion_R_one) != len_R or len(d) != npt_true or len(
            vsv) != npt_true or len(xi) != npt_true or len(vp) != npt_true:
        return True, {}, {}

    if len(dispersion_L_one) == 0 or len(dispersion_R_one) == 0 or len(d) == 0 or len(vsv) == 0 or len(xi) == 0 or len(
            vp) == 0:
        return True, {}, {}

    if model['Ad_R'] == 0 or model['Ad_L'] == 0:
        return True, {}, {}

    return True, dispersion_one, model


def get_model_ref(filename='Modified_PREM_GLOBAL.in'):
    '''
    reads the reference model, puts into a dict
    contains:
        npt_true: number of points
        depth: depth of interfaces
        vpv, vph, vsv, vsh: velocities

    Parameters
    ----------
    filename : str
        file containing the reference model. Default 'Model_PREM_SIMPLE.in'.

    Returns
    -------
    model_ref : dict
        contains the data of the reference model.

    '''

    file = filename

    model_ref = {}

    f = open(file, 'r')
    lines = f.readlines()
    f.close()

    data = lines[0].split()
    ntot = int(data[0])
    nic = int(data[1])
    noc = int(data[2])

    npt_true = ntot - noc

    rearth = float(lines[-1].split()[0])

    d = np.zeros((npt_true))
    vsv = np.zeros_like(d)
    vsh = np.zeros_like(d)
    vpv = np.zeros_like(d)
    vph = np.zeros_like(d)
    for i in range(npt_true):
        data = lines[ntot - i].split()
        d[i] = (rearth - float(data[0])) / 1000.
        vpv[i] = float(data[2]) / 1000.
        vsv[i] = float(data[3]) / 1000.
        vph[i] = float(data[6]) / 1000.
        vsh[i] = float(data[7]) / 1000.

    model_ref['npt_true'] = npt_true
    model_ref['depth'] = d
    model_ref['vpv'] = vpv
    model_ref['vph'] = vph
    model_ref['vsv'] = vsv
    model_ref['vsh'] = vsh

    return model_ref


def process_one_file_all(file, cluster, widenings_in, functions, params_dispersion, dispersion_ref,
                         dispersion_all, model_ref, outputs_all={}, file_probsum=None, thin_target=None):
    '''
    Processes one file, reading the podels in the file, applying the functions to them and stacking the results

    Parameters
    ----------
    input_dict : dict
        input dictionary containing the file name (file), the list of functions (functions),
        params_inversion, params_dispersion, dispersion_ref, dispersion_all and model_ref.

    Returns
    -------
    outputs : dict
        has one subdict for each function, called with the function name.
        Each of those has a subdict 'stack' and a subdict 'nostack' which will respectively stacked and kept as in the first model

    '''

    time_read = 0
    time_alpha = 0
    time_apply = 0

    numdis = dispersion_all['numdis']
    alpha_sum = np.zeros((numdis + 1))
    exp_sum = np.zeros((numdis + 1))

    t1 = time.time()

    f = open(file, 'r')
    lines = f.readlines()
    params_inversion = read_header(lines[:9])

    # thinning to get all files to same 'effective thinning'
    if thin_target is None:
        thin_true = 1
    else:
        thin_true = int(round(thin_target / params_inversion['thin']))

    widening = params_inversion['widening']

    t2 = time.time()
    time_read += t2 - t1

    t1 = time.time()

    if file_probsum is None:
        pass
    elif os.path.isfile(file_probsum):
        compute_probsum = False
        f_probsum = open(file_probsum, 'r')
        lines_probsum = f_probsum.readlines()
        f_probsum.close()
        i_line_probsum = 0
        line = lines_probsum[i_line_probsum]
        num_models_prop = {}
        while len(line.split()) == 3:
            data = line.split()
            cluster_tmp = int(data[0])
            widening_tmp = float(data[1])
            num_models_tmp = int(data[2])
            if cluster_tmp in num_models_prop:
                num_models_prop[cluster_tmp][widening_tmp] = num_models_tmp
            else:
                num_models_prop[cluster_tmp] = {}
                num_models_prop[cluster_tmp][widening_tmp] = num_models_tmp
            i_line_probsum += 1
            line = lines_probsum[i_line_probsum]

        compute_probsum = params_dispersion['num_models'] != num_models_prop
    else:
        compute_probsum = True
        dir = file_probsum.rpartition('/')[0]
        os.makedirs(dir, exist_ok=True)
    probsums_out = []

    t2 = time.time()
    time_alpha += t2 - t1

    numtot = 0
    i_model = 0
    valid_model = True
    maxlines = len(lines)
    while valid_model:
        if 9 + 11 * (i_model) > maxlines:
            valid_model = False
            continue
        lines_model = lines[9 + 11 * i_model:9 + 11 * (i_model + 1)]
        i_model += 1
        # thinning to get all files to same 'effective thinning'
        if (i_model // params_inversion['num_proc']) % thin_true != 0:
            continue
        t1 = time.time()
        valid_model, dispersion_one, model = read_model(lines_model, widening, cluster)
        #        if numtot>1000:
        #            valid_model=False
        if not valid_model:
            t2 = time.time()
            time_read += t2 - t1    
            continue
        elif len(dispersion_one.keys()) == 0:
            t2 = time.time()
            time_read += t2 - t1
            continue

        t2 = time.time()
        time_read += t2 - t1

        t1 = time.time()

        if file_probsum is None:
            probsum_prop = None
        elif compute_probsum:
            probsum_prop = None
        else:
            probsum_prop = float(lines_probsum[i_line_probsum])
            i_line_probsum += 1

        alpha, probsum = get_alpha_all(dispersion_one, widenings_in, params_dispersion, dispersion_ref,
                                       dispersion_all,
                                       num_models_cur=params_dispersion['num_models'][cluster][widening],
                                       probsum=probsum_prop)

        probsums_out.append(probsum)

        t2 = time.time()
        time_alpha += t2 - t1
        t1 = time.time()

        # apply functions
        for function in functions:
            function(file, model, model_ref, dispersion_one, params_dispersion, params_inversion, dispersion_all, alpha,
                     outputs=outputs_all[function.__name__])

        numtot += 1
        alpha_sum[:numdis] += np.minimum(np.zeros_like(alpha), alpha - dispersion_all['alpha_max'])

        exp_sum[:numdis] += np.exp(np.minimum(np.zeros_like(alpha), alpha - dispersion_all['alpha_max']))
        exp_sum[numdis] += 1.

        t2 = time.time()
        time_apply += t2 - t1
        # if numtot%10000==0:
        # print(numtot)

        # testing
        # if numtot==10:
        #     f.close()
        #     print('read: ',time_read)
        #     print('get alpha: ',time_alpha)
        #     print('apply: ',time_apply)

        #     return alpha_sum,numtot,exp_sum,exp_ssq
    f.close()

    if file_probsum is None:
        pass
    elif compute_probsum:
        dir = file_probsum.rpartition('/')[0]
        os.makedirs(dir, exist_ok=True)
        f_probsum = open(file_probsum, 'w')
        for cluster_probsum in params_dispersion['num_models']:
            for widening_probsum in params_dispersion['num_models'][cluster_probsum]:
                f_probsum.write(str(cluster_probsum) + ' ' + str(widening_probsum) + ' ' + str(
                    params_dispersion['num_models'][cluster_probsum][widening_probsum]) + '\n')
        for probsum in probsums_out:
            f_probsum.write(str(probsum) + '\n')
        f_probsum.close()

    print('read: ', time_read)
    print('get alpha: ', time_alpha)
    print('apply: ', time_apply)

    return alpha_sum, numtot, exp_sum


def get_alpha_norm(dispersion_one, widenings, params_dispersion, dispersion_ref, dispersion_all):
    '''
    Gets alpha for one dispersion curve

    Parameters
    ----------
    dispersion_one : dict
        dispersion curve for one model
    params_dispersion : dict

    dispersion_ref : dict

    dispersion_all : dict


    Returns
    -------
    alpha_ref : float
        alpha for the reference dispersion curve.
    alpha : numpy array, length numdis
        alpha for all cluster members.

    '''
    clusters = params_dispersion['clusterlist']

    ndatad_R = params_dispersion['R']['ndatad']
    ndatad_L = params_dispersion['L']['ndatad']
    dispersion_R = dispersion_all['R']['dispersion']
    dispersion_L = dispersion_all['L']['dispersion']
    dispersion_R_one = dispersion_one['R']['dispersion']
    dispersion_L_one = dispersion_one['L']['dispersion']
    Ad_R = dispersion_one['R']['Ad']
    Ad_L = dispersion_one['L']['Ad']
    error_R = dispersion_all['R']['error']
    error_L = dispersion_all['L']['error']
    error_R_ref_sum = dispersion_ref['R']['error_sum']
    error_L_ref_sum = dispersion_ref['L']['error_sum']
    numdis = dispersion_all['numdis']

    like_alt = -np.sum(np.divide(
        np.square(dispersion_R - np.transpose(np.tile(dispersion_R_one, (numdis, 1)))),
        2 * Ad_R ** 2 *
        np.square(
            np.multiply(
                error_R, dispersion_R / 100.)))
        , axis=0)
    like_alt -= np.sum(np.divide(
        np.square(dispersion_L - np.transpose(np.tile(dispersion_L_one, (numdis, 1)))),
        2 * Ad_L ** 2 *
        np.square(
            np.multiply(
                error_L, dispersion_L / 100.)))
        , axis=0)
    # static shift
    # like_alt-=error_R_sum
    # like_alt-=error_L_sum
    like_alt -= ndatad_R * np.log(Ad_R)
    like_alt -= ndatad_L * np.log(Ad_L)

    probsum = 0.
    shift = 0
    while probsum == 0.:
        try:
            numprob = 0.
            probsum = 0.
            for i in range(len(clusters)):
                for widening in widenings:

                    dispersion_R_ref = dispersion_ref['R']['dispersion'][i, :]
                    dispersion_L_ref = dispersion_ref['L']['dispersion'][i, :]

                    error_R_ref = dispersion_ref['R']['error'][i, :]
                    error_L_ref = dispersion_ref['L']['error'][i, :]

                    if np.shape(dispersion_R_ref) == (0,):
                        print('here')
                    if np.shape(dispersion_R_one) == (0,):
                        print('here2')
                        print(dispersion_R_one)

                    # print('new')
                    like_one = -np.sum(np.divide(np.square(dispersion_R_ref - dispersion_R_one),
                                                 2 * Ad_R ** 2 * np.square(
                                                     np.multiply(error_R_ref / 100., dispersion_R_ref))), axis=0)
                    # print(like_one)
                    like_one -= np.sum(np.divide(np.square(dispersion_L_ref - dispersion_L_one),
                                                 2 * Ad_L ** 2 * np.square(
                                                     np.multiply(error_L_ref / 100., dispersion_L_ref))), axis=0)
                    # print(like_one)
                    like_one -= error_R_ref_sum[i] + 150  # static shift to avoid exponential explosion
                    # print(like_one)
                    like_one -= error_L_ref_sum[i] + 150  # static shift to avoid exponential explosion
                    # print(like_one)
                    like_one -= ndatad_R * np.log(Ad_R)
                    # print(like_one)
                    like_one -= ndatad_L * np.log(Ad_L)
                    # print(like_one)
                    like_one /= widening
                    # print(like_one)
                    # like_one+=1000
                    # print(like_one)

                    like_one += shift

                    like_one -= np.log(params_dispersion['num_models'][clusters[i]][widening])

                    probsum += np.exp(like_one)

                    numprob += 1
            # print(probsum)
            if probsum == 0.:
                raise Probsum1
        except RuntimeWarning:
            probsum = 0.
            shift -= 100
            print('too big probsum')
        except Probsum1:
            probsum = 0.
            shift += 1000
            # print('too small probsum')
    # print('done')
    probsum /= numprob
    # if probsum==0.:
    #     alpha=np.ones_like(like_alt)*(-float('inf'))
    # else:
    alpha = like_alt - (np.log(probsum) - shift)

    return alpha


def count_models(directories, widenings, params_dispersion, comm, thin_target=None):
    size = comm.Get_size()
    rank = comm.Get_rank()
    rank2 = rank

    clusters = params_dispersion['clusterlist']

    files_all = []
    num_models_dict = {}
    # initiate stuff
    for i, cluster in enumerate(clusters):
        num_models_dict[cluster] = {}
        for directory in params_dispersion['directories'][cluster]:
            for widening in widenings:
                num_models_dict[cluster][widening] = 0
                # list files
                files_all_tmp = glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening)
                files_all_tmp.sort()
                files_all.extend(files_all_tmp)

    whodoeswhat = {}
    for i in range(len(files_all)):
        whodoeswhat[files_all[i]] = i % size

    num_models_dict = {}
    for i in np.unique(clusters):
        num_models_dict[i] = {}
        for widening in widenings:
            num_models_dict[i][widening] = 0

    for i, cluster in enumerate(clusters):
        for widening in widenings:
            files_all_tmp = []
            for directory in params_dispersion['directories'][cluster]:
            
                files_all_tmp.extend(glob.glob(directory + '/All_models_processed_prepare_*%4.2f.out' % widening))
                # files_all_tmp=['OUT_CLUSTER0_START/OUT/All_models_processed_prepare_0000_ 1.00.out']
            files_all_tmp.sort()
            files = []
            for i in range(len(files_all_tmp)):
                if i % size == rank2:
                    files.append(files_all_tmp[i])

            # if more cores than files: go to the next iteration of the loop
            # but keep core 0 in, it's the core that gathers everything
            computing = []
            if len(files) == 0:
                rank2 -= len(files_all_tmp)
                if rank != 0:
                    continue
            for f in files_all_tmp:
                try:
                    computing.append(whodoeswhat[f])
                except:
                    pass

            num_models_tmp = 0

            if rank == 0:
                print(computing)

            # loop over the files
            for i in range(len(files)):
                file = files[i]

                # process one file
                print('counting for file', file, 'rank', rank)
                num_models_prop = count_one_file(files[i], thin_target=thin_target)
                num_models_tmp += num_models_prop

            computing[:] = [x for x in computing if x != 0]
            # merge results between cores on core 0
            if rank in np.unique(computing):
                comm.Ssend([np.array([num_models_tmp]), MPI.INT], dest=0,
                            tag=3 + int(widening) * 100 + int(cluster) * 1000)
            if rank == 0:
                # recieve and merge data from all other cores
                for i in np.unique(computing):
                    num_models_prop = np.zeros_like(num_models_tmp)
                    comm.Recv([num_models_prop, MPI.INT], i, 3 + int(widening) * 100 + int(cluster) * 1000)
                    num_models_tmp += num_models_prop
                    print('recieved', i)

                num_models_dict[cluster][widening] += num_models_tmp

            if rank != 0 or len(files) != 0:
                rank2 += size * len(files) - len(files_all_tmp)

    num_models_dict = comm.bcast(num_models_dict, root=0)
    params_dispersion['num_models'] = num_models_dict

    return


def count_one_file(file, thin_target=None):
    f = open(file, 'r')
    lines = f.readlines()

    params_inversion = read_header(lines[:9])

    # thinning to get all files to same 'effective thinning'
    if thin_target is None:
        thin_true = 1
    else:
        thin_true = int(round(thin_target / params_inversion['thin']))

    if not params_inversion:
        return

    valid_model = True
    cnt = 0
    i_model = 0
    maxlines = len(lines)
    print('decimation: ',thin_true)
    while valid_model:
        if 9 + 11 * (i_model) > maxlines:
            valid_model = False
            continue
        lines_model = lines[9 + 11 * i_model:9 + 11 * (i_model + 1)]
        i_model += 1
        # thinning to get all files to same 'effective thinning'
        if (i_model // params_inversion['num_proc']) % thin_true != 0:
            continue
        valid_model, dispersion_one, model = read_model(lines_model, 1, 0)
        if not valid_model:
            continue
        elif len(dispersion_one.keys()) == 0:
            continue

        cnt += 1
    f.close()
    return cnt
