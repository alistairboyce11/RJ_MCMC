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
import matplotlib.pyplot as plt
import concurrent.futures
import multiprocessing
import sys
from mpi4py import MPI


#########################################################
# Read all true dispersion curves
#########################################################

def get_dispersion(directory, filename='dispersion_all.in'):
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

    dispersion_ref={}
    dispersion_ref['R']={}
    dispersion_ref['L']={}
    dispersion_all={}
    dispersion_all['R']={}
    dispersion_all['L']={}
    params_dispersion={}
    params_dispersion['R']={}
    params_dispersion['L']={}
    file=open(directory+'/'+filename,'r')

    # for simultaneous inversions
    numdis=int(file.readline())
    dispersion_all['numdis']=numdis


    dispersion_all['lat']=np.array(file.readline().split()).astype('float')
    dispersion_all['lon']=np.array(file.readline().split()).astype('float')
    dispersion_all['cluster']=np.array(file.readline().split()).astype('int')

    # Read rayleigh data
    ndatad_R=int(file.readline())
    params_dispersion['R']['ndatad']=ndatad_R
    mode_rayl=int(file.readline())

    mode_R=np.zeros((ndatad_R))
    period_R=np.zeros((ndatad_R))
    dispersion_R=np.zeros((ndatad_R,numdis))
    error_R=np.zeros((ndatad_R,numdis))
    dispersion_R_ref=np.zeros((ndatad_R))
    error_R_ref=np.zeros((ndatad_R))
    j=0
    for i in range(mode_rayl):
        file.readline()
        ind_k=int(file.readline())+1
        for k in range(ind_k):
            data=file.readline().split()
            mode_R[j]=int(data[0])
            period_R[j]=float(data[1])
            dispersion_R_ref[j]=float(data[2])
            error_R_ref[j]=float(data[3])

            #for simultaneous inversions
            for l in range(numdis):
                data=file.readline().split()
                dispersion_R[j,l]=float(data[0])
                error_R[j,l]=float(data[1])
            j+=1

    dispersion_ref['R']['dispersion']=dispersion_R_ref
    dispersion_ref['R']['error']=error_R_ref
    dispersion_all['R']['dispersion']=dispersion_R
    dispersion_all['R']['error']=error_R
    params_dispersion['R']['periods']=period_R
    params_dispersion['R']['modes']=mode_R

    # read love data
    ndatad_L=int(file.readline())
    params_dispersion['L']['ndatad']=ndatad_L
    mode_love=int(file.readline())

    mode_L=np.zeros((ndatad_L))
    period_L=np.zeros((ndatad_L))
    dispersion_L=np.zeros((ndatad_L,numdis))
    error_L=np.zeros((ndatad_L,numdis))
    dispersion_L_ref=np.zeros((ndatad_L))
    error_L_ref=np.zeros((ndatad_L))
    j=0
    for i in range(mode_love):
        file.readline()
        ind_k=int(file.readline())+1
        for k in range(ind_k):
            data=file.readline().split()
            mode_L[j]=int(data[0])
            period_L[j]=float(data[1])
            dispersion_L_ref[j]=float(data[2])
            error_L_ref[j]=float(data[3])

            # for simultaneous inversions
            for l in range(numdis):
                data=file.readline().split()
                dispersion_L[j,l]=float(data[0])
                error_L[j,l]=float(data[1])
            j+=1

    dispersion_ref['L']['dispersion']=dispersion_L_ref
    dispersion_ref['L']['error']=error_L_ref
    dispersion_all['L']['dispersion']=dispersion_L
    dispersion_all['L']['error']=error_L
    params_dispersion['L']['periods']=period_L
    params_dispersion['L']['modes']=mode_L

    file.close()

    return params_dispersion, dispersion_ref, dispersion_all

def get_metadata(directory, prepare=False,widening=1.0):
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

    params_inversion={}

    if prepare:
        files=glob.glob(directory+'/All_models_prepare*%4.2f.out'%widening)
    else:
        files=glob.glob(directory+'/All_models_invert*.out')

    files.sort()

    ind=0
    file=files[ind]

    f=open(file,'r')



    d=f.readline() # contains rank of processor, number of file for the processor and number of models stored at most. Not used or tested currently
    while not d: # if some of the files are empty, take the first one that isn't
        f.close()
        ind+=1
        file=files[ind]

        f=open(file,'r')
        d=f.readline()
        print(ind)

    params_inversion['everyall']=int(float(d.split()[2]))
    data=f.readline().split() # contains a forth information, the thinning. Not used currently, but can be added easily
    params_inversion['burn-in']=float(data[0])
    params_inversion['nsample']=float(data[1])
    params_inversion['widening']=float(data[2])
    params_inversion['thin']=float(data[3])
    data=f.readline().split()
    params_inversion['d_min']=float(data[0])
    params_inversion['d_max']=float(data[1])
    params_inversion['width_vsv']=float(f.readline())
    data=f.readline().split()
    params_inversion['xi_min']=float(data[0])
    params_inversion['xi_max']=float(data[1])
    data=f.readline().split()
    params_inversion['vpvs']=float(data[0])
    params_inversion['vpvs_min']=float(data[1])
    params_inversion['vpvs_max']=float(data[2])
    data=f.readline().split()
    params_inversion['Ad_R_min']=float(data[0])
    params_inversion['Ad_R_max']=float(data[1])
    data=f.readline().split()
    params_inversion['Ad_L_min']=float(data[0])
    params_inversion['Ad_L_max']=float(data[1])
    data=f.readline().split()
    params_inversion['milay']=int(data[0])
    params_inversion['malay']=int(data[1])
    f.close()

    return params_inversion

def get_alpha_max(comm,directory,params_inversion, params_dispersion,dispersion_ref,dispersion_all,maxpercent=0.,prepare=False,widening=1.):
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
    prepare : bool, optional
        says if we are looking at results from the preprocessing or the main inversion. The default is False.
    widening : float, optional
        if prepare, which part of the preprocessing are we looking at? The default is 1..

    Returns
    -------
    None.
    Modifies dispersion_all and dispersion_ref in place, adding 'alpha_max' keys to them with the corresponding values

    '''

    if prepare:
        files_all=glob.glob(directory+'/All_models_prepare*%4.2f.out'%widening)
    else:
        files_all=glob.glob(directory+'/All_models_invert*.out')

    files_all.sort()
    files_all=files_all

    rank = comm.Get_rank()

    size = comm.Get_size()
    files=[]
    for i in range(len(files_all)):
        if i%size==rank:
            files.append(files_all[i])

    numdis=dispersion_all['numdis']

    if maxpercent==0.:
        alpha_ref_max=float('-inf')
        alpha_max=np.ones((numdis))*float('-inf')
    else:
        alpha_ref_max=0
        alpha_max=np.zeros((numdis))

    print(rank,files)

    for file in files:

        # input_dict={}
        # input_dict['file']=file
        # input_dict['params_inversion']=params_inversion
        # input_dict['params_dispersion']=params_dispersion
        # input_dict['dispersion_ref']=dispersion_ref
        # input_dict['dispersion_all']=dispersion_all
        # input_dict['maxpercent']=maxpercent
        # input_dict['prepare']=prepare

        alpha_ref_max_prop,alpha_max_prop=process_file_alphamax(file,params_inversion,params_dispersion,dispersion_ref,dispersion_all,maxpercent,prepare)

        if maxpercent==0.:
            alpha_ref_max=max(alpha_ref_max_prop,alpha_ref_max)
            alpha_max=np.maximum.reduce([alpha_max,alpha_max_prop])
        else:
            alpha_ref_max+=alpha_ref_max_prop/len(files_all)
            alpha_max+=alpha_max_prop/len(files_all)


    #print(rank,alpha_ref_max)

    # for i in range(size)[1:]:
    #     if rank==i:
    #         comm.send(obj=(alpha_ref_max,alpha_max),dest=0,tag=rank)
    #     elif rank==0:
    #         print(i)
    #         alpha_ref_max_prop,alpha_max_prop=comm.recv(source=i,tag=rank)
    #         #alpha_max_prop=comm.recv(i,1)
    #         if maxpercent==0.:
    #             alpha_ref_max=max(alpha_ref_max_prop,alpha_ref_max)
    #             alpha_max=np.maximum.reduce([alpha_max,alpha_max_prop])
    #         else:
    #             alpha_ref_max+=alpha_ref_max_prop/len(files_all)
    #             alpha_max+=alpha_max_prop/len(files_all)
    #         print(alpha_ref_max)

    #alpha_ref_max=rank**2

    if rank!=0:
        #print('sending',rank)
        comm.send(obj=(alpha_ref_max,alpha_max),dest=0,tag=1)

        #print('sent',rank)
    elif rank==0:

        for i in range(size)[1:]:
            #print('recieving',i)
            (alpha_ref_max_prop,alpha_max_prop)=comm.recv(source=i,tag=1)
            #print('recieved',i)
            #alpha_max_prop=comm.recv(i,1)
            if maxpercent==0.:
                alpha_ref_max=max(alpha_ref_max_prop,alpha_ref_max)
                alpha_max=np.maximum.reduce([alpha_max,alpha_max_prop])
            else:
                alpha_ref_max+=alpha_ref_max_prop
                alpha_max+=alpha_max_prop


    #if rank==0:
    #time.sleep(10)

    alpha_max = comm.bcast(alpha_max, root=0)
    alpha_ref_max = comm.bcast(alpha_ref_max, root=0)

    #print(rank,alpha_ref_max)
    dispersion_all['alpha_max']=alpha_max
    dispersion_ref['alpha_max']=alpha_ref_max

def process_file_alphamax(file,params_inversion,params_dispersion,dispersion_ref,dispersion_all,maxpercent,prepare):
    '''processes one file to get alphamax
    input: input_dict
    input_dict contains:
        file: filename of the file to process
        params_inversion


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

    # file=input_dict['file']
    # params_inversion=input_dict['params_inversion']
    # params_dispersion=input_dict['params_dispersion']
    # dispersion_ref=input_dict['dispersion_ref']
    # dispersion_all=input_dict['dispersion_all']
    # maxpercent=input_dict['maxpercent']
    # prepare=input_dict['prepare']

    # reading the file
    f=open(file,'r')
    print(file)

    numdis=dispersion_all['numdis']
    alpha_ref_max=float('-inf')
    alpha_max=np.ones((numdis))*float('-inf')

    if prepare:
        num_file=int(float(file.split('/')[-1].split('_')[4]))

    else:
        num_file=int(float(file.split('/')[-1].split('_')[-1].split('.')[0]))

    num_models=int(min(params_inversion['everyall'],params_inversion['nsample']//params_inversion['thin']-num_file*params_inversion['everyall']))

    alpha_ref_all=np.ones((int(num_models*maxpercent/100)))*float('-inf')
    alpha_all=np.ones((numdis,int(num_models*maxpercent/100)))*float('-inf')

    fline=f.readline()
    if not fline:
        return
    widening=float(f.readline().split()[2])
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    num_model=0

    line=f.readline()
    while line:
        #print(num_models)
        dispersion_one={}
        dispersion_one['R']={}
        dispersion_one['L']={}
        data=line.split()
        npt_true=int(data[0])
        data=f.readline().split()
        dispersion_one['R']['Ad']=float(data[0])
        dispersion_one['L']['Ad']=float(data[1])

        for i in range(npt_true):
            f.readline()

        dispersion_one['like_w']=float(f.readline())
        dispersion_one['widening']=widening

        f.readline()
        dispersion_one['R']['dispersion']=np.array(f.readline().split()).astype('float')

        f.readline()
        dispersion_one['L']['dispersion']=np.array(f.readline().split()).astype('float')

        alpha_ref,alpha=get_alpha(dispersion_one,params_dispersion,dispersion_ref,dispersion_all)

        if maxpercent==0.:
            alpha_ref_max=max(alpha_ref,alpha_ref_max)
            alpha_max=np.maximum.reduce([alpha_max,alpha])

        else:
            min_indices=np.argmin(alpha_all,axis=1)
            alpha_all[np.arange(numdis),min_indices]=np.maximum(alpha,np.amin(alpha_all,axis=1))
            if alpha_ref>np.amin(alpha_ref_all):
                alpha_ref_all[np.argmin(alpha_ref_all)]=alpha_ref

        num_model+=1
        line=f.readline()

        #if num_model==10:
        #    break
    f.close()

    if maxpercent>0.:
        alpha_max=np.amin(alpha_all,axis=1)
        alpha_ref_max=np.amin(alpha_ref_all)

    return alpha_ref_max,alpha_max

def get_alpha(dispersion_one,params_dispersion,dispersion_ref,dispersion_all):
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
    widening=dispersion_one['widening']
    ndatad_R=params_dispersion['R']['ndatad']
    ndatad_L=params_dispersion['L']['ndatad']
    dispersion_R_ref=dispersion_ref['R']['dispersion']
    dispersion_L_ref=dispersion_ref['L']['dispersion']
    dispersion_R=dispersion_all['R']['dispersion']
    dispersion_L=dispersion_all['L']['dispersion']
    dispersion_R_one=dispersion_one['R']['dispersion']
    dispersion_L_one=dispersion_one['L']['dispersion']
    Ad_R=dispersion_one['R']['Ad']
    Ad_L=dispersion_one['L']['Ad']
    error_R_ref=dispersion_ref['R']['error']
    error_L_ref=dispersion_ref['L']['error']
    error_R=dispersion_all['R']['error']
    error_L=dispersion_all['L']['error']
    like_w=dispersion_one['like_w']
    numdis=dispersion_all['numdis']

    like_alt_R=0
    like_alt_L=0
    #print(np.shape(dispersion_R_ref))
    #print(np.shape(dispersion_R_one))
    #print(np.shape(error_R_ref))
    if np.shape(dispersion_R_ref)==(0,):
        print('here')
    if np.shape(dispersion_R_one)==(0,):
        print('here2')
        print(dispersion_R_one)
    #print(np.shape(dispersion_R_ref-dispersion_R_one))
    #print(np.shape(np.square(dispersion_R_ref-dispersion_R_one)))
    #print(np.shape(2*Ad_R**2*np.square(error_R_ref)))
    #print(np.shape(np.square(dispersion_R_ref-dispersion_R_one)),np.shape(2*Ad_R**2*np.square(error_R_ref)))
    like_alt_R=np.sum(np.divide(np.square(dispersion_R_ref-dispersion_R_one),2*Ad_R**2*np.square(error_R_ref)))
    like_alt_L=np.sum(np.divide(np.square(dispersion_L_ref-dispersion_L_one),2*Ad_L**2*np.square(error_L_ref)))
    like_alt=like_alt_R+like_alt_L

    alpha_ref=ndatad_R*(1/widening-1)*np.log(Ad_R)+ndatad_L*(1/widening-1)*np.log(Ad_L)+like_w-like_alt

    alpha=np.zeros((numdis))

    like_alt_R=np.sum(np.divide(np.square(dispersion_R-np.transpose(np.tile(dispersion_R_one,(numdis,1)))),2*Ad_R**2*np.square(error_R)),axis=0)

    like_alt_L=np.sum(np.divide(np.square(dispersion_L-np.transpose(np.tile(dispersion_L_one,(numdis,1)))),2*Ad_L**2*np.square(error_L)),axis=0)

    alpha=ndatad_R*(1/widening-1)*np.log(Ad_R)+ndatad_L*(1/widening-1)*np.log(Ad_L)+like_w-(like_alt_R+like_alt_L)

    return alpha_ref,alpha

def get_model_ref(filename='Model_PREM_SIMPLE.in'):
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

    file=filename

    model_ref={}

    f=open(file,'r')
    lines=f.readlines()
    f.close()

    data=lines[0].split()
    ntot=int(data[0])
    nic=int(data[1])
    noc=int(data[2])

    npt_true=ntot-noc

    rearth=float(lines[-1].split()[0])

    d=np.zeros((npt_true))
    vsv=np.zeros_like(d)
    vsh=np.zeros_like(d)
    vpv=np.zeros_like(d)
    vph=np.zeros_like(d)
    for i in range(npt_true):
        data=lines[ntot-i].split()
        d[i]=(rearth-float(data[0]))/1000.
        vpv[i]=float(data[2])
        vsv[i]=float(data[3])
        vph[i]=float(data[6])
        vsh[i]=float(data[7])

    model_ref['npt_true']=npt_true
    model_ref['depth']=d
    model_ref['vpv']=vpv
    model_ref['vph']=vph
    model_ref['vsv']=vsv
    model_ref['vsh']=vsh

    return model_ref


def apply_stuff(comm,directory,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,prepare=False,widening=1.):
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

    if prepare:
        files_all=glob.glob(directory+'/All_models_prepare*%4.2f.out'%widening)
    else:
        files_all=glob.glob(directory+'/All_models_invert*.out')

    files_all.sort()
    files_all=files_all

    rank = comm.Get_rank()

    size = comm.Get_size()
    files=[]
    for i in range(len(files_all)):
        if i%size==rank:
            files.append(files_all[i])

    files.sort()


    outputs_all={}
    #for function in functions:
    #    outputs_all[function.__name__]={}
    #    outputs_all[function.__name__]['stack']={}
    #    outputs_all[function.__name__]['nostack']={}

    #parallel processing needs a list of single inputs, so we put all input into a dict and create a list of dicts
    #input_dicts=[]
    # for file in files:
    #     input_dict={}
    #     input_dict['file']=file
    #     input_dict['functions']=functions
    #     input_dict['params_inversion']=params_inversion
    #     input_dict['params_dispersion']=params_dispersion
    #     input_dict['dispersion_ref']=dispersion_ref
    #     input_dict['dispersion_all']=dispersion_all
    #     input_dict['model_ref']=model_ref
    #     input_dicts.append(input_dict)

    # # Parallel processing
    # with concurrent.futures.ProcessPoolExecutor(max_workers=28) as executor:

    #     results = executor.map(process_one_file, input_dicts)

    #     # collapsing of results
    #     for res in results:
    #         for function in res:
    #             if function in outputs_all:
    #                 for key in res[function]['stack']: # stack are stacked
    #                     outputs_all[function]['stack'][key]+=res[function]['stack'][key]
    #             else:
    #                 outputs_all[function]={}
    #                 outputs_all[function]['stack']={}
    #                 outputs_all[function]['nostack']={}
    #                 for key in res[function]['stack']:
    #                     outputs_all[function]['stack'][key]=res[function]['stack'][key]
    #                 for key in res[function]['nostack']:
    #                     outputs_all[function]['nostack'][key]=res[function]['nostack'][key] # nostack are kept as in the first one

    for file in files:
        output=process_one_file(file,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref)
        print('done',rank,file)
        for function in output:
            if function in outputs_all:
                for key in output[function]['stack']: # stack are stacked
                    outputs_all[function]['stack'][key]+=output[function]['stack'][key]
            else:
                outputs_all[function]={}
                outputs_all[function]['stack']={}
                outputs_all[function]['nostack']={}

                for key in output[function]['stack']:
                    outputs_all[function]['stack'][key]=output[function]['stack'][key]
                for key in output[function]['nostack']:
                    outputs_all[function]['nostack'][key]=output[function]['nostack'][key]
                for key in output[function]:
                    if key=='stack' or key=='nostack':
                        continue
                    outputs_all[function][key]=output[function][key]
    for function in outputs_all:
        for key in outputs_all[function]['stack']:

            if rank!=0:
                print(function,key)
                if type(outputs_all[function]['stack'][key])==float:
                    outputs_all[function]['stack'][key]=np.array([outputs_all[function]['stack'][key]])
                #outputs_all[function]['stack'][key]=outputs_all[function]['stack'][key].astype('f')
                comm.Send([outputs_all[function]['stack'][key],MPI.FLOAT],dest=0,tag=3)

            if rank==0:
                for i in range(size)[1:]:
                    if type(outputs_all[function]['stack'][key])==float:
                        outputs_all_prop=np.array([0.])
                    else:
                        outputs_all_prop=np.zeros_like(outputs_all[function]['stack'][key])
                    comm.Recv([outputs_all_prop,MPI.FLOAT],i,3)
                    if type(outputs_all[function]['stack'][key])==float:
                        outputs_all_prop=outputs_all_prop[0]
                    outputs_all[function]['stack'][key]+=outputs_all_prop

    #         for function in outputs_all_prop:
    #             for key in outputs_all_prop[function]['stack']: # stack are stacked
    #                     outputs_all[function]['stack'][key]+=outputs_all_prop[function]['stack'][key]
    if rank==0:
        for function in outputs_all:
            if 'meansum' in outputs_all[function].keys():
            # for each member of the cluster
                for i in range(len(outputs_all[function]['meansum'])): # calculate average
                    key_mean=outputs_all[function]['meansum'][i]
                    key_sum=outputs_all[function]['meansum_weights'][i]

                    weights=np.transpose(np.tile(outputs_all[function]['stack'][key_sum],(outputs_all[function]['nostack']['ndata'],1))) # sum of weights

                    outputs_all[function]['stack'][key_mean]=np.divide(outputs_all[function]['stack'][key_mean],weights)

            if 'meansum_sq' in outputs_all[function].keys():
                for i in range(len(outputs_all[function]['meansum_sq'])): # calculate standart deviation
                    key_sum=outputs_all[function]['meansum_sq_weights'][i]

                    key_mean_sq=outputs_all[function]['meansum_sq'][i]

                    key_mean_shift=outputs_all[function]['meansum_sq_shift'][i]

                    print(key_sum,key_mean_sq,key_mean_shift)

                    weights=np.transpose(np.tile(outputs_all[function]['stack'][key_sum],(outputs_all[function]['nostack']['ndata'],1))) # get weights
                    #weights=outputs_all[function]['stack'][key_sum] # get weights

                    print(np.shape(weights),np.shape(outputs_all[function]['stack'][key_sum]),outputs_all[function]['nostack']['ndata'],np.shape(np.tile(outputs_all[function]['stack'][key_sum],(outputs_all[function]['nostack']['ndata'],1))),np.shape(outputs_all[function]['stack'][key_mean_sq]))

                    mean_sq=np.divide(outputs_all[function]['stack'][key_mean_sq],weights) # mean of squares
                    mean_shifted=np.square(np.divide(outputs_all[function]['stack'][key_mean_shift],weights)) # square of means, shifted by reference model

                    outputs_all[function]['stack'][key_mean_sq]=np.sqrt(mean_sq-mean_shifted)

                    del outputs_all[function]['stack'][key_mean_shift]

            if 'meansum_ref' in outputs_all[function].keys():
                # for the cluster centre
                for i in range(len(outputs_all[function]['meansum_ref'])):
                    key_mean=outputs_all[function]['meansum_ref'][i]
                    key_sum=outputs_all[function]['meansum_ref_weights'][i]

                    weights=outputs_all[function]['stack'][key_sum] # sum of weights

                    outputs_all[function]['stack'][key_mean]=outputs_all[function]['stack'][key_mean]/weights # mean
            if 'meansum_sq_ref' in outputs_all[function].keys():
                for i in range(len(outputs_all[function]['meansum_sq_ref'])):
                    key_sum=outputs_all[function]['meansum_sq_ref_weights'][i]

                    key_mean_sq=outputs_all[function]['meansum_sq_ref'][i]

                    key_mean_shift=outputs_all[function]['meansum_sq_ref_shift'][i]

                    weights=outputs_all[function]['stack'][key_sum] # sum of weights

                    mean_sq=outputs_all[function]['stack'][key_mean_sq]/weights # mean of squares
                    mean_shifted=(outputs_all[function]['stack'][key_mean_shift]/weights)**2 # square of mean, shifted

                    outputs_all[function]['stack'][key_mean_sq]=np.sqrt(mean_sq-mean_shifted)

                    del outputs_all[function]['stack'][key_mean_shift]

    # outputs_all = comm.bcast(outputs_all, root=0)

    return outputs_all



def process_one_file(file,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref):
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

    time_read=0
    time_alpha=0
    time_apply=0
    # file=input_dict['file']
    # functions=input_dict['functions']
    # params_inversion=input_dict['params_inversion']
    # params_dispersion=input_dict['params_dispersion']
    # dispersion_ref=input_dict['dispersion_ref']
    # dispersion_all=input_dict['dispersion_all']
    # model_ref=input_dict['model_ref']
    outputs={}

    numtot=0

    # read file
    f=open(file,'r')
    print(file)

    #file2='/'.join(file.split('/')[:-1])+'/alpha_'+file.split('/')[-1]
    #
    #if os.path.isfile(file2):
    #    g=open(file2,'r')
    #    reading=True
    #else:
    #    g=open(file2,'w')
    #    reading=False

    fline=f.readline()
    if not fline:
        print('file ' +file + ' empty')
        return
    widening=float(f.readline().split()[2])
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    line=f.readline()
    while line:

        t1=time.time()
        dispersion_one={}
        dispersion_one['R']={}
        dispersion_one['L']={}

        model={}

        data=line.split()
        npt_true=int(data[0])
        model['npt_true']=npt_true
        model['npt']=int(data[1])
        model['npt_ani']=int(data[2])
        data=f.readline().split()
        dispersion_one['R']['Ad']=float(data[0])
        dispersion_one['L']['Ad']=float(data[1])
        model['Ad_R']=dispersion_one['R']['Ad']
        model['Ad_L']=dispersion_one['L']['Ad']

        d=np.zeros(npt_true)
        vsv=np.zeros_like(d)
        xi=np.zeros_like(d)
        vp=np.zeros_like(d)

        for i in range(npt_true):
            data=f.readline().split()
            d[i]=float(data[0])
            vsv[i]=float(data[1])
            xi[i]=float(data[2])
            vp[i]=float(data[3])
        model['depth']=d
        model['vsv']=vsv
        model['xi']=xi
        model['vp']=vp

        dispersion_one['like_w']=float(f.readline())
        model['like_w']=dispersion_one['like_w']
        dispersion_one['widening']=widening

        f.readline()
        dispersion_R_one=np.array(f.readline().split()).astype('float')
        dispersion_one['R']['dispersion']=dispersion_R_one

        f.readline()
        dispersion_L_one=np.array(f.readline().split()).astype('float')

        dispersion_one['L']['dispersion']=dispersion_L_one

        t2=time.time()
        time_read+=t2-t1

        t1=time.time()


        #if not reading:
            # calculating alphas
        alpha_ref,alpha=get_alpha(dispersion_one,params_dispersion,dispersion_ref,dispersion_all)


       #     g.write(str(alpha_ref)+'\n')
       #     g.write(' '.join(map(str,alpha))+'\n')
       # else:
       #     alpha_ref=float(g.readline())
       #     alpha=np.array(g.readline().split()).astype('float')

        t2=time.time()
        time_alpha+=t2-t1

        t1=time.time()

        # apply functions
        for function in functions:
            output=function(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=(numtot==0))

            # stack outputs
            if function.__name__ in outputs:
                for key in output['stack']:

                    outputs[function.__name__]['stack'][key]+=output['stack'][key]
            else:
                outputs[function.__name__]={}
                outputs[function.__name__]['stack']={}
                outputs[function.__name__]['nostack']={}
                for key in output['stack']:
                    outputs[function.__name__]['stack'][key]=output['stack'][key]
                for key in output['nostack']:
                    outputs[function.__name__]['nostack'][key]=output['nostack'][key]
                for key in output:
                    if key=='stack' or key=='nostack':
                        continue
                    outputs[function.__name__][key]=output[key]
        line=f.readline()
        #print(numtot)
        numtot+=1

        t2=time.time()
        time_apply+=t2-t1

        #if numtot==10:
        #    break

    # normalize outputs
    #for function in functions:
    #    for key in outputs[function.__name__]['stack']:
    #        outputs[function.__name__]['stack'][key]/=numtot
    f.close()

    #g.close()

    print('read: ',time_read)
    print('get alpha: ',time_alpha)
    print('apply: ',time_apply)

    return outputs

def create_posterior(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=True):
    '''
    example of a function to be applied to the data
    creates a 2D histogram of vsv for the reference and each member of the cluster
    The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

    Parameters
    ----------
    model : dict

    model_ref : dict

    dispersion_one : dict

    params_dispersion : dict

    params_inversion : dict

    dispersion_ref : dict

    dispersion_all : dict

    alpha_ref : float

    alpha : numpy array, length numdis

    first : bool, optional
        whether this model is the first of its file. May provide a minor speed boost. The default is True.

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    ndatad=50
    ndatav=50
    outputs={}
    outputs['stack']={}
    outputs['nostack']={}
    vsv_ref=np.zeros((ndatad,ndatav))
    xi_ref=np.zeros_like(vsv_ref)
    vp_ref=np.zeros_like(vsv_ref)

    vsv_wide=np.zeros((ndatad,ndatav))
    xi_wide=np.zeros_like(vsv_wide)
    vp_wide=np.zeros_like(vsv_wide)

    numdis=dispersion_all['numdis']
    vsv_all=np.zeros((ndatad,ndatav,numdis))
    xi_all=np.zeros_like(vsv_all)
    vp_all=np.zeros_like(vsv_all)

    depths=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)
    vels_vsv=np.linspace(np.amin(model_ref['vsv']*(1-params_inversion['width_vsv'])),
                     np.amax(model_ref['vsv']*(1+params_inversion['width_vsv'])),
                     ndatav)
    vels_vp=np.linspace(params_inversion['vpvs_min'],params_inversion['vpvs_max'],ndatav)
    vels_xi=np.linspace(params_inversion['xi_min'],params_inversion['xi_max'],ndatav)

    outputs['nostack']['depths']=depths
    outputs['nostack']['vels_vsv']=vels_vsv
    outputs['nostack']['vels_vp']=vels_vp
    outputs['nostack']['vels_xi']=vels_xi
    outputs['nostack']['ndatad']=ndatad
    outputs['nostack']['ndatav']=ndatav

    vsv_model=model['vsv']
    xi_model=model['xi']
    vp_model=model['vp']
    depth_model=model['depth']
    alpha_ref_max=dispersion_ref['alpha_max']
    alpha_max=dispersion_all['alpha_max']

    alpha_ref=min(0.,alpha_ref-alpha_ref_max)
    alpha_ref=np.exp(alpha_ref)

    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    #alpha=np.ma.array(alpha,mask=alpha<-20.)
    alpha=np.exp(alpha)

    ind_vsv=np.digitize(np.interp(depths,depth_model,vsv_model),bins=vels_vsv,right=True)
    vsv_ref[np.arange(ndatad),ind_vsv]+=np.reshape(np.tile(alpha_ref,(ndatad,1)),(ndatad))
    vsv_wide[np.arange(ndatad),ind_vsv]+=np.reshape(np.tile(1.,(ndatad,1)),(ndatad))

    vsv_all[np.arange(ndatad),ind_vsv,:]+=np.tile(np.ma.exp(alpha),(ndatad,1))

    ind_xi=np.digitize(np.interp(depths,depth_model,xi_model),bins=vels_xi,right=True)
    xi_ref[np.arange(ndatad),ind_xi]+=np.reshape(np.tile(alpha_ref,(ndatad,1)),(ndatad))
    xi_wide[np.arange(ndatad),ind_xi]+=np.reshape(np.tile(1.,(ndatad,1)),(ndatad))

    xi_all[np.arange(ndatad),ind_xi,:]+=np.tile(np.ma.exp(alpha),(ndatad,1))

    ind_vp=np.digitize(np.interp(depths,depth_model,vp_model),bins=vels_vp,right=True)


    vp_ref[np.arange(ndatad),ind_vp]+=np.reshape(np.tile(alpha_ref,(ndatad,1)),(ndatad))
    vp_wide[np.arange(ndatad),ind_vp]+=np.reshape(np.tile(1.,(ndatad,1)),(ndatad))

    vp_all[np.arange(ndatad),ind_vp,:]+=np.tile(np.ma.exp(alpha),(ndatad,1))

    outputs['stack']['vsv_ref']=vsv_ref
    outputs['stack']['xi_ref']=xi_ref
    outputs['stack']['vp_ref']=vp_ref
    outputs['stack']['vsv_all']=vsv_all
    outputs['stack']['xi_all']=xi_all
    outputs['stack']['vp_all']=vp_all
    outputs['stack']['vsv_wide']=vsv_wide
    outputs['stack']['xi_wide']=xi_wide
    outputs['stack']['vp_wide']=vp_wide
    #outputs['stack']['alpha_sum']=alpha
    #outputs['stack']['alpha_ref_sum']=alpha_ref
    #outputs['stack']['models_sum_vsv']=np.multiply(np.tile(alpha,(ndatad,1)),np.transpose(np.tile(np.interp(depths,depth_model,vsv_model),(numdis,1))))
    #outputs['stack']['models_sum_vsv_ref']=alpha_ref*np.interp(depths,depth_model,vsv_model)
    #outputs['stack']['models_sum_vp']=np.multiply(np.tile(alpha,(ndatad,1)),np.transpose(np.tile(np.interp(depths,depth_model,vp_model),(numdis,1))))
    #outputs['stack']['models_sum_vp_ref']=alpha_ref*np.interp(depths,depth_model,vp_model)
    #outputs['stack']['models_sum_xi']=np.multiply(np.tile(alpha,(ndatad,1)),np.transpose(np.tile(np.interp(depths,depth_model,xi_model),(numdis,1))))
    #outputs['stack']['models_sum_xi_ref']=alpha_ref*np.interp(depths,depth_model,xi_model)



    return outputs

def get_average(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=True):

    ndatad=200
    outputs={}
    outputs['stack']={}
    outputs['nostack']={}
    # vsv_ref=np.zeros((ndatad))
    # xi_ref=np.zeros_like(vsv_ref)
    # vp_ref=np.zeros_like(vsv_ref)
    # probani_ref=np.zeros_like(vsv_ref)

    # vsv_ref_sq=np.zeros((ndatad))
    # xi_ref_sq=np.zeros_like(vsv_ref)
    # vp_ref_sq=np.zeros_like(vsv_ref)

    # vsv_wide=np.zeros((ndatad))
    # xi_wide=np.zeros_like(vsv_ref)
    # vp_wide=np.zeros_like(vsv_ref)
    # probani_wide=np.zeros_like(vsv_ref)

    # vsv_wide_sq=np.zeros((ndatad))
    # xi_wide_sq=np.zeros_like(vsv_ref)
    # vp_wide_sq=np.zeros_like(vsv_ref)

    # numdis=dispersion_all['numdis']
    # vsv_all=np.zeros((ndatad,numdis))
    # xi_all=np.zeros_like(vsv_all)
    # vp_all=np.zeros_like(vsv_all)
    # probani_all=np.zeros_like(vsv_all)

    # vsv_all_sq=np.zeros((ndatad,numdis))
    # xi_all_sq=np.zeros_like(vsv_all)
    # vp_all_sq=np.zeros_like(vsv_all)



    depths=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)

    outputs['nostack']['depths']=depths
    outputs['nostack']['ndata']=ndatad

    # vsv_model=model['vsv']
    xi_model=model['xi']
    # vp_model=model['vp']
    depth_model=model['depth']

    alpha_ref_max=dispersion_ref['alpha_max']
    alpha_max=dispersion_all['alpha_max']

    alpha_ref=min(0.,alpha_ref-alpha_ref_max)
    alpha_ref=np.exp(alpha_ref)

    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    #alpha=np.ma.array(alpha,mask=alpha<-20.)
    alpha=np.exp(alpha)

    # vsv=np.interp(depths,depth_model,vsv_model)
    # vsv_model_ref=np.interp(depths,model_ref['depth'],model_ref['vsv']) # get reference model

    # vsv_ref=alpha_ref*vsv # reference vsv
    # vsv_ref_sq=alpha_ref*np.square(vsv-vsv_model_ref) # squared reference vsv

    # vsv_wide=vsv # widened vsv
    # vsv_wide_sq=np.square(vsv-vsv_model_ref) # squared widened vsv

    # vsv_alt,alpha_alt=np.meshgrid(vsv,alpha)
    # vsv_alt2,alpha_alt2=np.meshgrid(np.square(vsv-vsv_model_ref),alpha)
    # vsv_all=np.multiply(vsv_alt,alpha_alt) # cluster vsv
    # vsv_all_sq=np.multiply(vsv_alt2,alpha_alt2) # squared cluster vsv

    xi=np.interp(depths,depth_model,xi_model)

    # xi_ref=alpha_ref*xi
    # xi_ref_sq=alpha_ref*(np.square(xi)-np.ones_like(xi))

    # xi_wide=xi
    # xi_wide_sq=np.square(xi)-np.ones_like(xi)

    # xi_alt,alpha_alt=np.meshgrid(np.square(xi-np.ones_like(xi)),alpha)
    # xi_alt2,alpha_alt2=np.meshgrid(xi,alpha)
    # xi_all=np.multiply(xi_alt2,alpha_alt2)
    # xi_all_sq=np.multiply(xi_alt,alpha_alt)


    # vp=np.interp(depths,depth_model,vp_model)

    # vp_ref=alpha_ref*vp
    # vp_ref_sq=alpha_ref*np.square(vp)

    # vp_wide=vp
    # vp_wide_sq=np.square(vp)

    # vp_alt,alpha_alt=np.meshgrid(np.square(vp),alpha)
    # vp_all=np.multiply(vp_alt,alpha_alt)
    # vp_all_sq=np.multiply(vp_alt,alpha_alt)

    probani_wide=(xi!=1.)

    probani_ref=alpha_ref*(xi!=1.)

    probani_alt,alpha_alt=np.meshgrid((xi!=1.),alpha)
    probani_all=np.multiply(probani_alt,alpha_alt)

    outputs['meansum']=[]
    outputs['meansum_weights']=[]
    # outputs['meansum_sq']=[]
    # outputs['meansum_sq_weights']=[]
    # outputs['meansum_sq_shift']=[]
    outputs['meansum_ref']=[]
    outputs['meansum_ref_weights']=[]
    # outputs['meansum_sq_ref']=[]
    # outputs['meansum_sq_ref_weights']=[]
    # outputs['meansum_sq_ref_shift']=[]
    outputs['meansum_wide']=[]
    outputs['meansum_wide_weights']=[]

    outputs['stack']['alpha_sum']=alpha
    outputs['stack']['alpha_ref_sum']=alpha_ref
    outputs['stack']['alpha_wide_sum']=1.

    # outputs['stack']['models_mean_vsv']=vsv_all
    # outputs['meansum'].append('models_mean_vsv')
    # outputs['meansum_weights'].append('alpha_sum')
    # outputs['stack']['models_mean_vsv_ref']=vsv_ref
    # outputs['meansum_ref'].append('models_mean_vsv_ref')
    # outputs['meansum_ref_weights'].append('alpha_ref_sum')
    # outputs['stack']['models_mean_vsv_wide']=vsv_wide
    # outputs['stack']['models_mean_vsv_sq']=vsv_all_sq
    # outputs['stack']['models_mean_vsv_sq_ref']=vsv_ref_sq
    # outputs['stack']['models_mean_vsv_sq_wide']=vsv_wide_sq
    # outputs['nostack']['models_mean_vsv_mean_sq']=vsv_model_ref
    # outputs['stack']['models_mean_vp']=vp_all
    # outputs['stack']['models_mean_vp_ref']=vp_ref
    # outputs['stack']['models_mean_vp_wide']=vp_wide
    # outputs['stack']['models_mean_vp_sq']=vp_all_sq
    # outputs['stack']['models_mean_vp_sq_ref']=vp_ref_sq
    # outputs['stack']['models_mean_vp_sq_wide']=vp_wide_sq
    # outputs['nostack']['models_mean_vp_mean_sq']=np.zeros_like(vp_ref)
    # outputs['stack']['models_mean_xi']=xi_all
    # outputs['stack']['models_mean_xi_ref']=xi_ref
    # outputs['stack']['models_mean_xi_wide']=xi_wide
    # outputs['stack']['models_mean_xi_sq']=xi_all_sq
    # outputs['stack']['models_mean_xi_sq_ref']=xi_ref_sq
    # outputs['stack']['models_mean_xi_sq_wide']=xi_wide_sq
    # outputs['nostack']['models_mean_xi_mean_sq']=np.ones_like(xi_ref)
    outputs['stack']['probani']=probani_all
    outputs['meansum'].append('probani')
    outputs['meansum_weights'].append('alpha_sum')
    outputs['stack']['probani_ref']=probani_ref
    outputs['meansum_ref'].append('probani_ref')
    outputs['meansum_ref_weights'].append('alpha_ref_sum')
    outputs['stack']['probani_wide']=probani_wide
    outputs['meansum_wide'].append('probani_wide')
    outputs['meansum_wide_weights'].append('alpha_wide_sum')




    return outputs

def get_histograms(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=True):
    '''
    example of a function to be applied to the data
    creates a 2D histogram of vsv for the reference and each member of the cluster
    The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

    Parameters
    ----------
    model : dict

    model_ref : dict

    dispersion_one : dict

    params_dispersion : dict

    params_inversion : dict

    dispersion_ref : dict

    dispersion_all : dict

    alpha_ref : float

    alpha : numpy array, length numdis

    first : bool, optional
        whether this model is the first of its file. May provide a minor speed boost. The default is True.

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    ndata=200
    outputs={}
    outputs['stack']={}
    outputs['nostack']={}

    #params_inversion['milay']=5
    #params_inversion['malay']=40

    range_nlay=np.arange(params_inversion['milay'],params_inversion['malay']+1)
    range_sigmaR=np.linspace(params_inversion['Ad_R_min'],
                     params_inversion['Ad_R_max'],
                     ndata)
    range_sigmaL=np.linspace(params_inversion['Ad_L_min'],
                     params_inversion['Ad_L_max'],
                     ndata)
    range_alpha=np.linspace(-100,
                     10,
                     ndata)

    numdis=dispersion_all['numdis']
    nlay=np.zeros((params_inversion['malay']-params_inversion['milay']+1,numdis))
    nlay_ref=np.zeros((params_inversion['malay']-params_inversion['milay']+1))
    nlay_wide=np.zeros_like(nlay_ref)
    sigmaR=np.zeros((ndata,numdis))
    sigmaR_ref=np.zeros((ndata))
    sigmaR_wide=np.zeros_like(sigmaR_ref)
    sigmaL=np.zeros_like(sigmaR)
    sigmaL_ref=np.zeros_like(sigmaR_ref)
    sigmaL_wide=np.zeros_like(sigmaR_ref)
    alphahist=np.zeros_like(sigmaR)
    alphahist_ref=np.zeros_like(sigmaR_ref)

    outputs['nostack']['range_nlay']=range_nlay
    outputs['nostack']['range_sigmaR']=range_sigmaR
    outputs['nostack']['range_sigmaL']=range_sigmaL
    outputs['nostack']['range_alpha']=range_alpha

    nlay_model=model['npt']

    alpha_ref_max=dispersion_ref['alpha_max']
    alpha_max=dispersion_all['alpha_max']
    alpha_ref=min(0.,alpha_ref-alpha_ref_max)
    alpha_ref_old=alpha_ref
    alpha_ref=np.exp(alpha_ref)

    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    alpha_old=alpha
    #alpha=np.ma.array(alpha,mask=alpha<-20.)

    alpha=np.exp(alpha)

    ind_nlay=nlay_model-params_inversion['milay']
    nlay[ind_nlay,:]=alpha
    nlay_ref[ind_nlay]=alpha_ref
    nlay_wide[ind_nlay]+=1

    ind_sigmaR=np.digitize(model['Ad_R'],bins=range_sigmaR)
    sigmaR[ind_sigmaR,:]=alpha
    sigmaR_ref[ind_sigmaR]=alpha_ref
    sigmaR_wide[ind_sigmaR]+=1

    ind_sigmaL=np.digitize(model['Ad_L'],bins=range_sigmaL)
    sigmaL[ind_sigmaL,:]=alpha
    sigmaL_ref[ind_sigmaL]=alpha_ref
    sigmaL_wide[ind_sigmaL]+=1

    ind_alpha=np.digitize(alpha_old,bins=range_alpha)
    alphahist[ind_alpha,np.arange(numdis)]=1
    ind_alpha=np.digitize(alpha_ref_old,bins=range_alpha)
    alphahist_ref[ind_alpha]=1

    outputs['stack']['nlay_hist']=nlay
    outputs['stack']['nlay_ref_hist']=nlay_ref
    outputs['stack']['nlay_wide_hist']=nlay_wide
    outputs['stack']['sigmaR_hist']=sigmaR
    outputs['stack']['sigmaR_ref_hist']=sigmaR_ref
    outputs['stack']['sigmaR_wide_hist']=sigmaR_wide
    outputs['stack']['sigmaL_hist']=sigmaL
    outputs['stack']['sigmaL_ref_hist']=sigmaL_ref
    outputs['stack']['sigmaL_wide_hist']=sigmaL_wide
    outputs['stack']['alpha_hist']=alphahist
    outputs['stack']['alpha_ref_hist']=alphahist_ref


    return outputs

def get_dispersion_mean(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=True):

    outputs={}
    outputs['stack']={}
    outputs['nostack']={}



    dispersion_R_model=dispersion_one['R']['dispersion']
    dispersion_L_model=dispersion_one['L']['dispersion']

    outputs['nostack']['ndata']=len(dispersion_R_model)

    alpha_ref_max=dispersion_ref['alpha_max']
    alpha_max=dispersion_all['alpha_max']
    alpha_ref=min(0.,alpha_ref-alpha_ref_max)
    alpha_ref=np.exp(alpha_ref)
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    #alpha=np.ma.array(alpha,mask=alpha<-20.)
    alpha=np.exp(alpha)
    dispersion_R_true=dispersion_ref['R']['dispersion']
    dispersion_L_true=dispersion_ref['L']['dispersion']

    #numdis=dispersion_all['numdis']

    dis_alt,alpha_alt=np.meshgrid(dispersion_R_model,alpha)
    dispersions_R=np.multiply(dis_alt,alpha_alt)
    dis_alt,alpha_alt=np.meshgrid(dispersion_R_model-dispersion_R_true,alpha)
    dispersions_R_shift=np.multiply(dis_alt,alpha_alt)
    dis_alt,alpha_alt=np.meshgrid(np.square(dispersion_R_model-dispersion_R_true),alpha)
    dispersions_R_sq=np.multiply(dis_alt,alpha_alt)

    dis_alt,alpha_alt=np.meshgrid(dispersion_L_model,alpha)
    dispersions_L=np.multiply(dis_alt,alpha_alt)
    dis_alt,alpha_alt=np.meshgrid(dispersion_L_model-dispersion_L_true,alpha)
    dispersions_L_shift=np.multiply(dis_alt,alpha_alt)
    dis_alt,alpha_alt=np.meshgrid(np.square(dispersion_L_model-dispersion_L_true),alpha)
    dispersions_L_sq=np.multiply(dis_alt,alpha_alt)

    outputs['meansum']=[]
    outputs['meansum_weights']=[]
    outputs['meansum_sq']=[]
    outputs['meansum_sq_weights']=[]
    outputs['meansum_sq_shift']=[]
    outputs['meansum_ref']=[]
    outputs['meansum_ref_weights']=[]
    outputs['meansum_sq_ref']=[]
    outputs['meansum_sq_ref_weights']=[]
    outputs['meansum_sq_ref_shift']=[]
    outputs['meansum_wide']=[]
    outputs['meansum_wide_weights']=[]

    outputs['stack']['dispersion_R']=dispersions_R
    outputs['stack']['dispersion_R_sq']=dispersions_R_sq
    outputs['stack']['dispersion_R_shift']=dispersions_R_shift
    outputs['meansum'].append('dispersion_R')
    outputs['meansum_weights'].append('alpha_sum')
    outputs['meansum_sq'].append('dispersion_R_sq')
    outputs['meansum_sq_weights'].append('alpha_sum')
    outputs['meansum_sq_shift'].append('dispersion_R_shift')

    outputs['stack']['dispersion_R_ref']=alpha_ref*dispersion_R_model
    outputs['stack']['dispersion_R_sq_ref']=alpha_ref*np.square(dispersion_R_model-dispersion_R_true)
    outputs['stack']['dispersion_R_shift_ref']=alpha_ref*(dispersion_R_model-dispersion_R_true)
    outputs['meansum_ref'].append('dispersion_R_ref')
    outputs['meansum_ref_weights'].append('alpha_ref_sum')
    outputs['meansum_sq_ref'].append('dispersion_R_sq_ref')
    outputs['meansum_sq_ref_weights'].append('alpha_ref_sum')
    outputs['meansum_sq_ref_shift'].append('dispersion_R_shift_ref')

    outputs['stack']['dispersion_R_wide']=dispersion_R_model
    outputs['stack']['dispersion_R_sq_wide']=np.square(dispersion_R_model-dispersion_R_true)
    outputs['stack']['dispersion_R_shift_wide']=dispersion_R_model-dispersion_R_true
    outputs['meansum_ref'].append('dispersion_R_wide')
    outputs['meansum_ref_weights'].append('alpha_wide_sum')
    outputs['meansum_sq_ref'].append('dispersion_R_sq_wide')
    outputs['meansum_sq_ref_weights'].append('alpha_wide_sum')
    outputs['meansum_sq_ref_shift'].append('dispersion_R_shift_wide')

    outputs['stack']['dispersion_L']=dispersions_L
    outputs['stack']['dispersion_L_sq']=dispersions_L_sq
    outputs['stack']['dispersion_L_shift']=dispersions_L_shift
    outputs['meansum'].append('dispersion_L')
    outputs['meansum_weights'].append('alpha_sum')
    outputs['meansum_sq'].append('dispersion_L_sq')
    outputs['meansum_sq_weights'].append('alpha_sum')
    outputs['meansum_sq_shift'].append('dispersion_L_shift')

    outputs['stack']['dispersion_L_ref']=alpha_ref*dispersion_L_model
    outputs['stack']['dispersion_L_sq_ref']=alpha_ref*np.square(dispersion_L_model-dispersion_L_true)
    outputs['stack']['dispersion_L_shift_ref']=alpha_ref*(dispersion_L_model-dispersion_L_true)
    outputs['meansum_ref'].append('dispersion_L_ref')
    outputs['meansum_ref_weights'].append('alpha_ref_sum')
    outputs['meansum_sq_ref'].append('dispersion_L_sq_ref')
    outputs['meansum_sq_ref_weights'].append('alpha_ref_sum')
    outputs['meansum_sq_ref_shift'].append('dispersion_L_shift_ref')

    outputs['stack']['dispersion_L_wide']=dispersion_L_model
    outputs['stack']['dispersion_L_sq_wide']=np.square(dispersion_R_model-dispersion_R_true)
    outputs['stack']['dispersion_L_shift_wide']=dispersion_L_model-dispersion_L_true
    outputs['meansum_ref'].append('dispersion_L_wide')
    outputs['meansum_ref_weights'].append('alpha_wide_sum')
    outputs['meansum_sq_ref'].append('dispersion_L_sq_wide')
    outputs['meansum_sq_ref_weights'].append('alpha_wide_sum')
    outputs['meansum_sq_ref_shift'].append('dispersion_L_shift_wide')

    outputs['stack']['alpha_sum']=alpha
    outputs['stack']['alpha_ref_sum']=alpha_ref
    outputs['stack']['alpha_wide_sum']=1.

    return outputs

def get_tradeoff(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,first=True):

    outputs={}
    outputs['stack']={}
    outputs['nostack']={}

    numdis=dispersion_all['numdis']

    alpha_ref_max=dispersion_ref['alpha_max']
    alpha_max=dispersion_all['alpha_max']
    alpha_ref=min(0.,alpha_ref-alpha_ref_max)
    alpha_ref=np.exp(alpha_ref)
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    #alpha=np.ma.array(alpha,mask=alpha<-20.)
    alpha=np.exp(alpha)

    #params_inversion['malay']=40
    malay=params_inversion['malay']

    nlay=model['npt']
    nlay_ani=model['npt_ani']

    tradeoff_ref=np.zeros((malay+1,malay+1))
    tradeoff_wide=np.zeros((malay+1,malay+1))
    tradeoff=np.zeros((malay+1,malay+1,numdis))

    tradeoff_ref[nlay,nlay_ani]=alpha_ref
    tradeoff_wide[nlay,nlay_ani]=1.
    tradeoff[nlay,nlay_ani,:]=alpha

    outputs['nostack']['range']=np.arange(malay)

    outputs['stack']['tradeoff_ref']=tradeoff_ref
    outputs['stack']['tradeoff_wide']=tradeoff_wide
    outputs['stack']['tradeoff']=tradeoff

    return outputs


############################################################################
#        Example
############################################################################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


directory='OUT_CLUSTER1_ALL/OUT'

functions=[create_posterior,get_average,get_histograms,get_dispersion_mean] #create_posterior,get_average,get_histograms,,get_tradeoff

prepare=True
maxpercent=1.
print('start')
params_dispersion, dispersion_ref, dispersion_all=get_dispersion(directory)
print('get dispersion')
model_ref=get_model_ref()
print('get model ref')

if rank==0:
    filename='dispersion.h5'
    f = h5py.File(directory+'/Processing/'+filename,'w')
    grp=f.create_group('cluster_params')
    grp.create_dataset('numdis',data=dispersion_all['numdis'])
    grp.create_dataset('lat',data=dispersion_all['lat'])
    grp.create_dataset('lon',data=dispersion_all['lon'])
    grp.create_dataset('cluster',data=dispersion_all['cluster'])
    grp=f.create_group('dispersion_params')
    grp2=grp.create_group('L')
    grp2.create_dataset('periods',data=params_dispersion['L']['periods'])
    grp2.create_dataset('modes',data=params_dispersion['L']['modes'])
    grp2=grp.create_group('R')
    grp2.create_dataset('periods',data=params_dispersion['R']['periods'])
    grp2.create_dataset('modes',data=params_dispersion['R']['modes'])
    grp=f.create_group('reference')
    grp2=grp.create_group('L')
    grp2.create_dataset('dispersion',data=dispersion_ref['L']['dispersion'])
    grp2.create_dataset('error',data=dispersion_ref['L']['error'])
    grp2=grp.create_group('R')
    grp2.create_dataset('dispersion',data=dispersion_ref['R']['dispersion'])
    grp2.create_dataset('error',data=dispersion_ref['R']['error'])
    grp=f.create_group('cluster')
    grp2=grp.create_group('L')
    grp2.create_dataset('dispersion',data=dispersion_all['L']['dispersion'])
    grp2.create_dataset('error',data=dispersion_all['L']['error'])
    grp2=grp.create_group('R')
    grp2.create_dataset('dispersion',data=dispersion_all['R']['dispersion'])
    grp2.create_dataset('error',data=dispersion_all['R']['error'])
    f.close()

# #print(np.amax(dispersion_all['R']['dispersion']))

# #sys.exit()
for widening in [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]:
#widening=1.0
    for maxpercent in [0.,0.1,1.,0.5]:

        output={}
        params_inversion={}
        params_inversion=get_metadata(directory,prepare=prepare,widening=widening)
        print('get metadata')
        alphafile=directory+'/Processing/alphamax_'+str(prepare)+'_'+str(widening)+'_'+str(maxpercent)+'.txt'
        if os.path.isfile(alphafile):
            file=open(alphafile,'r')
            dispersion_ref['alpha_max']=float(file.readline())
            dispersion_all['alpha_max']=np.array(file.readline().split()).astype('float')
            file.close()
        else:
            get_alpha_max(comm,directory, params_inversion, params_dispersion, dispersion_ref, dispersion_all, prepare=prepare, widening=widening, maxpercent=maxpercent)
            file=open(alphafile,'w')
            file.write(str(dispersion_ref['alpha_max'])+'\n')
            file.write(' '.join(map(str,dispersion_all['alpha_max']))+'\n')
            file.close()

        print('get alphamax')
        output=apply_stuff(comm,directory,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,prepare=prepare,widening=widening) #
        print('apply functions')

        if rank==0:

            ############################################################
            # print to files
            ############################################################
            for function in output:
                filename='processing_'+str(widening)+'_'+str(prepare)+'_'+str(maxpercent)+'_'+function+'_outputs.h5'
                f = h5py.File(directory+'/Processing/'+filename,'w')
                f.create_dataset('widening',data=params_inversion['widening'])
                f.create_dataset('prepare',data=prepare)
                f.create_dataset('burn-in',data=params_inversion['burn-in'])
                f.create_dataset('nsample',data=params_inversion['nsample'])


                grp=f.create_group(function)
                grp_stack=grp.create_group('stack')
                grp_nostack=grp.create_group('nostack')
                for key in output[function]['stack']:
                    grp_stack.create_dataset(key,data=output[function]['stack'][key])
                    output[function]['stack'][key]=[]
                for key in output[function]['nostack']:
                    grp_nostack.create_dataset(key,data=output[function]['nostack'][key])
                    output[function]['nostack'][key]=[]
                f.close()
                output[function]={}

        #del params_inversion,output,data,weights

for maxpercent in [0.,0.1,1.,0.5]:
    output={}
    prepare=False
    params_inversion=get_metadata(directory,prepare=prepare)
    print('get metadata')
    alphafile=directory+'/Processing/alphamax_'+str(prepare)+'_'+str(maxpercent)+'.txt'
    if os.path.isfile(alphafile):
        file=open(alphafile,'r')
        dispersion_ref['alpha_max']=float(file.readline())
        dispersion_all['alpha_max']=np.array(file.readline().split()).astype('float')
        file.close()
    else:
        get_alpha_max(comm,directory,params_inversion,params_dispersion,dispersion_ref,dispersion_all,prepare=prepare,maxpercent=maxpercent)
        if rank==0:
            file=open(alphafile,'w')
            file.write(str(dispersion_ref['alpha_max'])+'\n')
            file.write(' '.join(map(str,dispersion_all['alpha_max']))+'\n')
            file.close()

    continue
    widening=params_inversion['widening']
    print('get alphamax')
    output=apply_stuff(comm,directory,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,prepare=prepare,widening=widening)
    print('apply functions')

    ##########################################################
    # Print to file
    ##########################################################
    if rank==0:
        for function in output:
            filename='processing_'+'_'+str(maxpercent)+'_'+function+'_outputs.h5'
            f = h5py.File(directory+'/Processing/'+filename,'w')
            f.create_dataset('widening',data=params_inversion['widening'])
            f.create_dataset('prepare',data=prepare)
            f.create_dataset('burn-in',data=params_inversion['burn-in'])
            f.create_dataset('nsample',data=params_inversion['nsample'])


            grp=f.create_group(function)
            grp_stack=grp.create_group('stack')
            grp_nostack=grp.create_group('nostack')
            for key in output[function]['stack']:
                grp_stack.create_dataset(key,data=output[function]['stack'][key])
            for key in output[function]['nostack']:
                grp_nostack.create_dataset(key,data=output[function]['nostack'][key])
            f.close()

    #del params_inversion,output

