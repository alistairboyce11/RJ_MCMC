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
warnings.filterwarnings("error")

class Probsum1(Exception):
    pass

def get_dispersion(directories,clusters,points=[],filename='dispersion_all.in'):
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
    dispersion_ref={}
    dispersion_ref['R']={}
    dispersion_ref['L']={}
    dispersion_all={}
    dispersion_all['R']={}
    dispersion_all['L']={}
    params_dispersion={}
    params_dispersion['R']={}
    params_dispersion['L']={}
    params_dispersion['clusters']=np.array(clusters)

    for ind_cluster in range(len(clusters)):

        # different references for different clusters

        file=open(directories[ind_cluster]+'/'+filename,'r')

        numdis_tmp=int(file.readline())
        if type(points)==str and points=='all':
            numdis=numdis_tmp
        elif len(points)==0:
            numdis=0
        else:
            numdis=len(points)
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

                for l in range(numdis_tmp):
                    line=file.readline()
                    #print(line)
                    data=line.split()
                    #print(data)
                    if type(points)==str and points=='all':
                        dispersion_R[j,l]=float(data[0])
                        error_R[j,l]=float(data[1])
                    else:
                        if l not in points:
                            pass
                        else:
                            l_ind=list(points).index(l)
                            dispersion_R[j,l_ind]=float(data[0])
                            error_R[j,l_ind]=float(data[1])
                j+=1

        # fill reference dispersion curves: one line per cluster


        if ind_cluster==0:
            dispersion_ref['R']['dispersion']=np.atleast_2d(dispersion_R_ref)
            dispersion_ref['R']['error']=np.atleast_2d(error_R_ref)
            dispersion_ref['R']['error_sum']=[np.sum(np.log(error_R_ref))]
        else:
            ar1,ar2=np.atleast_2d(dispersion_ref['R']['dispersion'],dispersion_R_ref)
            dispersion_ref['R']['dispersion']=np.append(ar1,ar2,axis=0)
            ar1,ar2=np.atleast_2d(dispersion_ref['R']['error'],error_R_ref)
            dispersion_ref['R']['error']=np.append(ar1,ar2,axis=0)
            dispersion_ref['R']['error_sum'].append(np.sum(np.log(error_R_ref)))

        # individual data, identical every time
        dispersion_all['R']['dispersion']=dispersion_R
        dispersion_all['R']['error']=error_R
        dispersion_all['R']['error_sum']=np.sum(np.log(error_R),axis=0)
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

                for l in range(numdis_tmp):
                    data=file.readline().split()
                    if type(points)==str and points=='all':
                        dispersion_L[j,l]=float(data[0])
                        error_L[j,l]=float(data[1])
                    else:
                        if l not in points:
                            continue
                        else:
                            l_ind=list(points).index(l)
                            dispersion_L[j,l_ind]=float(data[0])
                            error_L[j,l_ind]=float(data[1])

                j+=1

        if ind_cluster==0:
            dispersion_ref['L']['dispersion']=np.atleast_2d(dispersion_L_ref)
            dispersion_ref['L']['error']=np.atleast_2d(error_L_ref)
            dispersion_ref['L']['error_sum']=[np.sum(np.log(error_L_ref))]
        else:
            ar1,ar2=np.atleast_2d(dispersion_ref['L']['dispersion'],dispersion_L_ref)
            dispersion_ref['L']['dispersion']=np.append(ar1,ar2,axis=0)
            ar1,ar2=np.atleast_2d(dispersion_ref['L']['error'],error_L_ref)
            dispersion_ref['L']['error']=np.append(ar1,ar2,axis=0)
            dispersion_ref['L']['error_sum'].append(np.sum(np.log(error_L_ref)))
        dispersion_all['L']['dispersion']=dispersion_L
        dispersion_all['L']['error']=error_L
        dispersion_all['L']['error_sum']=np.sum(np.log(error_L),axis=0)
        params_dispersion['L']['periods']=period_L
        params_dispersion['L']['modes']=mode_L

        file.close()
    dispersion_ref['R']['error_sum']=np.array(dispersion_ref['R']['error_sum'])
    dispersion_ref['L']['error_sum']=np.array(dispersion_ref['L']['error_sum'])
    return params_dispersion, dispersion_ref, dispersion_all

def get_metadata(directories,widenings):
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

    files_all=[]
    for i in range(len(directories)):
        directory=directories[i]
        for j in range(len(widenings)):

            files_all.extend(glob.glob(directory+'/All_models_processed_prepare_*%4.2f.out'%widenings[j]))

    for file in files_all:

        f=open(file,'r')

        #print(file)

        params_inversion=read_header(f)
        f.close()

        if not params_inversion:
            continue
        else:
            #print('found data! ')
            return params_inversion
    return


def process_file_alphamax(file,cluster,params_dispersion,dispersion_ref,dispersion_all,maxpercent,num_to_store):
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
    f=open(file,'r')
    #print(file)

    numdis=dispersion_all['numdis']
    alpha_ref_max=np.ones((num_to_store))*float('-inf')
    alpha_max=np.ones((num_to_store,numdis))*float('-inf')

    # fline=f.readline()
    # if not fline:
    #     return
    # widening=float(f.readline().split()[2])
    # f.readline()
    # f.readline()
    # f.readline()
    # f.readline()
    # f.readline()
    # f.readline()
    # f.readline()

    params_inversion=read_header(f)
    if not params_inversion:
        return

    widening=params_inversion['widening']

    num_model=0
    valid_model=True
    while valid_model:
        valid_model,dispersion_one,model=read_model(f,widening,cluster)
        if not valid_model:
            continue
        elif len(dispersion_one.keys())==0:
            continue

        # calculate alpha
        alpha_ref,alpha=get_alpha(dispersion_one,cluster,params_dispersion,dispersion_ref,dispersion_all)

        if maxpercent==0.:
            alpha_ref_max=max(alpha_ref,alpha_ref_max)
            alpha_max=np.maximum.reduce([alpha_max,np.atleast_2d(alpha)])

        else:

            if np.any(np.amin(alpha_ref_max,axis=0)<alpha):
                min_indices=np.argmin(alpha_ref_max)
                alpha_ref_max[min_indices]=np.maximum(alpha_ref,np.amin(alpha_ref_max,axis=0))

            if np.any(np.amin(alpha_max,axis=0)<alpha):
                min_indices=np.argmin(alpha_max,axis=0)
                alpha_max[min_indices,np.arange(numdis)]=np.maximum(alpha,np.amin(alpha_max,axis=0))

        num_model+=1
        #print(num_model)

    f.close()

    return alpha_max,alpha_ref_max,num_model

def get_alpha(dispersion_one,cluster,params_dispersion,dispersion_ref,dispersion_all):
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

    #i_cluster=np.where(params_dispersion['clusters']==cluster)[0][0]
    widening=dispersion_one['widening']
    ndatad_R=params_dispersion['R']['ndatad']
    ndatad_L=params_dispersion['L']['ndatad']
    dispersion_R_ref=dispersion_ref['R']['dispersion']
    #dispersion_L_ref=dispersion_ref['L']['dispersion']
    dispersion_R=dispersion_all['R']['dispersion']
    dispersion_L=dispersion_all['L']['dispersion']
    dispersion_R_one=dispersion_one['R']['dispersion']
    dispersion_L_one=dispersion_one['L']['dispersion']
    Ad_R=dispersion_one['R']['Ad']
    Ad_L=dispersion_one['L']['Ad']
    #error_R_ref=dispersion_ref['R']['error']
    #error_L_ref=dispersion_ref['L']['error']
    error_R=dispersion_all['R']['error']
    error_L=dispersion_all['L']['error']
    #error_R_ref_sum=dispersion_ref['R']['error_sum']
    #error_L_ref_sum=dispersion_ref['L']['error_sum']
    #error_R_sum=dispersion_all['R']['error_sum']
    #error_L_sum=dispersion_all['L']['error_sum']
    like_w=dispersion_one['like_w']
    numdis=dispersion_all['numdis']

    if np.shape(dispersion_R_ref)==(0,):
        print('here')
    if np.shape(dispersion_R_one)==(0,):
        print('here2')
        print(dispersion_R_one)

    like_alt=np.sum(np.divide(np.square(dispersion_R-np.transpose(np.tile(dispersion_R_one,(numdis,1)))),2*Ad_R**2*np.square(error_R)),axis=0)
    like_alt-=np.sum(np.divide(np.square(dispersion_L-np.transpose(np.tile(dispersion_L_one,(numdis,1)))),2*Ad_L**2*np.square(error_L)),axis=0)
    # static shift
    #like_alt-=error_R_sum
    #like_alt-=error_L_sum
    like_alt-=ndatad_R*np.log(Ad_R)
    like_alt-=ndatad_L*np.log(Ad_L)

    # both contained in like_w
    #like_one=np.sum(np.divide(np.square(dispersion_R-dispersion_R_ref),2*Ad_R**2*np.square(error_R_ref)),axis=0)
    #like_one-=np.sum(np.divide(np.square(dispersion_L-dispersion_L_ref),2*Ad_L**2*np.square(error_L_ref)),axis=0)
    # static
    #like_one-=error_R_ref_sum
    #like_one-=error_L_ref_sum
    # contained in like_w after post-processing
    #like_one-=ndatad_R*Ad_R
    #like_one-=ndatad_L*Ad_L
    #like_one/=widening

    alpha=like_alt-like_w
    alpha_ref=like_w*(widening-1)

    return alpha_ref,alpha

def process_file_alphamax_all(file,clusters,widenings,cluster,widening,params_dispersion,dispersion_ref,dispersion_all,maxpercent,num_to_store):
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
    f=open(file,'r')


    num_model=0

    numdis=dispersion_all['numdis']
    alpha_ref_max=np.ones((num_to_store))*float('-inf')
    alpha_max=np.ones((num_to_store,numdis))*float('-inf')

    params_inversion=read_header(f)
    if not params_inversion:
        return

    widening=params_inversion['widening']

    # loop on all models
    valid_model=True
    cnt=0
    #print(file)
    while valid_model:
        valid_model,dispersion_one,model=read_model(f,widening,cluster)
        if not valid_model:
            continue
        elif len(dispersion_one.keys())==0:
            continue
        elif len(dispersion_one['R']['dispersion'])==0:
            continue
        #elif model['npt_true']<8:
        #    continue
        cnt+=1

        #if num_model==999:
        #    valid_model=False

        # calculate alpha
        alpha_ref,alpha=get_alpha_all(dispersion_one,clusters,widenings,params_dispersion,dispersion_ref,dispersion_all)

        if maxpercent==0.:
            alpha_ref_max=max(alpha_ref,alpha_ref_max)
            #print(alpha_max,np.atleast_2d(alpha))
            alpha_max=np.maximum.reduce([alpha_max,np.atleast_2d(alpha)])

        else:

            if np.any(np.amin(alpha_ref_max,axis=0)<alpha):
                min_indices=np.argmin(alpha_ref_max)
                alpha_ref_max[min_indices]=np.maximum(alpha_ref,np.amin(alpha_ref_max,axis=0))

            if np.any(np.amin(alpha_max,axis=0)<alpha):
                min_indices=np.argmin(alpha_max,axis=0)
                alpha_max[min_indices,np.arange(numdis)]=np.maximum(alpha,np.amin(alpha_max,axis=0))

        num_model+=1

        # if (num_model%10000)==0:
        #     print(num_model)

    f.close()

    return alpha_max,alpha_ref_max,num_model

def get_alpha_all(dispersion_one,clusters,widenings,params_dispersion,dispersion_ref,dispersion_all):
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

    ndatad_R=params_dispersion['R']['ndatad']
    ndatad_L=params_dispersion['L']['ndatad']
    dispersion_R=dispersion_all['R']['dispersion']
    dispersion_L=dispersion_all['L']['dispersion']
    dispersion_R_one=dispersion_one['R']['dispersion']
    dispersion_L_one=dispersion_one['L']['dispersion']
    Ad_R=dispersion_one['R']['Ad']
    Ad_L=dispersion_one['L']['Ad']
    error_R=dispersion_all['R']['error']
    error_L=dispersion_all['L']['error']
    error_R_ref_sum=dispersion_ref['R']['error_sum']
    error_L_ref_sum=dispersion_ref['L']['error_sum']
    #error_R_sum=dispersion_all['R']['error_sum']
    #error_L_sum=dispersion_all['L']['error_sum']
    numdis=dispersion_all['numdis']

    like_alt=-np.sum(np.divide(
        np.square(dispersion_R-np.transpose(np.tile(dispersion_R_one,(numdis,1)))),
        2*Ad_R**2*
        np.square(
            np.multiply(
                error_R,dispersion_R/100.)))
        ,axis=0)
    #print(len(dispersion_L),len(dispersion_L_one),len(error_L))
    like_alt-=np.sum(np.divide(
        np.square(dispersion_L-np.transpose(np.tile(dispersion_L_one,(numdis,1)))),
        2*Ad_L**2*
        np.square(
            np.multiply(
                error_L,dispersion_L/100.)))
        ,axis=0)
    # static shift
    #like_alt-=error_R_sum
    #like_alt-=error_L_sum
    like_alt-=ndatad_R*np.log(Ad_R)
    like_alt-=ndatad_L*np.log(Ad_L)

    #print(numdis)



    # alpha=log(like_alt/sum(cluster,widening like_one) )
    # alpha=-log(sum(cluster,widening like_one/like_alt))

    probsum=0.
    shift=0
    while probsum==0.:
        try:
            numprob=0.
            probsum=0.
            for i in range(len(clusters)):
                for widening in widenings:

                    dispersion_R_ref=dispersion_ref['R']['dispersion'][i,:]
                    dispersion_L_ref=dispersion_ref['L']['dispersion'][i,:]

                    error_R_ref=dispersion_ref['R']['error'][i,:]
                    error_L_ref=dispersion_ref['L']['error'][i,:]

                    if np.shape(dispersion_R_ref)==(0,):
                        print('here')
                    if np.shape(dispersion_R_one)==(0,):
                        print('here2')
                        print(dispersion_R_one)

                    #print('new')
                    like_one=-np.sum(np.divide(np.square(dispersion_R_ref-dispersion_R_one),2*Ad_R**2*np.square(np.multiply(error_R_ref/100.,dispersion_R_ref))),axis=0)
                    #print(like_one)
                    like_one-=np.sum(np.divide(np.square(dispersion_L_ref-dispersion_L_one),2*Ad_L**2*np.square(np.multiply(error_L_ref/100.,dispersion_L_ref))),axis=0)
                    #print(like_one)
                    like_one-=error_R_ref_sum[i] + 150 # static shift to avoid exponential explosion
                    #print(like_one)
                    like_one-=error_L_ref_sum[i] + 150 # static shift to avoid exponential explosion
                    #print(like_one)
                    like_one-=ndatad_R*np.log(Ad_R)
                    #print(like_one)
                    like_one-=ndatad_L*np.log(Ad_L)
                    #print(like_one)
                    like_one/=widening
                    #print(like_one)
                    #like_one+=1000
                    #print(like_one)

                    like_one+=shift

                    probsum+=np.exp(like_one)

                    numprob+=1
            #print(probsum)
            if probsum==0.:
                raise Probsum1
        except RuntimeWarning:
            probsum=0.
            shift-=100
            print('too big probsum')
        except Probsum1:
            probsum=0.
            shift+=1000
            #print('too small probsum')
    #print('done')
    probsum/=numprob
    # if probsum==0.:
    #     alpha=np.ones_like(like_alt)*(-float('inf'))
    # else:
    alpha=like_alt-(np.log(probsum)-shift)

    alpha_ref=0.

    return alpha_ref,alpha

def read_header(f):

    params_inversion={}

    d=f.readline() # contains rank of processor, number of file for the processor and number of models stored at most. Not used or tested currently
    params_inversion['everyall']=int(float(d.split()[1]))
    data=f.readline().split()
    params_inversion['burn-in']=float(data[0])
    params_inversion['widening']=float(data[1])
    params_inversion['thin']=float(data[2])
    data=f.readline().split()
    params_inversion['d_min']=float(data[0])
    params_inversion['d_max']=float(data[1])
    params_inversion['width_vsv']=float(f.readline())
    data=f.readline().split()
    params_inversion['xi_min']=float(data[0])
    params_inversion['xi_max']=float(data[1])
    data=f.readline().split()
    params_inversion['vp_min']=float(data[0])
    params_inversion['vp_max']=float(data[1])
    data=f.readline().split()
    params_inversion['Ad_R_min']=float(data[0])
    params_inversion['Ad_R_max']=float(data[1])
    data=f.readline().split()
    params_inversion['Ad_L_min']=float(data[0])
    params_inversion['Ad_L_max']=float(data[1])
    data=f.readline().split()
    if len(data)==0:
        return
    params_inversion['milay']=int(data[0])
    params_inversion['malay']=int(data[1])

    return params_inversion

def read_model(f,widening,cluster):
    dispersion_one={}
    dispersion_one['R']={}
    dispersion_one['L']={}

    model={}

    #try:

    line=f.readline()
    if not line:
        #print('file end')
        return False,{},{}

    data=line.split()
    # if len(data)<3:
    #     print(line)
    #     print('not much data')
    #     return False,{},{}
    npt_true=int(data[0])
    if npt_true==0:
        #print(line)
        #print('no model')
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        return True,{},{}
    #print(data)
    model['npt_true']=npt_true
    model['npt']=int(data[1])
    model['npt_ani']=int(data[2])
    data=f.readline().split()
    dispersion_one['R']['Ad']=float(data[0])
    dispersion_one['L']['Ad']=float(data[1])
    model['Ad_R']=dispersion_one['R']['Ad']
    model['Ad_L']=dispersion_one['L']['Ad']




    # new output
    d=np.array(f.readline().split()).astype('float')
    vsv=np.array(f.readline().split()).astype('float')
    xi=np.array(f.readline().split()).astype('float')
    vp=np.array(f.readline().split()).astype('float')

    # old output
    # d=np.zeros((npt_true))
    # vsv=np.zeros_like(d)
    # xi=np.zeros_like(d)
    # vp=np.zeros_like(d)
    # for i in range(npt_true):
    #     data=f.readline().split()
    #     d[i]=float(data[0])
    #     vsv[i]=float(data[1])
    #     xi[i]=float(data[2])
    #     vp[i]=float(data[3])
    model['depth']=d
    model['vsv']=vsv
    model['xi']=xi
    model['vp']=vp

    dispersion_one['like_w']=float(f.readline())
    model['like_w']=dispersion_one['like_w']
    dispersion_one['widening']=widening
    dispersion_one['cluster']=cluster

    len_R=int(f.readline())
    dispersion_R_one=np.array(f.readline().split()).astype('float')
    dispersion_one['R']['dispersion']=dispersion_R_one

    len_L=int(f.readline())
    dispersion_L_one=np.array(f.readline().split()).astype('float')

    dispersion_one['L']['dispersion']=dispersion_L_one

    if len(dispersion_L_one)<len_L or len(dispersion_R_one)<len_R or len(d)<npt_true or len(vsv)<npt_true or len(xi)<npt_true or len(vp)<npt_true:
        return True,{},{}

    if len(dispersion_L_one)==0 or len(dispersion_R_one)==0 or len(d)==0 or len(vsv)==0 or len(xi)==0 or len(vp)==0:
        return True,{},{}

    #print(len(dispersion_L_one),len_L)

    return True,dispersion_one,model
    # except:
    #     print('file finished prematurely')
    #     return False,{},{}

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

    file=filename

    model_ref={}

    f=open(file,'r')
    lines=f.readlines()
    f.close()

    data=lines[0].split()
    ntot=int(data[0])
    #nic=int(data[1])
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
        vpv[i]=float(data[2])/1000.
        vsv[i]=float(data[3])/1000.
        vph[i]=float(data[6])/1000.
        vsh[i]=float(data[7])/1000.

    model_ref['npt_true']=npt_true
    model_ref['depth']=d
    model_ref['vpv']=vpv
    model_ref['vph']=vph
    model_ref['vsv']=vsv
    model_ref['vsh']=vsh

    return model_ref

def process_one_file(file,cluster,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref,outputs_all={}):
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

    numdis=dispersion_all['numdis']
    numtot=0
    alpha_sum=np.zeros((numdis))
    exp_sum=np.zeros((numdis))

    # read file
    f=open(file,'r')
    #print(file)

    params_inversion=read_header(f)

    widening=params_inversion['widening']

    valid_model=True
    while valid_model:

        t1=time.time()
        valid_model,dispersion_one,model=read_model(f,widening,cluster)
        if not valid_model:
            continue
        elif len(dispersion_one.keys())==0:
            continue

        t2=time.time()
        time_read+=t2-t1

        t1=time.time()

        alpha_ref,alpha=get_alpha(dispersion_one,cluster,params_dispersion,dispersion_ref,dispersion_all)

        t2=time.time()
        time_alpha+=t2-t1

        t1=time.time()

        # apply functions
        for function in functions:
            function(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha_ref,alpha,outputs=outputs_all[function.__name__])

        numtot+=1
        alpha_sum+=alpha
        exp_sum+=np.exp(alpha)

        t2=time.time()
        time_apply+=t2-t1

        # testing
        # if numtot==10:
        #     f.close()
        #     print('read: ',time_read)
        #     print('get alpha: ',time_alpha)
        #     print('apply: ',time_apply)
        #     return alpha_sum,numtot,exp_sum,exp_ssq

    f.close()

    print('read: ',time_read)
    print('get alpha: ',time_alpha)
    print('apply: ',time_apply)

    return alpha_sum,numtot,exp_sum


def process_one_file_all(file,cluster,clusters_in,widenings_in,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref,outputs_all={}):
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

    numdis=dispersion_all['numdis']
    alpha_sum=np.zeros((numdis+1))
    exp_sum=np.zeros((numdis+1))

    # read file
    f=open(file,'r')
    #print(file)

    params_inversion=read_header(f)
    widening=params_inversion['widening']

    numtot=0
    valid_model=True
    while valid_model:

        t1=time.time()
        valid_model,dispersion_one,model=read_model(f,widening,cluster)
#        if numtot>1000:
#            valid_model=False
        if not valid_model:
            continue
        elif len(dispersion_one.keys())==0:
            continue
       # elif model['npt_true']<8:
       #     continue

        t2=time.time()
        time_read+=t2-t1

        t1=time.time()


        alpha_ref,alpha=get_alpha_all(dispersion_one,clusters_in,widenings_in,params_dispersion,dispersion_ref,dispersion_all)
        alpha-=np.log(params_dispersion['num_models'][cluster][widening])

        #print(np.shape(alpha))

#        print(np.amax(alpha-dispersion_all['alpha_max']))

        t2=time.time()
        time_alpha+=t2-t1
        t1=time.time()

        # apply functions
        for function in functions:
            function(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha_ref,alpha,outputs=outputs_all[function.__name__])

        numtot+=1
        alpha_sum[:numdis]+=np.minimum(np.zeros_like(alpha),alpha-dispersion_all['alpha_max'])

        exp_sum[:numdis]+=np.exp(np.minimum(np.zeros_like(alpha),alpha-dispersion_all['alpha_max']))
        exp_sum[numdis]+=1.

        t2=time.time()
        time_apply+=t2-t1
        #if numtot%10000==0:
            #print(numtot)

        # testing
        # if numtot==10:
        #     f.close()
        #     print('read: ',time_read)
        #     print('get alpha: ',time_alpha)
        #     print('apply: ',time_apply)

        #     return alpha_sum,numtot,exp_sum,exp_ssq

    f.close()

    print('read: ',time_read)
    print('get alpha: ',time_alpha)
    print('apply: ',time_apply)

    return alpha_sum,numtot,exp_sum
