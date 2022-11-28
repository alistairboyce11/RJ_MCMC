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
from mpi4py import MPI
import postprocess_util
import postprocess_functions
import copy


def get_alpha_max(comm,directories,widenings,clusters_in,params_inversion, params_dispersion,dispersion_ref,dispersion_all,maxpercent=0.,widening=1.,alphafile='OUT_KL/Processing/alpha_max.txt'):
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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if os.path.isfile(alphafile):
        file=open(alphafile,'r')
        lines=f.readlines()
        file.close()
        dispersion_ref['alpha_max']={}
        dispersion_all['alpha_max']={}
        for line in lines:
            data=line.split()
            cluster_cur=int(data[0])
            widening=float(data[1])
            alpha_ref_max=float(data[2])
            alpha_max=np.array(line[3:]).astype('float')


            if ['alpha_max'] in dispersion_all:
                if cluster_cur in dispersion_all['alpha_max']:
                    dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max
                else:
                    dispersion_all['alpha_max'][cluster_cur]={}
                    dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max
            else:
                dispersion_all['alpha_max']={}
                dispersion_all['alpha_max'][cluster_cur]={}
                dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max

            if ['alpha_max'] in dispersion_ref:
                if cluster_cur in dispersion_ref['alpha_max']:
                    dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max
                else:
                    dispersion_ref['alpha_max'][cluster_cur]={}
                    dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max
            else:
                dispersion_ref['alpha_max']={}
                dispersion_ref['alpha_max'][cluster_cur]={}
                dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max

    return

    for i_dir in range(len(directories)):
        directory=directories[i_dir]
        cluster_cur=clusters_in[i_dir]
        for widening in widenings:
            # get one alpha max per cluster and per widening
            files_all=[]

            # list files
            files_all.extend(glob.glob(directory+'/All_models_prepare_*%4.2f.out'%widening))

            rank = comm.Get_rank()

            numdis=dispersion_all['numdis']

            # split it over the cores
            size = comm.Get_size()
            files=[]
            clusters=[]
            for i in range(len(files_all)):
                if i%size==rank:
                    files.append(files_all[i])
                    clusters.append(cluster_cur)

            # create array for max alpha, potentially more than one if maxpercent is bigger than zero
            num_to_store=0
            if maxpercent==0.:

                num_to_store=1
            else:
                for file in files_all:
                    num_to_store+=maxpercent/100.*params_inversion[cluster_cur][widening]['everyall']
            num_to_store=int(num_to_store)

            alpha_ref_max=np.ones((num_to_store))*float('-inf')
            alpha_max=np.ones((num_to_store,numdis))*float('-inf')

            print(rank,files)
            print(num_to_store)
            i=0
            num_models=0

            # for this cluster and widening, get the reference dispersion curves
            dispersion_ref_cluster={}
            dispersion_ref_cluster['R']={}
            dispersion_ref_cluster['L']={}

            i_cluster=np.where(params_dispersion['clusters']==cluster_cur)[0][0]

            dispersion_ref_cluster['R']['dispersion']=dispersion_ref['R']['dispersion'][i_cluster,:]
            dispersion_ref_cluster['R']['error']=dispersion_ref['R']['error'][i_cluster,:]
            dispersion_ref_cluster['R']['error_sum']=dispersion_ref['R']['error_sum'][i_cluster]

            dispersion_ref_cluster['L']['dispersion']=dispersion_ref['L']['dispersion'][i_cluster,:]
            dispersion_ref_cluster['L']['error']=dispersion_ref['L']['error'][i_cluster,:]
            dispersion_ref_cluster['L']['error_sum']=dispersion_ref['L']['error_sum'][i_cluster]


            for i in range(len(files)):

                file=files[i]
                cluster=cluster_cur



                # process one file
                alphas_max_prop,alpha_ref_max_prop,num_models_prop=postprocess_util.process_file_alphamax(file,cluster,params_dispersion,dispersion_ref_cluster,dispersion_all,maxpercent,num_to_store)

                # merge results on one core
                alpha_ref_max=np.append(alpha_ref_max,alpha_ref_max_prop)
                alpha_ref_max.sort()
                alpha_ref_max=alpha_ref_max[-num_to_store:]

                alpha_max=np.append(alpha_max,alphas_max_prop,axis=0)
                alpha_max.sort(axis=0)
                alpha_max=alpha_max[-num_to_store:,:]

                num_models+=num_models_prop

            # merge results between cores on core 0
            if rank!=0:

                comm.Send([alpha_max,MPI.FLOAT],dest=0,tag=1)
                comm.Send([alpha_ref_max,MPI.FLOAT],dest=0,tag=2)
                comm.Send([np.array([num_models]),MPI.FLOAT],dest=0,tag=3)
            if rank==0:
                for i in range(size)[1:]:
                    alpha_max_prop=np.zeros_like(alpha_max)
                    comm.Recv([alpha_max_prop,MPI.FLOAT],i,1)
                    for j in range(num_to_store):
                        if np.any(np.amin(alpha_max,axis=0)<alpha_max_prop[j,:]):
                            min_indices=np.argmin(alpha_max,axis=0)
                            alpha_max[min_indices,np.arange(numdis)]=np.maximum(alpha_max_prop[j,:],np.amin(alpha_max,axis=0))

                    alpha_max_ref_prop=np.zeros_like(alpha_ref_max)
                    comm.Recv([alpha_max_ref_prop,MPI.FLOAT],i,2)

                    for j in range(num_to_store):
                        if np.amin(alpha_ref_max)<alpha_max_ref_prop[j]:
                            min_indices=np.argmin(alpha_ref_max)
                            alpha_ref_max[min_indices]=alpha_max_ref_prop[j]

                    num_models_prop=np.zeros_like(num_models)
                    comm.Recv([num_models_prop,MPI.FLOAT],i,3)
                    num_models+=num_models_prop

            # real size of alpha_max
            real_num_to_store=int(num_models*maxpercent/100.)

            alpha_max.sort(axis=0)
            alpha_max=alpha_max[-real_num_to_store,:]

            alpha_ref_max.sort()
            alpha_ref_max=alpha_ref_max[-real_num_to_store]

            alpha_max = comm.bcast(alpha_max, root=0)
            alpha_ref_max = comm.bcast(alpha_ref_max, root=0)

            # put alpha_max in dispersion_all and dispersion_ref for its cluster and widening
            if rank==0:
                if i_dir==0 and widening==widenings[0]:
                    file_alpha=open(alphafile,'w')
                file_alpha.write(str(cluster_cur)+' '+str(widening)+' '+alpha_ref_max+' '+' '.join(np.array(alpha_max).astype('str'))+'\n')
            if ['alpha_max'] in dispersion_all:
                if cluster_cur in dispersion_all['alpha_max']:
                    dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max
                else:
                    dispersion_all['alpha_max'][cluster_cur]={}
                    dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max
            else:
                dispersion_all['alpha_max']={}
                dispersion_all['alpha_max'][cluster_cur]={}
                dispersion_all['alpha_max'][cluster_cur][widening]=alpha_max

            if ['alpha_max'] in dispersion_ref:
                if cluster_cur in dispersion_ref['alpha_max']:
                    dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max
                else:
                    dispersion_ref['alpha_max'][cluster_cur]={}
                    dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max
            else:
                dispersion_ref['alpha_max']={}
                dispersion_ref['alpha_max'][cluster_cur]={}
                dispersion_ref['alpha_max'][cluster_cur][widening]=alpha_ref_max

    if rank==0:
        file_alpha.close()
    return

def apply_stuff(comm,directories,clusters_in,widenings,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,alphafile='OUT_KL/Processing/alpha_max.txt',kl_file='OUT_KL/Processing/fk_dist.txt'):
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

    The function must give the following output:
    output: dict
        contains the output of the function. Must have 2 keys, 'stack' and 'nostack' that each have sub-dicts
        values in output['stack'] will be added for key and each model
        values in output['nostack'] will be kept as in the first model
    '''
    outputs_all_final={}
    for function in functions:
        outputs_all_final[function.__name__]={}
        outputs_all_final[function.__name__]['stack']={}
        outputs_all_final[function.__name__]['nostack']={}
    numdis=dispersion_all['numdis']

    rank = comm.Get_rank()

    size = comm.Get_size()

    if rank==0:
        file_kl=open(kl_file,'w')

    kl_sum=np.zeros((numdis))

    for i_dir in range(len(directories)):
        directory=directories[i_dir]
        cluster_cur=clusters_in[i_dir]

        for i_wide in range(len(widenings)):
            widening=widenings[i_wide]

            files_all=glob.glob(directory+'/All_models_prepare*%4.2f.out'%widening)

            files=[]
            for i in range(len(files_all)):
                if i%size==rank:
                    files.append(files_all[i])

            files.sort()

            file_alpha=open(alphafile,'r')
            line=file_alpha.readline()
            data=line.split()
            while int(data[0])!=cluster_cur and float(data[1])!=widening:
                line=file_alpha.readline()
                data=line.split()
            dispersion_all['alpha_max']=np.array(data[2:]).astype('float')
            file_alpha.close()

            outputs_all={}
            numdis=dispersion_all['numdis']
            alpha_sum=np.zeros((numdis))
            num_models=0
            exp_sum=np.zeros_like(alpha_sum)
            for function in functions:
                outputs_all[function.__name__]=function('',{},model_ref,{'cluster':0},params_dispersion,params_inversion,dispersion_ref,dispersion_all,0.,np.zeros((1)),outputs={},init=True)

            # for this cluster and widening, get the reference dispersion curves
            dispersion_ref_cluster={}
            dispersion_ref_cluster['R']={}
            dispersion_ref_cluster['L']={}

            i_cluster=np.where(params_dispersion['clusters']==cluster_cur)[0][0]

            dispersion_ref_cluster['R']['dispersion']=dispersion_ref['R']['dispersion'][i_cluster,:]
            dispersion_ref_cluster['R']['error']=dispersion_ref['R']['error'][i_cluster,:]
            dispersion_ref_cluster['R']['error_sum']=dispersion_ref['R']['error_sum'][i_cluster]

            dispersion_ref_cluster['L']['dispersion']=dispersion_ref['L']['dispersion'][i_cluster,:]
            dispersion_ref_cluster['L']['error']=dispersion_ref['L']['error'][i_cluster,:]
            dispersion_ref_cluster['L']['error_sum']=dispersion_ref['L']['error_sum'][i_cluster]

            dispersion_ref_cluster['alpha_max']=dispersion_ref['alpha_max'][cluster_cur][widening]


            dispersion_all_cluster = copy.deepcopy(dispersion_all)
            dispersion_all_cluster['alpha_max']=dispersion_all['alpha_max'][cluster_cur][widening]


            for i in range(len(files)):

                file=files[i]
                cluster=cluster_cur
                (alpha_sum_prop,num_models_prop,exp_sum_prop)=postprocess_util.process_one_file(file,cluster,functions,params_inversion,params_dispersion,dispersion_ref_cluster,dispersion_all_cluster,model_ref,outputs_all)

                alpha_sum+=alpha_sum_prop
                num_models+=num_models_prop
                exp_sum+=exp_sum_prop
                print('done',rank,file)

            if rank!=0:
                comm.Send([alpha_sum,MPI.FLOAT],dest=0,tag=101)
                comm.Send([np.array([num_models]),MPI.INT],dest=0,tag=103)
                comm.Send([exp_sum,MPI.FLOAT],dest=0,tag=104)

            elif rank==0:
                for i in range(size)[1:]:
                    alpha_sum_prop=np.zeros_like(alpha_sum)
                    comm.Recv([alpha_sum_prop,MPI.FLOAT],i,101)
                    alpha_sum+=alpha_sum_prop

                    num_models_prop=np.zeros_like(num_models)
                    comm.Recv([num_models_prop,MPI.FLOAT],i,103)
                    num_models+=num_models_prop

                    exp_sum_prop=np.zeros_like(exp_sum)
                    comm.Recv([exp_sum_prop,MPI.FLOAT],i,104)
                    exp_sum+=exp_sum_prop

            if rank==0:
                kl_dist=alpha_sum/num_models
                file_kl.write(str(cluster_cur)+' '+str(widening)+' '+' '.join(kl_dist[:10].astype(str))+'\n')
                kl_weight=postprocess_util.kl_weights(kl_dist)
                kl_sum+=kl_weight
                outputs_all_final['kl_dists'][cluster_cur][widening]=kl_dist
                outputs_all_final['kl_weight'][cluster_cur][widening]=kl_weight



            for function in outputs_all:
                # copying over the nonstack keys
                for key in outputs_all[function]['stack']:
                    if rank!=0:
                        if type(outputs_all[function]['stack'][key])==float:
                            outputs_all[function]['stack'][key]=np.array([outputs_all[function]['stack'][key]])
                        comm.Send([outputs_all[function]['stack'][key],MPI.FLOAT],dest=0,tag=3)

                    if rank==0:
                        outputs_all[function]['stack'][key][...,:]/=exp_sum
                        for i in range(size)[1:]:

                            outputs_all_prop=np.zeros_like(outputs_all[function]['stack'][key])
                            comm.Recv([outputs_all_prop,MPI.FLOAT],i,3)
                            outputs_all_prop[...,:]/=exp_sum
                            if key in outputs_all[function]['stack']:
                                outputs_all[function]['stack'][key]+=outputs_all_prop
                            else:
                                outputs_all[function]['stack'][key]=outputs_all_prop

            if rank==0:
                for function in outputs_all:
                    if 'meansum_sq' in outputs_all[function].keys():
                        for i in range(len(outputs_all[function]['meansum_sq'])): # calculate standart deviation

                            key_mean_sq=outputs_all[function]['meansum_sq'][i]

                            key_mean_shift=outputs_all[function]['meansum_sq_shift'][i]

                            mean_sq=outputs_all[function]['stack'][key_mean_sq] # mean of squares
                            mean_shifted=np.square(outputs_all[function]['stack'][key_mean_shift]) # square of means, shifted by reference model
                            outputs_all[function]['stack'][key_mean_sq]=np.sqrt(mean_sq-mean_shifted)

                            del outputs_all[function]['stack'][key_mean_shift]

            for function in outputs_all:
                for key in outputs_all[function]['nostack']:
                    outputs_all_final[function]['nostack'][key]=outputs_all['function']['nostack'][key]
                for key in outputs_all[function]['stack']:
                    if key in outputs_all_final[function]['stack']:
                        outputs_all_final[function]['stack'][key]+=kl_weight*outputs_all[function]['stack'][key]

    if rank==0:
        for i_dir in range(len(directories)):
            directory=directories[i_dir]
            cluster_cur=clusters_in[i_dir]
            for i_wide in range(len(widenings)):
                widening=widenings[i_wide]
                outputs_all_final['kl_weight'][cluster_cur][widening]/=kl_sum

        for function in outputs_all:
            for key in outputs_all[function]['stack']:
                if key in outputs_all_final[function]['stack']:
                    outputs_all_final[function]['stack'][key]/=kl_sum

    #outputs_all_final = comm.bcast(outputs_all_final, root=0)

    return outputs_all_final



# def process_one_file(file,cluster,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,outputs_all={}):
#     '''
#     Processes one file, reading the podels in the file, applying the functions to them and stacking the results

#     Parameters
#     ----------
#     input_dict : dict
#         input dictionary containing the file name (file), the list of functions (functions),
#         params_inversion, params_dispersion, dispersion_ref, dispersion_all and model_ref.

#     Returns
#     -------
#     outputs : dict
#         has one subdict for each function, called with the function name.
#         Each of those has a subdict 'stack' and a subdict 'nostack' which will respectively stacked and kept as in the first model

#     '''

#     time_read=0
#     time_alpha=0
#     time_apply=0

#     numdis=dispersion_all['numdis']
#     numtot=0
#     alpha_sum=np.zeros((numdis))
#     exp_sum=np.zeros((numdis))

#     # read file
#     f=open(file,'r')
#     print(file)

#     fline=f.readline()
#     if not fline:
#         print('file ' +file + ' empty')
#         return
#     widening=float(f.readline().split()[2])
#     f.readline()
#     f.readline()
#     f.readline()
#     f.readline()
#     f.readline()
#     f.readline()
#     f.readline()

#     line=f.readline()
#     while line:

#         t1=time.time()
#         dispersion_one={}
#         dispersion_one['R']={}
#         dispersion_one['L']={}

#         model={}

#         data=line.split()
#         npt_true=int(data[0])
#         model['npt_true']=npt_true
#         model['npt']=int(data[1])
#         model['npt_ani']=int(data[2])
#         data=f.readline().split()
#         dispersion_one['R']['Ad']=float(data[0])
#         dispersion_one['L']['Ad']=float(data[1])
#         model['Ad_R']=dispersion_one['R']['Ad']
#         model['Ad_L']=dispersion_one['L']['Ad']

#         # new output
#         # d=np.array(f.readline().split()).astype('float')
#         # vsv=np.array(f.readline().split()).astype('float')
#         # xi=np.array(f.readline().split()).astype('float')
#         # vp=np.array(f.readline().split()).astype('float')

#         # old output
#         d=np.zeros((npt_true))
#         vsv=np.zeros_like(d)
#         xi=np.zeros_like(d)
#         vp=np.zeros_like(d)
#         for i in range(npt_true):
#             data=f.readline().split()
#             d[i]=float(data[0])
#             vsv[i]=float(data[1])
#             xi[i]=float(data[2])
#             vp[i]=float(data[3])
#         model['depth']=d
#         model['vsv']=vsv
#         model['xi']=xi
#         model['vp']=vp

#         dispersion_one['like_w']=float(f.readline())
#         model['like_w']=dispersion_one['like_w']
#         dispersion_one['widening']=widening
#         dispersion_one['cluster']=cluster

#         f.readline()
#         dispersion_R_one=np.array(f.readline().split()).astype('float')
#         dispersion_one['R']['dispersion']=dispersion_R_one

#         f.readline()
#         dispersion_L_one=np.array(f.readline().split()).astype('float')

#         dispersion_one['L']['dispersion']=dispersion_L_one

#         t2=time.time()
#         time_read+=t2-t1

#         if len(dispersion_L_one)==0:
#             print('read: ',time_read)
#             print('get alpha: ',time_alpha)
#             print('apply: ',time_apply)
#             return
#         t1=time.time()

#         alpha_ref,alpha=get_alpha(dispersion_one,cluster,params_dispersion,dispersion_ref,dispersion_all)

#         t2=time.time()
#         time_alpha+=t2-t1

#         t1=time.time()

#         # apply functions
#         for function in functions:
#             function(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs=outputs_all[function.__name__])

#         line=f.readline()
#         numtot+=1
#         alpha_sum+=alpha
#         exp_sum+=np.exp(alpha)

#         t2=time.time()
#         time_apply+=t2-t1

#         # if numtot==10:
#         #     f.close()
#         #     print('read: ',time_read)
#         #     print('get alpha: ',time_alpha)
#         #     print('apply: ',time_apply)

#         #     return alpha_sum,numtot,exp_sum,exp_ssq

#     f.close()

#     print('read: ',time_read)
#     print('get alpha: ',time_alpha)
#     print('apply: ',time_apply)

#     return alpha_sum,numtot,exp_sum

# def kl_weights(kl):
#     return np.exp(kl)

# def create_posterior(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs={},init=False):
#     '''
#     example of a function to be applied to the data
#     creates a 2D histogram of vsv for the reference and each member of the cluster
#     The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

#     Parameters
#     ----------
#     model : dict

#     model_ref : dict

#     dispersion_one : dict

#     params_dispersion : dict

#     params_inversion : dict

#     dispersion_ref : dict

#     dispersion_all : dict

#     alpha_ref : float

#     alpha : numpy array, length numdis

#     Returns
#     -------
#     outputs : dict
#         has 2 subdicts, 'stack' and 'nostack'.

#     '''
#     ndatad=50
#     ndatav=50

#     if init:

#         outputs={}
#         outputs['stack']={}
#         outputs['nostack']={}

#         numdis=dispersion_all['numdis']
#         outputs['stack']['vsv_all']=np.zeros((ndatad,ndatav,numdis))
#         outputs['stack']['xi_all']=np.zeros((ndatad,ndatav,numdis))
#         outputs['stack']['vp_all']=np.zeros((ndatad,ndatav,numdis))

#         outputs['nostack']['depths']=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)
#         outputs['nostack']['vels_vsv']=np.linspace(np.amin(model_ref['vsv']*(1-params_inversion['width_vsv'])),
#                          np.amax(model_ref['vsv']*(1+params_inversion['width_vsv'])),
#                          ndatav)
#         outputs['nostack']['vels_vp']=np.linspace(params_inversion['vpvs_min'],params_inversion['vpvs_max'],ndatav)
#         outputs['nostack']['vels_xi']=np.linspace(params_inversion['xi_min'],params_inversion['xi_max'],ndatav)

#         outputs['nostack']['ndatav']=ndatav

#         return outputs

#     depths=outputs['nostack']['depths']
#     vels_vsv=outputs['nostack']['vels_vsv']
#     vels_vp=outputs['nostack']['vels_vp']
#     vels_xi=outputs['nostack']['vels_xi']

#     vsv_model=model['vsv']
#     xi_model=model['xi']
#     vp_model=model['vp']
#     depth_model=model['depth']
#     alpha_max=dispersion_all['alpha_max']
#     alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#     alpha=np.exp(alpha)

#     ind_vsv=np.digitize(np.interp(depths,depth_model,vsv_model),bins=vels_vsv,right=True)
#     outputs['stack']['vsv_all'][np.arange(ndatad),ind_vsv,:]+=np.tile(alpha,(ndatad,1))

#     ind_xi=np.digitize(np.interp(depths,depth_model,xi_model),bins=vels_xi,right=True)
#     outputs['stack']['xi_all'][np.arange(ndatad),ind_xi,:]+=np.tile(alpha,(ndatad,1))

#     ind_vp=np.digitize(np.interp(depths,depth_model,vp_model),bins=vels_vp,right=True)
#     outputs['stack']['vp_all'][np.arange(ndatad),ind_vp,:]+=np.tile(alpha,(ndatad,1))

#     return

# def get_average(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs={},init=False):

#     ndatad=200

#     if init:
#         outputs={}
#         outputs['stack']={}
#         outputs['nostack']={}
#         numdis=dispersion_all['numdis']

#         depths=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)

#         outputs['nostack']['depths']=depths
#         outputs['nostack']['ndata']=ndatad

#         outputs['meansum']=[]
#         outputs['meansum_weights']=[]
#         outputs['meansum_ndata']=[]
#         outputs['meansum_sq']=[]
#         outputs['meansum_sq_weights']=[]
#         outputs['meansum_sq_shift']=[]
#         outputs['stack']['alpha_sum']=np.zeros((numdis))
#         outputs['stack']['probani']=np.zeros((ndatad,numdis))
#         outputs['meansum'].append('probani')
#         outputs['meansum_weights'].append('alpha_sum')
#         outputs['meansum_ndata'].append('ndata')

#         return outputs

#     depths=outputs['nostack']['depths']
#     ndatad=outputs['nostack']['ndata']

#     xi_model=model['xi']
#     depth_model=model['depth']
#     alpha_max=dispersion_all['alpha_max']
#     alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#     alpha=np.exp(alpha)

#     xi=np.interp(depths,depth_model,xi_model)
#     alpha_alt,probani_alt=np.meshgrid(alpha,(xi!=1.))
#     outputs['stack']['probani']+=np.multiply(probani_alt,alpha_alt)
#     outputs['stack']['alpha_sum']+=alpha
#     return #outputs

# def get_histograms(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs={},init=False):
#     '''
#     example of a function to be applied to the data
#     creates a 2D histogram of vsv for the reference and each member of the cluster
#     The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

#     Parameters
#     ----------
#     model : dict

#     model_ref : dict

#     dispersion_one : dict

#     params_dispersion : dict

#     params_inversion : dict

#     dispersion_ref : dict

#     dispersion_all : dict

#     alpha_ref : float

#     alpha : numpy array, length numdis

#     Returns
#     -------
#     outputs : dict
#         has 2 subdicts, 'stack' and 'nostack'.

#     '''
#     ndata=200
#     numdis=dispersion_all['numdis']
#     if init:
#         outputs={}
#         outputs['stack']={}
#         outputs['nostack']={}

#         range_nlay=np.arange(params_inversion['milay'],params_inversion['malay']+1)
#         range_sigmaR=np.linspace(params_inversion['Ad_R_min'],
#                          params_inversion['Ad_R_max'],
#                          ndata)
#         range_sigmaL=np.linspace(params_inversion['Ad_L_min'],
#                          params_inversion['Ad_L_max'],
#                          ndata)
#         range_alpha=np.linspace(-100,
#                          10,
#                          ndata)

#         outputs['stack']['nlay_hist']=np.zeros((params_inversion['malay']-params_inversion['milay']+1,numdis))
#         outputs['stack']['sigmaR_hist']=np.zeros((ndata,numdis))
#         outputs['stack']['sigmaL_hist']=np.zeros((ndata,numdis))
#         outputs['stack']['alpha_hist']=np.zeros((ndata,numdis))

#         outputs['nostack']['range_nlay']=range_nlay
#         outputs['nostack']['range_sigmaR']=range_sigmaR
#         outputs['nostack']['range_sigmaL']=range_sigmaL
#         outputs['nostack']['range_alpha']=range_alpha
#         return outputs

#     range_nlay=outputs['nostack']['range_nlay']
#     range_sigmaR=outputs['nostack']['range_sigmaR']
#     range_sigmaL=outputs['nostack']['range_sigmaL']
#     range_alpha=outputs['nostack']['range_alpha']

#     nlay_model=model['npt']

#     alpha_max=dispersion_all['alpha_max']
#     alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#     alpha_old=alpha
#     alpha=np.exp(alpha)

#     ind_nlay=nlay_model-params_inversion['milay']
#     outputs['stack']['nlay_hist'][ind_nlay,:]+=alpha
#     ind_sigmaR=np.digitize(model['Ad_R'],bins=range_sigmaR)
#     outputs['stack']['sigmaR_hist'][ind_sigmaR,:]+=alpha
#     ind_sigmaL=np.digitize(model['Ad_L'],bins=range_sigmaL)
#     outputs['stack']['sigmaL_hist'][ind_sigmaL,:]+=alpha
#     ind_alpha=np.digitize(alpha_old,bins=range_alpha)
#     outputs['stack']['alpha_hist'][ind_alpha,np.arange(numdis)]+=1
#     return #outputs

# def get_dispersion_mean(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs={},init=False):

#     numdis=dispersion_all['numdis']
#     cluster=dispersion_one['cluster']
#     i_cluster=np.where(params_dispersion['clusters']==cluster)[0][0]
#     if init:
#         outputs={}
#         outputs['stack']={}
#         outputs['nostack']={}

#         outputs['nostack']['clusters']=params_dispersion['clusters']



#         dispersion_R_ref=dispersion_ref['R']['dispersion'][i_cluster,:]
#         dispersion_L_ref=dispersion_ref['L']['dispersion'][i_cluster,:]

#         outputs['nostack']['ndata_R']=len(dispersion_R_ref)
#         outputs['nostack']['ndata_L']=len(dispersion_L_ref)

#         outputs['meansum']=[]
#         outputs['meansum_weights']=[]
#         outputs['meansum_ndata']=[]
#         outputs['meansum_sq']=[]
#         outputs['meansum_sq_weights']=[]
#         outputs['meansum_sq_ndata']=[]
#         outputs['meansum_sq_shift']=[]

#         shape_all_R=np.shape(dispersion_all['R']['dispersion'])
#         shape_all_L=np.shape(dispersion_all['L']['dispersion'])

#         outputs['stack']['dispersion_R']=np.zeros(shape_all_R)
#         outputs['stack']['dispersion_R_sq']=np.zeros(shape_all_R)
#         outputs['stack']['dispersion_R_shift']=np.zeros(shape_all_R)
#         outputs['meansum'].append('dispersion_R')
#         outputs['meansum_ndata'].append('ndata_R')
#         outputs['meansum_sq'].append('dispersion_R_sq')
#         outputs['meansum_sq_ndata'].append('ndata_R')
#         outputs['meansum_sq_shift'].append('dispersion_R_shift')

#         outputs['stack']['dispersion_L']=np.zeros(shape_all_L)
#         outputs['stack']['dispersion_L_sq']=np.zeros(shape_all_L)
#         outputs['stack']['dispersion_L_shift']=np.zeros(shape_all_L)
#         outputs['meansum'].append('dispersion_L')
#         outputs['meansum_ndata'].append('ndata_L')
#         outputs['meansum_sq'].append('dispersion_L_sq')
#         outputs['meansum_sq_ndata'].append('ndata_L')
#         outputs['meansum_sq_shift'].append('dispersion_L_shift')

#         return outputs

#     dispersion_R_model=dispersion_one['R']['dispersion']
#     dispersion_L_model=dispersion_one['L']['dispersion']

#     alpha_max=dispersion_all['alpha_max']
#     alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#     alpha=np.exp(alpha)
#     dispersion_R_true_ref=dispersion_ref['R']['dispersion'][i_cluster,:]
#     dispersion_L_true_ref=dispersion_ref['L']['dispersion'][i_cluster,:]
#     dispersion_R_true_all=dispersion_all['R']['dispersion']
#     dispersion_L_true_all=dispersion_all['L']['dispersion']

#     alpha_alt=np.tile(alpha,(len(dispersion_R_true_ref),1))
#     outputs['stack']['dispersion_R']+=np.multiply(dispersion_R_true_all,alpha_alt)
#     dis_alt=np.transpose(np.tile(dispersion_R_model,(numdis,1)))-dispersion_R_true_all
#     outputs['stack']['dispersion_R_shift']+=np.multiply(dis_alt,alpha_alt)
#     dis_alt=np.square(np.transpose(np.tile(dispersion_R_model,(numdis,1)))-dispersion_R_true_all)
#     outputs['stack']['dispersion_R_sq']+=np.multiply(dis_alt,alpha_alt)

#     alpha_alt=np.tile(alpha,(len(dispersion_L_true_ref),1))
#     outputs['stack']['dispersion_L']+=np.multiply(dispersion_R_true_all,alpha_alt)
#     dis_alt=np.transpose(np.tile(dispersion_L_model,(numdis,1)))-dispersion_L_true_all
#     outputs['stack']['dispersion_L_shift']+=np.multiply(dis_alt,alpha_alt)
#     dis_alt=np.square(np.transpose(np.tile(dispersion_L_model,(numdis,1)))-dispersion_L_true_all)
#     outputs['stack']['dispersion_L_sq']+=np.multiply(dis_alt,alpha_alt)

#     return #outputs

# def get_tradeoff(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_ref,dispersion_all,alpha_ref,alpha,outputs={},init=False):

#     if init:
#         outputs={}
#         outputs['stack']={}
#         outputs['nostack']={}

#         numdis=dispersion_all['numdis']

#         malay=params_inversion['malay']

#         outputs['nostack']['range']=np.arange(malay)

#         outputs['stack']['tradeoff']=np.zeros((malay+1,malay+1,numdis))

#         return outputs


#     alpha_max=dispersion_all['alpha_max']
#     alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#     alpha=np.exp(alpha)

#     malay=params_inversion['malay']

#     nlay=model['npt']
#     nlay_ani=model['npt_ani']
#     outputs['stack']['tradeoff'][nlay,nlay_ani,:]+=alpha
#     return #outputs


############################################################################
#        Example
############################################################################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# directories=['OUT_CLUSTER0_ALL/OUT']
# clusters=[0]
# maxpercent=0.05
# prepare=True
# params_dispersion, dispersion_ref, dispersion_all=get_dispersion(directories,clusters)
# model_ref=get_model_ref()
# output={}
# params_inversion={}
# widening=1.
# params_inversion=get_metadata(directories,clusters,prepare=prepare,widening=widening)
# get_alpha_max(comm,directories,clusters, params_inversion, params_dispersion, dispersion_ref, dispersion_all, prepare=prepare, widening=widening, maxpercent=maxpercent)

# if rank==0:
#     print(dispersion_all['alpha_max'][:100])

# # #directory='OUT_CLUSTER0_ALL/OUT'
# # directory='OUT_CLUSTER0_ALL/OUT'
# # directories=['OUT_CLUSTER0_ALL/OUT','OUT_CLUSTER1_ALL/OUT']
# # clusters=[0,1]
# # maxpercent=0.05
# # prepare=True
# # params_dispersion, dispersion_ref, dispersion_all=get_dispersion(directories,clusters)
# # model_ref=get_model_ref()
# # output={}
# # params_inversion={}
# # widening=1.
# # params_inversion=get_metadata(directories,clusters,prepare=prepare,widening=widening)
# # get_alpha_max(comm,directories,clusters, params_inversion, params_dispersion, dispersion_ref, dispersion_all, prepare=prepare, widening=widening, maxpercent=maxpercent)
# # print(dispersion_all['alpha_max'][:100])

# sys.exit()


directories=['OUT_CLUSTER0_ALL/OUT','OUT_CLUSTER1_ALL/OUT']
clusters=[0,1]
widenings=np.arange(1.,11.)
datadir='OUT_KL/Processing/'

functions=[postprocess_functions.create_posterior,postprocess_functions.get_average,postprocess_functions.get_histograms,postprocess_functions.get_dispersion_mean,postprocess_functions.get_tradeoff] #create_posterior,get_average,get_histograms,,get_tradeoff

print('start')
params_dispersion, dispersion_ref, dispersion_all=postprocess_util.get_dispersion(directories,clusters)
print('get dispersion')
model_ref=postprocess_util.get_model_ref()
print('get model ref')

if rank==0:
    filename='dispersion.h5'
    f = h5py.File(datadir+filename,'w')
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

#for maxpercent in [0.,0.01,0.05,0.1,0.5,1.]:
for maxpercent in [0.,0.01,0.05,0.1,0.5,1.]:
    comm.Barrier()
    alphafile=datadir+'alpha_max'+str(maxpercent)+'.txt'
    kl_file=datadir+'kl_file'+str(maxpercent)+'.txt'
    output={}
    params_inversion={}
    params_inversion=postprocess_util.get_metadata(directories,widenings,clusters)
    print('get metadata')
    print('calculating alphamax')
    get_alpha_max(comm,directories,widenings,clusters, params_inversion, params_dispersion, dispersion_ref, dispersion_all, maxpercent=maxpercent,alphafile=alphafile)

for maxpercent in [0.,0.01,0.05,0.1,0.5,1.]:
    comm.Barrier()
    alphafile=datadir+'alpha_max'+str(maxpercent)+'.txt'
    kl_file=datadir+'kl_file'+str(maxpercent)+'.txt'
    output={}
    params_inversion={}
    params_inversion=postprocess_util.get_metadata(directories,widenings,clusters)
    print('get metadata')
    get_alpha_max(comm,directories,widenings,clusters, params_inversion, params_dispersion, dispersion_ref, dispersion_all, maxpercent=maxpercent,alphafile=alphafile)

    comm.Barrier()
    output=apply_stuff(comm,directories,clusters,widenings,functions,params_inversion,params_dispersion,dispersion_ref,dispersion_all,model_ref,alphafile=alphafile,kl_file=kl_file)
    print('apply functions')

    if rank==0:

        ############################################################
        # print to files
        ############################################################
        for function in output:
            if function=='kl_dists' or function=='kl_weight':
                filename='processing_'+str(maxpercent)+'_'+function+'_outputs.h5'
                print(filename)
                f = h5py.File(datadir+filename,'w')
                grp=f.create_group(function)
                for key in output[function]:
                    grp_stack=grp.create_group(key)
                    for key2 in output[function][key]:
                        grp_nostack=grp_stack.create_dataset(key2,data=output['kl_dists'][key][key2])
                f.close()
                output[function]={}

            else:
                filename='processing_'+str(maxpercent)+'_'+function+'_outputs.h5'
                print(filename)
                f = h5py.File(datadir+filename,'w')

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

    del params_inversion,output


