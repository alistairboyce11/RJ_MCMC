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


def get_alpha_max(comm,directories,widenings_in,clusters_in,params_dispersion,dispersion_ref,dispersion_all,alphafile='OUT_KL/Processing/alpha_max.txt',maxpercent=0.):
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

    numdis=dispersion_all['numdis']

    num_models_dict={}

    if os.path.isfile(alphafile):
        file=open(alphafile,'r')
        dispersion_ref['alpha_max']=0.
        dispersion_all['alpha_max']=np.array(file.readline().split()).astype('float')
        for line in file:
            #print(line)
            data=line.split()
            cluster=int(data[0])
            widening=float(data[1])
            num=int(data[2])
            if cluster in num_models_dict:
                num_models_dict[cluster][widening]=num
            else:
                num_models_dict[cluster]={}
                num_models_dict[cluster][widening]=num
        file.close()
        params_dispersion['num_models']=num_models_dict
        return

    files_all=[]

    # initiate stuff
    for i_dir in range(len(directories)):
        directory=directories[i_dir]
        cluster_cur=clusters_in[i_dir]
        num_models_dict[cluster_cur]={}
        for widening in widenings_in:
            num_models_dict[cluster_cur][widening]=0
            # list files
            files_all_tmp=glob.glob(directory+'/All_models_processed_prepare_*%4.2f.out'%widening)
            files_all_tmp.sort()
            files_all.extend(files_all_tmp)


    # create array for max alpha, potentially more than one if maxpercent is bigger than zero
    num_to_store=0
    if maxpercent==0.:
        num_to_store=1
    else:
        for file in files_all:
            f=open(file,'r')
            params_inversion=postprocess_util.read_header(f)
            f.close()
            num_to_store+=maxpercent/100.*params_inversion['everyall']
            #print(maxpercent/100.*params_inversion['everyall'],num_to_store)

        num_to_store=round(num_to_store)
        print('num to store',num_to_store)

    numdis=dispersion_all['numdis']

    #alpha_ref_max=np.ones((num_to_store))*float('-inf')
    alpha_max=np.ones((num_to_store,numdis))*float('-inf')

    # split it over the cores
    size = comm.Get_size()
    rank = comm.Get_rank()
    rank2=rank

    if rank==0:
        print(files_all)

    whodoeswhat={}
    for i in range(len(files_all)):
        whodoeswhat[files_all[i]]=i%size

    num_models_all=0

    for i_dir in range(len(directories)):
        directory=directories[i_dir]
        cluster=clusters_in[i_dir]
        for widening in widenings_in:

            files_all_tmp=glob.glob(directory+'/All_models_processed_prepare_*%4.2f.out'%widening)
            #files_all_tmp=['OUT_CLUSTER0_START/OUT/All_models_processed_prepare_0000_ 1.00.out']
            files_all_tmp.sort()
            files=[]
            for i in range(len(files_all_tmp)):
                if i%size==rank2:
                    files.append(files_all_tmp[i])

            # if more cores than files: go to the next iteration of the loop
            # but keep core 0 in, it's the core that gathers everything
            computing=[]
            if len(files)==0:
                rank2-=len(files_all_tmp)
                if rank!=0:
                    continue
            for f in files_all_tmp:
                try:
                    computing.append(whodoeswhat[f])
                except:
                    pass

            # alpha_ref_max_tmp=np.ones((num_to_store))*float('-inf')
            alpha_max_tmp=np.ones((num_to_store,numdis))*float('-inf')
            num_models_tmp=0

            if rank==0:
                print(computing)

            #print(files)


            # loop over the files
            for i in range(len(files)):

                file=files[i]

                # process one file
                print('computing alphamax for file',file,'rank',rank)
                alpha_max_prop,alpha_ref_max_prop,num_models_prop=postprocess_util.process_file_alphamax_all(file, clusters_in,widenings_in,cluster,widening,params_dispersion,dispersion_ref,dispersion_all,maxpercent, num_to_store)
                #print(num_models_prop)
                #print('here')

                # merge results on one core
                alpha_max_tmp=np.append(alpha_max_tmp,alpha_max_prop,axis=0)
                alpha_max_tmp.sort(axis=0)
                alpha_max_tmp=alpha_max_tmp[-num_to_store:,:]

                num_models_tmp+=num_models_prop
            # print(computing)
            computing[:] = [x for x in computing if x != 0]
# =============================================================================
#             started=np.zeros_like(files_all_tmp)
#             computing=np.zeros_like(files_all_tmp)
#
#             # alpha_ref_max_tmp=np.ones((num_to_store))*float('-inf')
#             alpha_max_tmp=np.ones((num_to_store,numdis))*float('-inf')
#             num_models_tmp=0
#
#             while np.any(started==0):
#
#                 index=np.where(started==0)[0]
#                 print(index)
#                 file=files_all_tmp[index]
#                 started[index]=1
#                 started=comm.bcast(started, root=rank)
#                 computing[index]=rank
#                 computing=comm.bcast(computing, root=rank)
#
#                 # process one file
#                 print('computing alphamax for file',file,'rank',rank)
#                 alpha_max_prop,alpha_ref_max_prop,num_models_prop=postprocess_util.process_file_alphamax_all(file, clusters_in,widenings_in,cluster,widening,params_dispersion,dispersion_ref,dispersion_all,maxpercent, num_to_store)
#                 #print(num_models_prop)
#                 #print('here')
#
#                 # merge results on one core
#                 alpha_max_tmp=np.append(alpha_max_tmp,alpha_max_prop,axis=0)
#                 alpha_max_tmp.sort(axis=0)
#                 alpha_max_tmp=alpha_max_tmp[-num_to_store:,:]
#
#                 num_models_tmp+=num_models_prop
#             # print(computing)
#             computing2=[]
#             computing2[:] = [x for x in computing if x != 0]
#             # if 0 in computing:
#             #     computing.remove(0) # we gather on 0, so don't send 0 to 0
#             #print(computing)
# =============================================================================
            # merge results between cores on core 0
            if rank in np.unique(computing):
                print('sending',rank)
                comm.Ssend([alpha_max_tmp,MPI.FLOAT],dest=0,tag=1+int(widening)*100+int(cluster)*1000)
                comm.Ssend([np.array([num_models_tmp]),MPI.INT],dest=0,tag=3+int(widening)*100+int(cluster)*1000)
                print('sent',rank)
            if rank==0:
                # recieve and merge data from all other cores
                #print(num_models_tmp)
                for i in np.unique(computing):
                    print('recieving',i)
                    alpha_max_prop=np.zeros_like(alpha_max_tmp)
                    comm.Recv([alpha_max_prop,MPI.FLOAT],i,1+int(widening)*100+int(cluster)*1000)
                    alpha_max_tmp=np.append(alpha_max_tmp,alpha_max_prop,axis=0)
                    alpha_max_tmp.sort(axis=0)
                    alpha_max_tmp=alpha_max_tmp[-num_to_store:,:]
                    num_models_prop=np.zeros_like(num_models_tmp)
                    comm.Recv([num_models_prop,MPI.INT],i,3+int(widening)*100+int(cluster)*1000)
                    #print(num_models_prop)
                    num_models_tmp+=num_models_prop

                    print('recieved',i)

                # normalise alphamax by number of models for this cluster and widening
                print(num_models_tmp)
                if widening in num_models_dict[cluster]:
                    num_models_dict[cluster][widening]+=num_models_tmp
                else:
                    num_models_dict[cluster][widening]=num_models_tmp
                #if len(files)>0:
                if num_models_tmp==0:
                    alpha_max_tmp=np.ones_like(alpha_max_tmp)*float('-inf')

                else:

                    alpha_max_tmp-=np.log(num_models_tmp)

                # take max over all clusters/widenings (incremental)
                alpha_max=np.append(alpha_max,alpha_max_tmp,axis=0)
                alpha_max.sort(axis=0)
                alpha_max=alpha_max[-num_to_store:,:]
                num_models_all+=num_models_tmp


                print('computed alphamax for cluster',cluster,'and widening',widening)
            # special treatement for core 0, else
            # shift to take into account we have more cores than files (potentially)
            if rank!=0 or len(files)!=0:
                rank2+=size*len(files)-len(files_all_tmp)
    num_models_all = comm.bcast(num_models_all, root=0)
    if num_models_all==0:
        num_models_dict = comm.bcast(num_models_dict, root=0)
        params_dispersion['num_models']=num_models_dict
        return
    if rank==0:
        # handle alphamax on core 0
        # update to real size of alpha_max
        if maxpercent==0.:
            num_to_store=1
        else:
            num_to_store=int(num_models_all*maxpercent/100.)

        # get true alpha_max: min of alpha_max
        print('getting true alpha_max')
        alpha_max.sort(axis=0)
        alpha_max=alpha_max[-num_to_store,:]

        # print alpha_max to file
        print('writing')
        file_alpha=open(alphafile,'w')
        file_alpha.write(' '.join(np.array(alpha_max).astype('str'))+'\n')
        for cluster in num_models_dict:
            for widening in num_models_dict[cluster]:
                    file_alpha.write(' '.join([str(cluster),str(widening),str(num_models_dict[cluster][widening])])+'\n')
        file_alpha.close()
        print('written')

        # broadcast to all other cores
        print('broadcasting')
    alpha_max = comm.bcast(alpha_max, root=0)
    #alpha_ref_max = comm.bcast(alpha_ref_max, root=0)
    num_models_dict = comm.bcast(num_models_dict, root=0)
    #print('broadcasted')
    comm.Barrier()
    # put alpha_max in dispersion_all and dispersion_ref for its cluster and widening
    dispersion_ref['alpha_max']=0.
    dispersion_all['alpha_max']=alpha_max
    params_dispersion['num_models']=num_models_dict

    return

def apply_stuff(comm,directories,widenings_in,clusters_in,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref):
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
    outputs_all={}
    params_inversion=postprocess_util.get_metadata(directories,widenings_in)
    for function in functions:
        #print(params_inversion)
        outputs_all[function.__name__]=function('',{},model_ref,{'cluster':0},params_dispersion,params_inversion,dispersion_all,0.,np.zeros((1)),outputs={},init=True)

    # initialise weights
    numdis=dispersion_all['numdis']
    alpha_sum=np.zeros((numdis+1))
    num_models=0
    exp_sum=np.zeros_like(alpha_sum)

    files_all=[]
    clusters_all=[]
    widenings_all=[]

    for i_dir in range(len(directories)):
        directory=directories[i_dir]
        cluster_cur=clusters_in[i_dir]
        for widening in widenings_in:
            # list files
            files_all.extend(glob.glob(directory+'/All_models_processed_prepare_*%4.2f.out'%widening))
            numfiles=len(glob.glob(directory+'/All_models_processed_prepare_*%4.2f.out'%widening))
            clusters_all+=numfiles*[cluster_cur]
            widenings_all+=numfiles*[widening]

    # split it over the cores
    size = comm.Get_size()
    rank = comm.Get_rank()
    files=[]
    clusters=[]
    widenings=[]
    for i in range(len(files_all)):
        if i%size==rank:
            files.append(files_all[i])
            clusters.append(clusters_all[i])
            widenings.append(widenings_all[i])

    # loop over files
    for i in range(len(files)):

        file=files[i]
        cluster=clusters[i]
        widening=widenings[i]
        print(file,cluster,widening)

        (alpha_sum_prop,num_models_prop,exp_sum_prop)=postprocess_util.process_one_file_all(file,cluster,clusters_in,widenings_in,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref,outputs_all)
        print('done',rank,file)

        alpha_sum+=alpha_sum_prop
        num_models+=num_models_prop
        exp_sum+=exp_sum_prop

    # send and gather alpha, alphasum, num_models
    if rank!=0:
        comm.Ssend([alpha_sum,MPI.FLOAT],dest=0,tag=101)
        comm.Ssend([np.array([num_models]),MPI.INT],dest=0,tag=103)
        comm.Ssend([exp_sum,MPI.FLOAT],dest=0,tag=104)

    elif rank==0:
        for i in range(size)[1:]:
            alpha_sum_prop=np.zeros_like(alpha_sum)
            comm.Recv([alpha_sum_prop,MPI.FLOAT],i,101)
            alpha_sum+=alpha_sum_prop

            num_models_prop=np.zeros_like(num_models)
            comm.Recv([num_models_prop,MPI.INT],i,103)
            num_models+=num_models_prop

            exp_sum_prop=np.zeros_like(exp_sum)
            comm.Recv([exp_sum_prop,MPI.FLOAT],i,104)
            exp_sum+=exp_sum_prop

    # send and gather function outputs
    for function in outputs_all:
        for key in outputs_all[function]['stack']:

            # make every float an array for sending
            if type(outputs_all[function]['stack'][key])==float:
                outputs_all[function]['stack'][key]=np.array([outputs_all[function]['stack'][key]])

            # send to core 0
            if rank!=0:
                # print(function,key)
                comm.Ssend([outputs_all[function]['stack'][key],MPI.FLOAT],dest=0,tag=3)

            # recieve and sum up
            if rank==0:
                for i in range(size)[1:]:
                    outputs_all_prop=np.zeros_like(outputs_all[function]['stack'][key])
                    comm.Recv([outputs_all_prop,MPI.FLOAT],i,3)
                    outputs_all[function]['stack'][key]+=outputs_all_prop

                # normalise by sum of weights
                outputs_all[function]['stack'][key][...,:]/=exp_sum

    # calculate standart deviations
    if rank==0:
        # print(np.amax(exp_sum))
        for function in outputs_all:

            if 'meansum_sq' in outputs_all[function].keys():
                for i in range(len(outputs_all[function]['meansum_sq'])):

                    # appropriate keys
                    key_mean_sq=outputs_all[function]['meansum_sq'][i]
                    key_mean_shift=outputs_all[function]['meansum_sq_shift'][i]

                    # mean of squares, shifted
                    mean_sq=outputs_all[function]['stack'][key_mean_sq]

                    # square of means, shifted
                    mean_shifted=np.square(outputs_all[function]['stack'][key_mean_shift])
                    # print(function,key_mean_sq)
                    # print(mean_sq)
                    # print(mean_shifted)
                    # print(mean_sq-mean_shifted)
                    # print(np.where(mean_sq-mean_shifted<0.))
                    # if key_mean_sq=='dispersion_R_sq':
                    #     continue
                    # print(np.amax(mean_shifted))
                    # print(np.amin(mean_sq-mean_shifted))
                    try:
                        outputs_all[function]['stack'][key_mean_sq]=np.sqrt(mean_sq-mean_shifted)
                    except:
                        pass

                    del outputs_all[function]['stack'][key_mean_shift]

    return outputs_all # careful, it is only on core 0!


############################################################################
#        Example
############################################################################

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
print(size,rank)

# directories=['OUT_CLUSTER0_JOINT/OUT','OUT_CLUSTER1_JOINT/OUT','OUT_CLUSTER2_JOINT/OUT','OUT_CLUSTER6_JOINT/OUT','OUT_CLUSTER8_JOINT/OUT','OUT_CLUSTER9_JOINT/OUT']
# clusters=[0,1,2,6,8,9]
# widenings=[1.,2.,3.,4.,5.]
# #points=[100,1000,10000]
# outdir='OUT_ALL_10'
# points='all'
# maxpercent=10./(len(widenings)*8.*100000.*len(clusters))*100.

# directories=['OUT_CLUSTER0_JOINT/OUT']
# clusters=[0]
# widenings=[1.,2.,3.,4.,5.]
# outdir='OUT_0_ALL'
# points=[]

# directories=['OUT_TEST_REAL_SIGMAV_CONSTANT/OUT']
# directories=['OUT_CLUSTER0_START/OUT','OUT_CLUSTER1_START/OUT','OUT_CLUSTER2_START/OUT','OUT_CLUSTER6_START/OUT','OUT_CLUSTER7_START/OUT','OUT_CLUSTER8_START/OUT','OUT_CLUSTER9_START/OUT']
# clusters=[0,1,2,6,7,8,9]
# directories=['OUT_CLUSTER0_NARROW/OUT','OUT_CLUSTER1_NARROW/OUT','OUT_CLUSTER2_NARROW/OUT','OUT_CLUSTER6_NARROW/OUT','OUT_CLUSTER7_NARROW/OUT','OUT_CLUSTER8_NARROW/OUT','OUT_CLUSTER9_NARROW/OUT','OUT_CLUSTER0_NARROW_HIGHT/OUT','OUT_CLUSTER1_NARROW_HIGHT/OUT','OUT_CLUSTER2_NARROW_HIGHT/OUT','OUT_CLUSTER6_NARROW_HIGHT/OUT','OUT_CLUSTER7_NARROW_HIGHT/OUT','OUT_CLUSTER8_NARROW_HIGHT/OUT','OUT_CLUSTER9_NARROW_HIGHT/OUT']
# clusters=[0,1,2,6,7,8,9,0,1,2,6,7,8,9]
# directories=['OUT_CLUSTER0_START/OUT','OUT_CLUSTER1_START/OUT','OUT_CLUSTER2_START/OUT','OUT_CLUSTER7_START/OUT','OUT_CLUSTER8_START/OUT','OUT_CLUSTER9_START/OUT']
# clusters=[0,1,2,7,8,9]
# directories=['OUT_CLUSTER6_START/OUT']
# clusters=[6]

# directories=['OUT_CLUSTER0_NARROW/OUT','OUT_CLUSTER1_NARROW/OUT','OUT_CLUSTER2_NARROW/OUT','OUT_CLUSTER6_NARROW/OUT','OUT_CLUSTER7_NARROW/OUT','OUT_CLUSTER8_NARROW/OUT','OUT_CLUSTER9_NARROW/OUT','OUT_CLUSTER0_NARROW_HIGHT/OUT','OUT_CLUSTER1_NARROW_HIGHT/OUT','OUT_CLUSTER2_NARROW_HIGHT/OUT','OUT_CLUSTER6_NARROW_HIGHT/OUT','OUT_CLUSTER7_NARROW_HIGHT/OUT','OUT_CLUSTER8_NARROW_HIGHT/OUT','OUT_CLUSTER9_NARROW_HIGHT/OUT']
# clusters=[0,1,2,6,7,8,9,0,1,2,6,7,8,9]

maxpercent=0.
suffix=''


for cluster in [1]:
    directories=['OUT_CLUSTER'+str(cluster)+'_NARROW/OUT','OUT_CLUSTER'+str(cluster)+'_NARROW_HIGHT/OUT']
    clusters=[cluster,cluster]

    widenings=[1.]
    outdir='OUT_NARROW_'+str(cluster)+'_1'+suffix
    if not os.path.exists(outdir) and rank==0:
        os.mkdir(outdir)
    points=np.arange(100)#'all'
    #maxpercent=0.
    batch=100

    functions=[postprocess_functions.create_posterior,postprocess_functions.get_average,postprocess_functions.get_histograms,postprocess_functions.get_dispersion_mean,postprocess_functions.get_tradeoff]

    if rank==0:
        print('start')
    params_dispersion, dispersion_ref, dispersion_all_ref=postprocess_util.get_dispersion(directories,clusters,points=points)
    if rank==0:
        print('get dispersion')
    model_ref=postprocess_util.get_model_ref(filename='Modified_PREM_GLOBAL.in')
    if rank==0:
        print('get model ref')

    # save dispersion to file
    if rank==0:
        filename='dispersion.h5'
        f = h5py.File(outdir+'/'+filename,'w')
        grp=f.create_group('cluster_params')
        grp.create_dataset('batchsize',data=batch)
        grp.create_dataset('numdis',data=dispersion_all_ref['numdis'])
        grp.create_dataset('lat',data=dispersion_all_ref['lat'])
        grp.create_dataset('lon',data=dispersion_all_ref['lon'])
        grp.create_dataset('cluster',data=dispersion_all_ref['cluster'])
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
        grp2.create_dataset('dispersion',data=dispersion_all_ref['L']['dispersion'])
        grp2.create_dataset('error',data=dispersion_all_ref['L']['error'])
        grp2=grp.create_group('R')
        grp2.create_dataset('dispersion',data=dispersion_all_ref['R']['dispersion'])
        grp2.create_dataset('error',data=dispersion_all_ref['R']['error'])
        f.close()

    comm.Barrier()
    #sys.exit()

    #widening=1.0
    i=0
    while i*batch<dispersion_all_ref['numdis']+1:
        points2=[]
        points2=range(i*batch,min([(i+1)*batch,dispersion_all_ref['numdis']]))
        params_dispersion, dispersion_ref, dispersion_all=postprocess_util.get_dispersion(directories,clusters,points=points2)
        output={}
        if rank==0:
            print('get alphamax')
        alphafile=outdir+'/alphamax_'+str(i)+'.txt'
        get_alpha_max(comm,directories,widenings,clusters,params_dispersion,dispersion_ref,dispersion_all,alphafile,maxpercent=maxpercent)
        num_model_all=0
        for widening in widenings:
            num_model_all+=params_dispersion['num_models'][cluster][widening]
        if num_model_all==0:
            i+=1
            continue

        if rank==0:
            print('apply functions')
        output=apply_stuff(comm,directories,widenings,clusters,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref) #

        if rank==0:
            print('here')
            ############################################################
            # print to files
            ############################################################
            for function in output:
                filename='processing_'+'_'+function+'_outputs_'+str(i)+'.h5'
                f = h5py.File(outdir+'/'+filename,'w')
                # f.create_dataset('burn-in',data=params_inversion['burn-in'])
                # f.create_dataset('nsample',data=params_inversion['nsample'])


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
        i+=1

directories=['OUT_CLUSTER0_NARROW/OUT','OUT_CLUSTER1_NARROW/OUT','OUT_CLUSTER2_NARROW/OUT','OUT_CLUSTER6_NARROW/OUT','OUT_CLUSTER7_NARROW/OUT','OUT_CLUSTER8_NARROW/OUT','OUT_CLUSTER9_NARROW/OUT','OUT_CLUSTER0_NARROW_HIGHT/OUT','OUT_CLUSTER1_NARROW_HIGHT/OUT','OUT_CLUSTER2_NARROW_HIGHT/OUT','OUT_CLUSTER6_NARROW_HIGHT/OUT','OUT_CLUSTER7_NARROW_HIGHT/OUT','OUT_CLUSTER8_NARROW_HIGHT/OUT','OUT_CLUSTER9_NARROW_HIGHT/OUT']
clusters=[0,1,2,6,7,8,9,0,1,2,6,7,8,9]
widenings=[1.]
outdir='OUT_NARROW_ALL'+suffix
if not os.path.exists(outdir) and rank==0:
    os.mkdir(outdir)
points=np.arange(100)#'all'
#maxpercent=0.
batch=100

functions=[postprocess_functions.create_posterior,postprocess_functions.get_average,postprocess_functions.get_histograms,postprocess_functions.get_dispersion_mean,postprocess_functions.get_tradeoff]

if rank==0:
    print('start')
params_dispersion, dispersion_ref, dispersion_all_ref=postprocess_util.get_dispersion(directories,clusters,points=points)
if rank==0:
    print('get dispersion')
model_ref=postprocess_util.get_model_ref(filename='Modified_PREM_GLOBAL.in')
if rank==0:
    print('get model ref')

# save dispersion to file
if rank==0:
    filename='dispersion.h5'
    f = h5py.File(outdir+'/'+filename,'w')
    grp=f.create_group('cluster_params')
    grp.create_dataset('batchsize',data=batch)
    grp.create_dataset('numdis',data=dispersion_all_ref['numdis'])
    grp.create_dataset('lat',data=dispersion_all_ref['lat'])
    grp.create_dataset('lon',data=dispersion_all_ref['lon'])
    grp.create_dataset('cluster',data=dispersion_all_ref['cluster'])
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
    grp2.create_dataset('dispersion',data=dispersion_all_ref['L']['dispersion'])
    grp2.create_dataset('error',data=dispersion_all_ref['L']['error'])
    grp2=grp.create_group('R')
    grp2.create_dataset('dispersion',data=dispersion_all_ref['R']['dispersion'])
    grp2.create_dataset('error',data=dispersion_all_ref['R']['error'])
    f.close()

comm.Barrier()
#sys.exit()

#widening=1.0
i=0
while i*batch<dispersion_all_ref['numdis']+1:
    points2=[]
    points2=range(i*batch,min([(i+1)*batch,dispersion_all_ref['numdis']]))
    params_dispersion, dispersion_ref, dispersion_all=postprocess_util.get_dispersion(directories,clusters,points=points2)
    output={}
    if rank==0:
        print('get alphamax')
    alphafile=outdir+'/alphamax_'+str(i)+'.txt'
    get_alpha_max(comm,directories,widenings,clusters,params_dispersion,dispersion_ref,dispersion_all,alphafile,maxpercent=maxpercent)

    num_model_all=0
    for cluster in clusters:
        for widening in widenings:
            num_model_all+=params_dispersion['num_models'][cluster][widening]
    if num_model_all==0:
        i+=1
        continue

    if rank==0:
        print('apply functions')
    output=apply_stuff(comm,directories,widenings,clusters,functions,params_dispersion,dispersion_ref,dispersion_all,model_ref) #

    if rank==0:

        ############################################################
        # print to files
        ############################################################
        for function in output:
            filename='processing_'+'_'+function+'_outputs_'+str(i)+'.h5'
            f = h5py.File(outdir+'/'+filename,'w')
            # f.create_dataset('burn-in',data=params_inversion['burn-in'])
            # f.create_dataset('nsample',data=params_inversion['nsample'])


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
    i+=1
