#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 14:31:30 2022

@author: dorian
"""
import numpy as np

def create_posterior(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha,outputs={},init=False):
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

    dispersion_all : dict

    alpha : numpy array, length numdis

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    ndatad=100
    ndatav=200

    numdis=dispersion_all['numdis']

    if init:

        outputs={}
        outputs['stack']={}
        outputs['nostack']={}


        outputs['stack']['vsv_all']=np.zeros((ndatad,ndatav,numdis+1))
        outputs['stack']['xi_all']=np.zeros((ndatad,ndatav,numdis+1))
        outputs['stack']['vp_all']=np.zeros((ndatad,ndatav,numdis+1))

        outputs['nostack']['depths']=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)
        outputs['nostack']['vels_vsv']=np.linspace(np.amin(model_ref['vsv']*(1-params_inversion['width_vsv'])),
                         np.amax(model_ref['vsv']*(1+params_inversion['width_vsv'])),
                         ndatav) # in km/s
        outputs['nostack']['vels_vp']=np.linspace(np.amin(model_ref['vpv'])*(1+params_inversion['vp_min']),
                         np.amax(model_ref['vpv'])*(1+params_inversion['vp_max']),
                         ndatav) # in km/s
        outputs['nostack']['vels_xi']=np.linspace(params_inversion['xi_min'],params_inversion['xi_max'],ndatav)

        outputs['nostack']['ndatav']=ndatav

        return outputs

    depths=outputs['nostack']['depths']
    vels_vsv=outputs['nostack']['vels_vsv']
    vels_vp=outputs['nostack']['vels_vp']
    vels_xi=outputs['nostack']['vels_xi']

    vsv_model=model['vsv']
    xi_model=model['xi']
    vp_model=model['vp']
    depth_model=model['depth']
    alpha_max=dispersion_all['alpha_max']
    #print(np.shape(alpha),np.shape(alpha_max))
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
#    if np.amax(alpha)<-50:
#        ind_vsv=np.digitize(np.interp(depths,depth_model,vsv_model),bins=vels_vsv,right=True)
#        outputs['stack']['vsv_all'][np.arange(ndatad),ind_vsv,numdis]+=1.
#        ind_xi=np.digitize(np.interp(depths,depth_model,xi_model),bins=vels_xi,right=True)
#        outputs['stack']['xi_all'][np.arange(ndatad),ind_xi,numdis]+=1.
#        ind_vp=np.digitize(np.interp(depths,depth_model,vp_model),bins=vels_vp,right=True)
#        outputs['stack']['vp_all'][np.arange(ndatad),ind_vp,numdis]+=1.
#        return
    alpha=np.exp(alpha)
    #print(numdis)

    ind_vsv=np.digitize(np.interp(depths,depth_model,vsv_model),bins=vels_vsv,right=True)
    outputs['stack']['vsv_all'][np.arange(ndatad),ind_vsv,:numdis]+=np.tile(alpha,(ndatad,1))

    outputs['stack']['vsv_all'][np.arange(ndatad),ind_vsv,numdis]+=1.

    ind_xi=np.digitize(np.interp(depths,depth_model,xi_model),bins=vels_xi,right=True)
    outputs['stack']['xi_all'][np.arange(ndatad),ind_xi,:numdis]+=np.tile(alpha,(ndatad,1))

    outputs['stack']['xi_all'][np.arange(ndatad),ind_xi,numdis]+=1.

    ind_vp=np.digitize(np.interp(depths,depth_model,vp_model),bins=vels_vp,right=True)
    #
    #if np.amax(ind_vp)>=50:
    #    print(params_inversion['vp_max'],np.amax(ind_vp),np.amax(vp_model),np.amax(vels_vp))
    #    pass
    #else:
    outputs['stack']['vp_all'][np.arange(ndatad),ind_vp,:numdis]+=np.tile(alpha,(ndatad,1))

    outputs['stack']['vp_all'][np.arange(ndatad),ind_vp,numdis]+=1.

    return

def get_average(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha,outputs={},init=False):

    ndatad=200
    numdis=dispersion_all['numdis']

    if init:
        outputs={}
        outputs['stack']={}
        outputs['nostack']={}


        depths=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)

        outputs['nostack']['depths']=depths
        outputs['nostack']['ndata']=ndatad

        outputs['stack']['probani']=np.zeros((ndatad,numdis+1))

        outputs['stack']['vsv']=np.zeros((ndatad,numdis+1))
        outputs['stack']['vp']=np.zeros((ndatad,numdis+1))
        outputs['stack']['xi']=np.zeros((ndatad,numdis+1))

        return outputs

    depths=outputs['nostack']['depths']
    ndatad=outputs['nostack']['ndata']

    vsv_model=model['vsv']
    xi_model=model['xi']
    vp_model=model['vp']
    depth_model=model['depth']
    alpha_max=dispersion_all['alpha_max']
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
 #   if np.amax(alpha)<-50:
 #       xi=np.interp(depths,depth_model,xi_model)
 #       outputs['stack']['probani'][:,numdis]+=(xi!=1).astype(int)
 #       outputs['stack']['xi'][:,numdis]+=xi
 #       vp=np.interp(depths,depth_model,vp_model)
 #       outputs['stack']['vp'][:,numdis]+=vp
 #       vsv=np.interp(depths,depth_model,vsv_model)
 #       outputs['stack']['vsv'][:,numdis]+=vsv
 #       return
    alpha=np.exp(alpha)

    xi=np.interp(depths,depth_model,xi_model)
    alpha_alt,probani_alt=np.meshgrid(alpha,(xi!=1.).astype(int))
    outputs['stack']['probani'][:,:numdis]+=np.multiply(probani_alt,alpha_alt)
    outputs['stack']['probani'][:,numdis]+=(xi!=1).astype(int)
    alpha_alt,xi_alt=np.meshgrid(alpha,xi)
    outputs['stack']['xi'][:,:numdis]+=np.multiply(xi_alt,alpha_alt)
    outputs['stack']['xi'][:,numdis]+=xi

    vp=np.interp(depths,depth_model,vp_model)
    alpha_alt,vp_alt=np.meshgrid(alpha,vp)
    outputs['stack']['vp'][:,:numdis]+=np.multiply(vp_alt,alpha_alt)
    outputs['stack']['vp'][:,numdis]+=vp

    vsv=np.interp(depths,depth_model,vsv_model)
    alpha_alt,vsv_alt=np.meshgrid(alpha,vsv)
    outputs['stack']['vsv'][:,:numdis]+=np.multiply(vsv_alt,alpha_alt)
    outputs['stack']['vsv'][:,numdis]+=vsv
    return #outputs

def get_histograms(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha,outputs={},init=False):
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

    dispersion_all : dict

    alpha : numpy array, length numdis

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    ndata=200
    numdis=dispersion_all['numdis']
    if init:
        outputs={}
        outputs['stack']={}
        outputs['nostack']={}

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

        outputs['stack']['nlay_hist']=np.zeros((params_inversion['malay']-params_inversion['milay']+1,numdis+1))
        outputs['stack']['sigmaR_hist']=np.zeros((ndata,numdis+1))
        outputs['stack']['sigmaL_hist']=np.zeros((ndata,numdis+1))
        outputs['stack']['alpha_hist']=np.zeros((ndata,numdis+1))

        outputs['nostack']['range_nlay']=range_nlay
        outputs['nostack']['range_sigmaR']=range_sigmaR
        outputs['nostack']['range_sigmaL']=range_sigmaL
        outputs['nostack']['range_alpha']=range_alpha
        return outputs

    range_nlay=outputs['nostack']['range_nlay']
    range_sigmaR=outputs['nostack']['range_sigmaR']
    range_sigmaL=outputs['nostack']['range_sigmaL']
    range_alpha=outputs['nostack']['range_alpha']

    nlay_model=model['npt']

    alpha_max=dispersion_all['alpha_max']
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    alpha_old=alpha
#    if np.amax(alpha)<-50:
#        ind_nlay=nlay_model-params_inversion['milay']
#        outputs['stack']['nlay_hist'][ind_nlay,numdis]+=1.
#        ind_sigmaR=np.digitize(model['Ad_R'],bins=range_sigmaR)
#        outputs['stack']['sigmaR_hist'][ind_sigmaR,numdis]+=1.
#        ind_sigmaL=np.digitize(model['Ad_L'],bins=range_sigmaL)
#        outputs['stack']['sigmaL_hist'][ind_sigmaL,numdis]+=1.
#        return
    alpha=np.exp(alpha)

    #print(nlay_model,params_inversion['milay'],params_inversion['malay'],nlay_model-params_inversion['milay'])
    ind_nlay=nlay_model-params_inversion['milay']
    outputs['stack']['nlay_hist'][ind_nlay,:numdis]+=alpha
    outputs['stack']['nlay_hist'][ind_nlay,numdis]+=1.

    ind_sigmaR=np.digitize(model['Ad_R'],bins=range_sigmaR)
    outputs['stack']['sigmaR_hist'][ind_sigmaR,:numdis]+=alpha
    outputs['stack']['sigmaR_hist'][ind_sigmaR,numdis]+=1.
    ind_sigmaL=np.digitize(model['Ad_L'],bins=range_sigmaL)
    outputs['stack']['sigmaL_hist'][ind_sigmaL,:numdis]+=alpha
    outputs['stack']['sigmaL_hist'][ind_sigmaL,numdis]+=1.
    ind_alpha=np.digitize(alpha_old,bins=range_alpha)
    outputs['stack']['alpha_hist'][ind_alpha,np.arange(numdis)]+=1.

    return #outputs

def get_dispersion_mean(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha,outputs={},init=False):

    numdis=dispersion_all['numdis']
    if init:
        outputs={}
        outputs['stack']={}
        outputs['nostack']={}

        outputs['nostack']['clusters']=params_dispersion['clusters']

        outputs['nostack']['ndata_R']=params_dispersion['R']['ndatad']
        outputs['nostack']['ndata_L']=params_dispersion['L']['ndatad']

        outputs['meansum_sq']=[]
        outputs['meansum_sq_shift']=[]

        ndataR=params_dispersion['R']['ndatad']
        ndataL=params_dispersion['L']['ndatad']

        outputs['stack']['dispersion_R']=np.zeros((ndataR,numdis+1))
        outputs['stack']['dispersion_R_sq']=np.zeros((ndataR,numdis+1))
        outputs['stack']['dispersion_R_shift']=np.zeros((ndataR,numdis+1))
        outputs['meansum_sq'].append('dispersion_R_sq')
        outputs['meansum_sq_shift'].append('dispersion_R_shift')

        outputs['stack']['dispersion_L']=np.zeros((ndataL,numdis+1))
        outputs['stack']['dispersion_L_sq']=np.zeros((ndataL,numdis+1))
        outputs['stack']['dispersion_L_shift']=np.zeros((ndataL,numdis+1))
        outputs['meansum_sq'].append('dispersion_L_sq')
        outputs['meansum_sq_shift'].append('dispersion_L_shift')

        return outputs

    dispersion_R_model=dispersion_one['R']['dispersion']
    dispersion_L_model=dispersion_one['L']['dispersion']
    #if np.amax(alpha)<-50:
    #    return
    alpha_max=dispersion_all['alpha_max']
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    alpha=np.exp(alpha)
    dispersion_R_true_all=dispersion_all['R']['dispersion']
    dispersion_L_true_all=dispersion_all['L']['dispersion']

    if numdis==0:
        return

    alpha_alt=np.tile(alpha,(params_dispersion['R']['ndatad'],1))
    outputs['stack']['dispersion_R'][:,:numdis]+=np.multiply(dispersion_R_true_all,alpha_alt)

    dis_alt=np.transpose(np.tile(dispersion_R_model,(numdis,1)))-dispersion_R_true_all
    outputs['stack']['dispersion_R_shift'][:,:numdis]+=np.multiply(dis_alt,alpha_alt)

    dis_alt=np.square(np.transpose(np.tile(dispersion_R_model,(numdis,1)))-dispersion_R_true_all)
    outputs['stack']['dispersion_R_sq'][:,:numdis]+=np.multiply(dis_alt,alpha_alt)

    alpha_alt=np.tile(alpha,(params_dispersion['L']['ndatad'],1))
    outputs['stack']['dispersion_L'][:,:numdis]+=np.multiply(dispersion_L_true_all,alpha_alt)

    dis_alt=np.transpose(np.tile(dispersion_L_model,(numdis,1)))-dispersion_L_true_all
    outputs['stack']['dispersion_L_shift'][:,:numdis]+=np.multiply(dis_alt,alpha_alt)

    dis_alt=np.square(np.transpose(np.tile(dispersion_L_model,(numdis,1)))-dispersion_L_true_all)
    outputs['stack']['dispersion_L_sq'][:,:numdis]+=np.multiply(dis_alt,alpha_alt)

    return #outputs

def get_tradeoff(file,model,model_ref,dispersion_one,params_dispersion,params_inversion,dispersion_all,alpha,outputs={},init=False):

    numdis=dispersion_all['numdis']
    if init:
        outputs={}
        outputs['stack']={}
        outputs['nostack']={}



        malay=params_inversion['malay']

        outputs['nostack']['range']=np.arange(malay)

        outputs['stack']['tradeoff']=np.zeros((malay+1,malay+1,numdis+1))

        return outputs
#    if np.amax(alpha)<-50:
#        nlay=model['npt']
#        nlay_ani=model['npt_ani']
#        outputs['stack']['tradeoff'][nlay,nlay_ani,numdis]+=1.
#        return

    alpha_max=dispersion_all['alpha_max']
    alpha=np.minimum(np.zeros_like(alpha),alpha-alpha_max)
    alpha=np.exp(alpha)

    malay=params_inversion['malay']

    nlay=model['npt']
    nlay_ani=model['npt_ani']
    outputs['stack']['tradeoff'][nlay,nlay_ani,:numdis]+=alpha
    outputs['stack']['tradeoff'][nlay,nlay_ani,numdis]+=1.
    return #outputs

def kl_weights(kl_dist):
    return np.exp(kl_dist)
