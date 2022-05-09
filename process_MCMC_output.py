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
# import matplotlib.pyplot as plt
import concurrent.futures
import sys

import json

# Make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str


# input_directory='OUT_TEST'

num_args=len(sys.argv)
if num_args != 4:
    print('python process_MCMC_output.py <INPUT_DIR> <OUTPUT_DIR> <NUM_CORES>')
    print('Number of arguments (' + str(num_args) +') too low... exit')
    exit('exiting....')
    
input_directory = str(sys.argv[1])
output_directory = str(sys.argv[2])
cores = int(sys.argv[3])
print('Processing results for: ' +str(input_directory)+' using '+str(cores)+' cores...')
print('Sending output to: ' +str(output_directory))



#########################################################
# Get MetaData from one file
#########################################################


def get_metadata(input_directory):
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
    input_directory : str

    Returns
    -------
    params_inversion : dict
        dict containing all the metadata of the inversion.

    '''
    
    params_inversion={}

    files=glob.glob(input_directory+'/All_models_invert*.out')
    
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
    data=f.readline().split() # contains a third parameter, the thinning. Not used currently, but can be added easily
    print(data)
    params_inversion['burn-in']=float(data[0])
    params_inversion['nsample']=float(data[1])
    params_inversion['thinning']=float(data[2])
    params_inversion['cores']=len(files)
    data=f.readline().split()
    params_inversion['d_min']=float(data[0])
    params_inversion['d_max']=float(data[1])
    data=f.readline().split()
    params_inversion['vs_min']=float(data[0])
    params_inversion['vs_max']=float(data[1])
    params_inversion['width_vsv']=params_inversion['vs_max']
    params_inversion['vsref_min']=float(data[2])
    params_inversion['vsref_max']=float(data[3])
    # params_inversion['width_vsv']=float(f.readline())
    data=f.readline().split()
    params_inversion['xi_min']=float(data[0])
    params_inversion['xi_max']=float(data[1])
    data=f.readline().split()
    params_inversion['vp_min']=float(data[0])
    params_inversion['vp_max']=float(data[1])
    params_inversion['vpref_min']=float(data[2])
    params_inversion['vpref_max']=float(data[3])
    data=f.readline().split()
    params_inversion['Ad_R_min']=float(data[0])
    params_inversion['Ad_R_max']=float(data[1])
    data=f.readline().split()
    params_inversion['Ad_L_min']=float(data[0])
    params_inversion['Ad_L_max']=float(data[1])
    f.close()
    
    print(params_inversion['burn-in'],params_inversion['nsample'],params_inversion['thinning'],params_inversion['cores'])
    print(params_inversion['d_min'],params_inversion['d_max'])
    print(params_inversion['vs_min'],params_inversion['vs_max'],params_inversion['vsref_min'],params_inversion['vsref_max'])
    print(params_inversion['xi_min'],params_inversion['xi_max'])
    print(params_inversion['vp_min'],params_inversion['vp_max'],params_inversion['vpref_min'],params_inversion['vpref_max'])
    print(params_inversion['Ad_R_min'],params_inversion['Ad_R_max'])
    print(params_inversion['Ad_L_min'],params_inversion['Ad_L_max'])


    return params_inversion


#########################################################
# Get Reference model values
#########################################################

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

#########################################################
# A function to return sensible bounds on velocities tied to reference model and prior with 5% buffer.
#########################################################

def get_layer_vel_bounds(model_ref, ud, ld, vs_min, vs_max, vp_min, vp_max):
    '''
    Takes reference model and prior on vsv and vph

    Parameters
    ----------
    model_ref : dict
    ud : (float) upper layer depth
    ld : (float) lower layer depth
    vs_min : (float)    Prior on min(vsv) in voro
    vs_max : (float)    Prior on max(vsv) in voro
    vp_min : (float)    Prior on min(vph) in voro
    vp_max : (float)    Prior on max(vph) in voro

    Returns:
    ---------
    lay_vs_min, lay_vsv_max, lay_vph_min, lay_vph_max
    '''
    if not ld > ud:
        print('Lower depth not > upper depth: '+str(ld)+' !> '+str(ud))
        exit('exiting.....')
    # Do upper depth bound
    if len(np.where(model_ref['depth']==ud)[0])>0:
        ud_ind=np.where(model_ref['depth']==ud)[0][0]
        lay_vsv_min=np.round_(model_ref['vsv'][ud_ind]*(1+(vs_min-0.05)),4)
        lay_vph_min=np.round_(model_ref['vph'][ud_ind]*(1+(vp_min-0.05)),4)
    else:
        ud_ind_1=np.where(model_ref['depth']<ud)[0][-1]
        ud_ind_2=np.where(model_ref['depth']>ud)[0][0]
        vsv=np.interp(ud,[model_ref['depth'][ud_ind_1],model_ref['depth'][ud_ind_2]],[model_ref['vsv'][ud_ind_1],model_ref['vsv'][ud_ind_2]])
        vph=np.interp(ud,[model_ref['depth'][ud_ind_1],model_ref['depth'][ud_ind_2]],[model_ref['vph'][ud_ind_1],model_ref['vph'][ud_ind_2]])
        lay_vsv_min=np.round_(vsv*(1+(vs_min)),4)
        lay_vph_min=np.round_(vph*(1+(vp_min-0.05)),4)

    # Do lower depth bound
    if len(np.where(model_ref['depth']==ld)[0])>0:
        ld_ind=np.where(model_ref['depth']==ld)[0][-1]
        lay_vsv_max=np.round_(model_ref['vsv'][ld_ind]*(1+(vs_max+0.05)),4)
        lay_vph_max=np.round_(model_ref['vph'][ld_ind]*(1+(vp_max+0.05)),4)
    else:
        ld_ind_1=np.where(model_ref['depth']<ld)[0][-1]
        ld_ind_2=np.where(model_ref['depth']>ld)[0][0]
        vsv=np.interp(ld,[model_ref['depth'][ld_ind_1],model_ref['depth'][ld_ind_2]],[model_ref['vsv'][ld_ind_1],model_ref['vsv'][ld_ind_2]])
        vph=np.interp(ld,[model_ref['depth'][ld_ind_1],model_ref['depth'][ld_ind_2]],[model_ref['vph'][ld_ind_1],model_ref['vph'][ld_ind_2]])
        lay_vsv_max=np.round_(vsv*(1+(vs_max+0.05)),4)
        lay_vph_max=np.round_(vph*(1+(vp_max+0.05)),4)

    return lay_vsv_min, lay_vsv_max, lay_vph_min, lay_vph_max


#########################################################
# Process the files in parallel.
#########################################################

def apply_stuff(input_directory,cores,functions,params_inversion,model_ref):
    '''
    Takes a list of functions, reads all models, applies each function to all of the models and stacks the results

    Parameters
    ----------
    input_directory: str
        input_directory name where all the data is.
    functions : list
        list of functions.
    params_inversion : dict
        
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
        like_prop: log_likelihood of the model
    
    params_inversion: dict
    
    first: bool
    whether it is the first model to be processes on the file. May provide a minor speed boost. 
    
    The function must give the following output: 
    output: dict
        contains the output of the function. Must have 2 keys, 'stack' and 'nostack' that each have sub-dicts
        values in output['stack'] will be added for key and each model
        values in output['nostack'] will be kept as in the first model
    '''
    
    files=glob.glob(input_directory+'/All_models_invert*.out')

    files.sort()
    
    l=len(files)
    
    outputs_all={}
    for function in functions:
        outputs_all[function.__name__]={}
        outputs_all[function.__name__]['stack']={}
        outputs_all[function.__name__]['nostack']={}

    #parallel processing needs a list of single inputs, so we put all input into a dict and create a list of dicts
    input_dicts=[]
    for file in files:
        input_dict={}
        input_dict['file']=file
        input_dict['functions']=functions
        input_dict['params_inversion']=params_inversion
        input_dict['model_ref']=model_ref
        input_dicts.append(input_dict)
            
    # Parallel processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
        
        results = executor.map(process_one_file, input_dicts)
        
        # collapsing of results
        for res in results:
            for function in res:
                if len(outputs_all[function]['stack'].keys()) > 0:
                    for key in res[function]['stack']: # stack are stacked
                        if type(res[function]['stack'][key])==dict:
                            #Then do something else
                            # print('Here #1 ...')
                            for key2 in res[function]['stack'][key]:
                                for key3 in res[function]['stack'][key][key2]:
                                    # Sum of all values in 2D histogram is 1.
                                    outputs_all[function]['stack'][str(key)][str(key2)][str(key3)]+=res[function]['stack'][str(key)][str(key2)][str(key3)]/l
                        else:
                            # Follow original structure for posterior stacking.
                            outputs_all[function]['stack'][key]+=res[function]['stack'][key] # /l

                else:
                    for key in res[function]['stack']:
                        if type(res[function]['stack'][key])==dict:
                            #Then do something else
                            # print('Here #2 ...')
                            outputs_all[function]['stack'][str(key)]={}
                            for key2 in res[function]['stack'][key]:
                                outputs_all[function]['stack'][str(key)][key2]={}
                                for key3 in res[function]['stack'][key][key2]:
                                    # Sum of all values in 2D histogram is 1.
                                    outputs_all[function]['stack'][str(key)][str(key2)][str(key3)]=res[function]['stack'][str(key)][str(key2)][str(key3)]/l

                        else:
                            # Follow original structure for posterior stacking.
                            outputs_all[function]['stack'][key]=res[function]['stack'][key] # /l

                    for key in res[function]['nostack']:
                        if type(res[function]['nostack'][key])==dict:
                            #Then do something else
                            # print('Here #3 ...')
                            outputs_all[function]['nostack'][str(key)]={}
                            for key2 in res[function]['nostack'][key]:
                                # print(type(key), type(key2), key, key2, res[function]['nostack'].keys())
                                # Need to create the dictionarys in outputs to receive these.
                                outputs_all[function]['nostack'][str(key)][str(key2)]=res[function]['nostack'][str(key)][str(key2)]

                        else:
                            # Follow original structure for posterior stacking.
                            outputs_all[function]['nostack'][key]=res[function]['nostack'][key] # nostack are kept as in the first one
                    
    return outputs_all
    
#########################################################
# Function to process one file.
#########################################################

def process_one_file(input_dict):
    '''
    Processes one file, reading the podels in the file, applying the functions to them and stacking the results

    Parameters
    ----------
    input_dict : dict
        input dictionary containing the file name (file), the list of functions (functions), 
        params_inversion and model_ref.

    Returns
    -------
    outputs : dict
        has one subdict for each function, called with the function name.
        Each of those has a subdict 'stack' and a subdict 'nostack' which will respectively stacked and kept as in the first model

    '''
    
    file=input_dict['file']
    functions=input_dict['functions']
    params_inversion=input_dict['params_inversion']
    model_ref=input_dict['model_ref']
    outputs={}
        
    numtot=0

    # read file
    f=open(file,'r')
    print(file)

    fline=f.readline()
    if not fline:
        print('file ' +file + ' empty')
        # return

    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    line=f.readline()
    while line:

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
        vph=np.zeros_like(d)

        for i in range(npt_true):
            data=f.readline().split()
            d[i]=float(data[0])
            vsv[i]=float(data[1])
            xi[i]=float(data[2])
            vph[i]=float(data[3])
        model['depth']=d
        model['vsv']=vsv
        model['xi']=xi
        model['vph']=vph

        dispersion_one['like_prop']=float(f.readline())
        model['like_prop']=dispersion_one['like_prop']
        # dispersion_one['widening']=widening

        f.readline()
        dispersion_R_one=np.fromiter(map(float,f.readline().split()),float)
        dispersion_one['R']['dispersion']=dispersion_R_one

        f.readline()
        dispersion_L_one=np.fromiter(map(float,f.readline().split()),float)
        dispersion_one['L']['dispersion']=dispersion_L_one

        # apply functions
        for function in functions:

            output=function(model,model_ref,params_inversion,first=(numtot==0))
            
            # stack outputs
            if function.__name__ in outputs:
                for key in output['stack']:
                    if type(output['stack'][key])==dict:
                        #Then do something else
                        # print('Here #4 ...')
                        for key2 in output['stack'][key]:
                            for key3 in output['stack'][key][key2]:
                                outputs[function.__name__]['stack'][str(key)][str(key2)][str(key3)]+=output['stack'][str(key)][str(key2)][str(key3)]

                    else:
                        # Follow original structure for posterior stacking.
                        outputs[function.__name__]['stack'][key]+=output['stack'][key]
            else:
                outputs[function.__name__]={}
                outputs[function.__name__]['stack']={}
                outputs[function.__name__]['nostack']={}
                for key in output['stack']:
                    if type(output['stack'][key])==dict:
                        #Then do something else
                        # print('Here #5 ...')
                        outputs[function.__name__]['stack'][str(key)]={}
                        for key2 in output['stack'][key]:
                            outputs[function.__name__]['stack'][str(key)][str(key2)]={}
                            for key3 in output['stack'][key][key2]:
                                # Need to create the dictionarys in outputs to receive these.
                                outputs[function.__name__]['stack'][str(key)][str(key2)][str(key3)]=output['stack'][str(key)][str(key2)][str(key3)]

                    else:
                        # Follow original structure for posterior stacking.
                        outputs[function.__name__]['stack'][key]=output['stack'][key]

                for key in output['nostack']:
                    if type(output['nostack'][key])==dict:
                        #Then do something else
                        # print('Here #6 ...')
                        outputs[function.__name__]['nostack'][str(key)]={}
                        for key2 in output['nostack'][key]:
                            # Need to create the dictionarys in outputs to receive these.
                            outputs[function.__name__]['nostack'][str(key)][str(key2)]=output['nostack'][str(key)][str(key2)]
                    else:
                        # Follow original structure for posterior stacking.
                        outputs[function.__name__]['nostack'][key]=output['nostack'][key]

        line=f.readline()
        numtot+=1

    # normalize outputs
    for function in functions:
        for key in outputs[function.__name__]['stack']:
            if type(outputs[function.__name__]['stack'][key])==dict:
                #Then do something else
                # print('Here #7 ...')
                for key2 in outputs[function.__name__]['stack'][key]:
                    for key3 in outputs[function.__name__]['stack'][key][key2]:
                        # Sum of all values in 2D histogram is 1.
                        outputs[function.__name__]['stack'][str(key)][str(key2)][str(key3)]/=numtot

            else:
                # Follow original structure for posterior normalising
                # outputs[function.__name__]['stack'][key] # /=numtot
                dummy=1
    f.close()
    return outputs

#########################################################
# Example to produce posterior for vsv.
#########################################################

def create_post_array_ref(model,model_ref,params_inversion,first=True):
    '''
    example of a function to be applied to the data
    creates a 2D histogram of vsv, xi and vph for the reference
    The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

    Parameters
    ----------
    model : dict

    model_ref : dict

    params_inversion : dict

    first : bool, optional
        whether this model is the first of its file. May provide a minor speed boost. The default is True.

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    ndatad=200
    ndatav=100
    outputs={}
    outputs['stack']={}
    if first:
        outputs['nostack']={}
    vsv_ref=np.zeros((ndatad,ndatav))
    xi_ref=np.zeros((ndatad,ndatav))
    vph_ref=np.zeros((ndatad,ndatav))

    depths=np.linspace(params_inversion['d_min'],params_inversion['d_max'],ndatad)

    vsv_vels=np.linspace(np.amin(params_inversion['vsref_min']),np.amax(params_inversion['vsref_max']),ndatav)

    xi_vels=np.linspace(np.amin(params_inversion['xi_min']),np.amax(params_inversion['xi_max']),ndatav)

    vph_vels=np.linspace(np.amin(params_inversion['vpref_min']),np.amax(params_inversion['vpref_max']),ndatav)

    if first:
        outputs['nostack']['depths']=depths
        outputs['nostack']['vsv_vels']=vsv_vels
        outputs['nostack']['xi_vels']=xi_vels
        outputs['nostack']['vph_vels']=vph_vels
        outputs['nostack']['ndatad']=ndatad
        outputs['nostack']['ndatav']=ndatav
    
    vsv_model=model['vsv']
    xi_model=model['xi']
    vph_model=model['vph']
    depth_model=model['depth']

    ind_vsv=np.digitize(np.interp(depths,depth_model,vsv_model),bins=vsv_vels,right=True)
    vsv_ref[np.arange(ndatad),ind_vsv]+= 1 

    ind_xi=np.digitize(np.interp(depths,depth_model,xi_model),bins=xi_vels,right=True)
    xi_ref[np.arange(ndatad),ind_xi]+= 1 
    
    ind_vph=np.digitize(np.interp(depths,depth_model,vph_model),bins=vph_vels,right=True)
    vph_ref[np.arange(ndatad),ind_vph]+= 1 

    outputs['stack']['vsv_ref']=vsv_ref
    outputs['stack']['xi_ref'] =xi_ref
    outputs['stack']['vph_ref']=vph_ref

    return outputs

def write_posterior_vsv_xi_vph(output_directory,input,params_inversion):
    '''
    Writes posterior for vsv, xi and vph to mimic the original fortran output.
    
    write(*,*)prof,disd,d_max
    write(*,*)vsref_min,vsref_max,disv,width,xi_min,xi_max,vpref_min,vpref_max
    do i=1,disd
        do j=1,disv
            write(71,*)postvss(i,j),postxis(i,j),postvps(i,j)
        enddo
    enddo

    Parameters
    ----------
    output_directory : str - location to save output

    input : dictionary of posterior to be saved: vsv, xi, vph

    params_inversion : dict

    Returns
    -------
    outputs : saved output file: filename_out

    '''
    filename_out=output_directory+'/Proc_Posterior.out'

    if not 'vsv_ref' in input['stack']:
        sys.exit('Input dictionary is incomplete (vsv)- cannot be saved')
    if not 'xi_ref' in input['stack']:
        sys.exit('Input dictionary is incomplete (xi) - cannot be saved')
    if not 'vph_ref' in input['stack']:
        sys.exit('Input dictionary is incomplete (vph)- cannot be saved')

    out_file = open(filename_out, 'w')

    # Write Depth metadata
    # print(int(params_inversion['d_max']),int(input['nostack']['ndatad']),int(params_inversion['d_max']))
    out_file.write(
        "%i %i %i %i %i %i %i\n" %
        (int(params_inversion['d_max']),int(input['nostack']['ndatad']),int(params_inversion['d_max']),int(params_inversion['burn-in']),int(params_inversion['nsample']),int(params_inversion['thinning']),int(params_inversion['cores'])))
    # Write Velocity metadata
    # print(params_inversion['vsref_min'],params_inversion['vsref_max'],int(input['nostack']['ndatav']),params_inversion['width_vsv'],params_inversion['xi_min'],params_inversion['xi_max'],params_inversion['vpref_min'],params_inversion['vpref_max'])
    out_file.write(
    "%f %f %i %f %f %f %f %f\n" %
    (params_inversion['vsref_min'],params_inversion['vsref_max'],int(input['nostack']['ndatav']),params_inversion['width_vsv'],params_inversion['xi_min'],params_inversion['xi_max'],params_inversion['vpref_min'],params_inversion['vpref_max']))
    
    # Write posterior
    for i in range(input['nostack']['ndatad']):
        for j in range(input['nostack']['ndatav']):
            # print(input['stack']['vsv_ref'][i,j],input['stack']['xi_ref'][i,j],input['stack']['vph_ref'][i,j])
            out_file.write(
                "%f %f %f\n" %
                (input['stack']['vsv_ref'][i,j],input['stack']['xi_ref'][i,j],input['stack']['vph_ref'][i,j]))

    out_file.close()

    print('Written Processed Posterior to: '+str(filename_out))
    return()

#########################################################
# Produce correlations between posterior vsv, xi, vph
#########################################################

def posterior_correlations(model,model_ref,params_inversion,first=True):
    '''
    Compute the discretized correlations between parameters averaged over depth intervals 
    creates a 2D histogram of correlation between each parameter at each depth
    The functions will be called a lot of times, so vectorisation is very important, optimize as much as possible

    Parameters
    ----------
    model : dict

    model_ref : dict

    params_inversion : dict

    first : bool, optional
        whether this model is the first of its file. May provide a minor speed boost. The default is True.

    Returns
    -------
    outputs : dict
        has 2 subdicts, 'stack' and 'nostack'.

    '''
    phi = 1 # P-wave radial anisotropy variable: phi= (vpv/vph)^2
    av_ints=[100,200]
    dhs    =[50,100,200]
    outputs={}
    outputs['stack']={}
    if first:
        outputs['nostack']={}

    # Read model into useable format

    a_mod=np.zeros([model['npt_true'],7])
    a_mod[:,0]=model['depth']
    a_mod[:,2]=model['vsv']
    a_mod[:,3]=model['vsv']*np.sqrt(model['xi']) # Make VSH
    a_mod[:,4]=model['vph']*np.sqrt(phi) # Make VPV
    a_mod[:,5]=model['vph']
    a_mod[:,6]=model['xi']

    # Add layer thickness.
    for i in range(model['npt_true']-1):
        a_mod[i,1]=a_mod[i+1,0]-a_mod[i,0]
    # If gradients we will need to do more here... None for now.

    for av_ind, av_int in enumerate(av_ints):

        dep_array=np.arange(params_inversion['d_min'],params_inversion['d_max']+0.1,av_int)

        # Make depth arrays for averaging that 
        av_u_dep=[]
        av_l_dep=[]
        av_vsv_min=[]
        av_vsv_max=[]
        av_vph_min=[]
        av_vph_max=[]

        for i in range(len(dep_array)-1):
            av_u_dep.append(dep_array[i])
            av_l_dep.append(dep_array[i+1])
            lay_vsv_min, lay_vsv_max, lay_vph_min, lay_vph_max=get_layer_vel_bounds(model_ref,dep_array[i],dep_array[i+1],-1.0*params_inversion['width_vsv'],params_inversion['width_vsv'],params_inversion['vp_min'],params_inversion['vp_max'])
            av_vsv_min.append(lay_vsv_min)
            av_vsv_max.append(lay_vsv_max)
            av_vph_min.append(lay_vph_min)
            av_vph_max.append(lay_vph_max)

        av_vsv_min=np.array(av_vsv_min)
        av_vsv_max=np.array(av_vsv_max)
        av_vph_min=np.array(av_vph_min)
        av_vph_max=np.array(av_vph_max)    
        
        outputs['stack'][str(av_int)]={}
        if first:
            outputs['nostack'][str(av_int)]={}
            outputs['nostack'][str(av_int)]['av_u_deps']=av_u_dep
            outputs['nostack'][str(av_int)]['av_l_deps']=av_l_dep
            outputs['nostack'][str(av_int)]['dhs']=dhs
            outputs['nostack'][str(av_int)]['av_vsv_min']=av_vsv_min
            outputs['nostack'][str(av_int)]['av_vsv_max']=av_vsv_max
            outputs['nostack'][str(av_int)]['av_vph_min']=av_vph_min
            outputs['nostack'][str(av_int)]['av_vph_max']=av_vph_max
        
        av_u_dep=np.array(av_u_dep)
        av_l_dep=np.array(av_l_dep)

        av_vsv=np.zeros(len(av_u_dep))
        # av_vsh=np.zeros(len(av_u_dep))
        # av_vpv=np.zeros(len(av_u_dep))
        av_vph=np.zeros(len(av_u_dep))
        av_xi =np.zeros(len(av_u_dep))

        # Loop over model layers
        for j in range(model['npt_true']-1):
            # Check we are not at end or have thickness = 0
            if a_mod[j,0] < params_inversion['d_max'] and a_mod[j,1] > 0:
                vsv1=a_mod[j,2]
                # vsh1=a_mod[j,3]
                # vpv1=a_mod[j,4]
                vph1=a_mod[j,5]
                xi1 =a_mod[j,6]
                u_dep=a_mod[j,0]   # Depth of upper interface of layer 
                l_dep=a_mod[j+1,0] # Depth of lower interface of layer
                l_th=a_mod[j,1]    # Layer thickness
                # print(u_dep, l_dep, l_th, vsv1, vsh1, vpv1, vph1, xi1)

                for k in range(len(av_u_dep)):
                    if u_dep < av_u_dep[k] and l_dep <= av_u_dep[k]:                     # CASE 1
                        ### DO NOTHING HERE
                        # av_vsv[k]=av_vsv[k]+0.0
                        dummy=1
                    elif u_dep >= av_l_dep[k] and l_dep > av_l_dep[k]:                 # CASE 2
                        ### DO NOTHING HERE
                        # av_vsv[k]=av_vsv[k]+0.0
                        dummy=1
                    elif u_dep < av_u_dep[k] and l_dep < av_l_dep[k]:                 # CASE 3.1
                        av_vsv[k]=av_vsv[k]+vsv1*(l_dep-av_u_dep[k])
                        # av_vsh[k]=av_vsh[k]+vsh1*(l_dep-av_u_dep[k])
                        # av_vpv[k]=av_vpv[k]+vpv1*(l_dep-av_u_dep[k])
                        av_vph[k]=av_vph[k]+vph1*(l_dep-av_u_dep[k])
                        av_xi[k]=av_xi[k]+xi1*(l_dep-av_u_dep[k])

                    elif u_dep < av_u_dep[k] and l_dep == av_l_dep[k]:                 # CASE 3.2
                        av_vsv[k]=av_vsv[k]+vsv1*av_int
                        # av_vsh[k]=av_vsh[k]+vsh1*av_int
                        # av_vpv[k]=av_vpv[k]+vpv1*av_int
                        av_vph[k]=av_vph[k]+vph1*av_int
                        av_xi[k]=av_xi[k]+xi1*av_int

                    elif u_dep < av_u_dep[k] and l_dep > av_l_dep[k]:                 # CASE 3.3
                        av_vsv[k]=av_vsv[k]+vsv1*av_int
                        # av_vsh[k]=av_vsh[k]+vsh1*av_int
                        # av_vpv[k]=av_vpv[k]+vpv1*av_int
                        av_vph[k]=av_vph[k]+vph1*av_int
                        av_xi[k]=av_xi[k]+xi1*av_int

                    elif u_dep == av_u_dep[k] and l_dep < av_l_dep[k]:                 # CASE 4.1
                        av_vsv[k]=av_vsv[k]+vsv1*l_th
                        # av_vsh[k]=av_vsh[k]+vsh1*l_th
                        # av_vpv[k]=av_vpv[k]+vpv1*l_th
                        av_vph[k]=av_vph[k]+vph1*l_th
                        av_xi[k]=av_xi[k]+xi1*l_th

                    elif u_dep == av_u_dep[k] and l_dep == av_l_dep[k]:                 # CASE 4.2
                        av_vsv[k]=av_vsv[k]+vsv1*l_th
                        # av_vsh[k]=av_vsh[k]+vsh1*l_th
                        # av_vpv[k]=av_vpv[k]+vpv1*l_th
                        av_vph[k]=av_vph[k]+vph1*l_th
                        av_xi[k]=av_xi[k]+xi1*l_th

                    elif u_dep == av_u_dep[k] and l_dep > av_l_dep[k]:                 # CASE 4.3
                        av_vsv[k]=av_vsv[k]+vsv1*av_int
                        # av_vsh[k]=av_vsh[k]+vsh1*av_int
                        # av_vpv[k]=av_vpv[k]+vpv1*av_int
                        av_vph[k]=av_vph[k]+vph1*av_int
                        av_xi[k]=av_xi[k]+xi1*av_int

                    elif u_dep > av_u_dep[k] and l_dep < av_l_dep[k]:                 # CASE 5.1
                        av_vsv[k]=av_vsv[k]+vsv1*l_th
                        # av_vsh[k]=av_vsh[k]+vsh1*l_th
                        # av_vpv[k]=av_vpv[k]+vpv1*l_th
                        av_vph[k]=av_vph[k]+vph1*l_th
                        av_xi[k]=av_xi[k]+xi1*l_th

                    elif u_dep > av_u_dep[k] and l_dep == av_l_dep[k]:                 # CASE 5.2
                        av_vsv[k]=av_vsv[k]+vsv1*l_th
                        # av_vsh[k]=av_vsh[k]+vsh1*l_th
                        # av_vpv[k]=av_vpv[k]+vpv1*l_th
                        av_vph[k]=av_vph[k]+vph1*l_th
                        av_xi[k]=av_xi[k]+xi1*l_th

                    elif u_dep > av_u_dep[k] and l_dep > av_l_dep[k]:                 # CASE 5.3
                        av_vsv[k]=av_vsv[k]+vsv1*(av_l_dep[k]-u_dep)
                        # av_vsh[k]=av_vsh[k]+vsh1*(av_l_dep[k]-u_dep)
                        # av_vpv[k]=av_vpv[k]+vpv1*(av_l_dep[k]-u_dep)
                        av_vph[k]=av_vph[k]+vph1*(av_l_dep[k]-u_dep)
                        av_xi[k]=av_xi[k]+xi1*(av_l_dep[k]-u_dep)
        # Normalise
        av_vsv=av_vsv/av_int
        # av_vsh=av_vsh/av_int
        # av_vpv=av_vpv/av_int
        av_vph=av_vph/av_int
        av_xi=av_xi/av_int



        # Loop over discretization intervals.
        for dh_ind, dh in enumerate(dhs):
            outputs['stack'][str(av_int)][str(dh)]={}
            # print(outputs['stack'][str(av_int)].keys())
            # Loop over model (a_mod) layers


            # Clip to catch strange errors : operate on arrays
            # i_av_vsv=np.clip(np.floor((av_vsv-params_inversion['vsref_min'])*dh/(params_inversion['vsref_max']-params_inversion['vsref_min'])),0,dh-1)
            # # i_av_vsh=np.clip(np.floor((av_vsh-params_inversion['vsref_min'])*dh/(params_inversion['vsref_max']-params_inversion['vsref_min'])),0,dh-1)
            # # i_av_vpv=np.clip(np.floor((av_vpv-params_inversion['vpref_min'])*dh/(params_inversion['vpref_max']-params_inversion['vpref_min'])),0,dh-1)
            # i_av_vph=np.clip(np.floor((av_vph-params_inversion['vpref_min'])*dh/(params_inversion['vpref_max']-params_inversion['vpref_min'])),0,dh-1)
            # i_av_xi= np.clip(np.floor((av_xi-params_inversion['xi_min'])*    dh/(params_inversion['xi_max']   -params_inversion['xi_min'])),0,dh-1)

            # Clip to catch strange errors : operate on arrays
            i_av_vsv=np.clip(np.floor((av_vsv-av_vsv_min)*dh/(av_vsv_max-av_vsv_min)),0,dh-1)
            i_av_vph=np.clip(np.floor((av_vph-av_vph_min)*dh/(av_vph_max-av_vph_min)),0,dh-1)
            i_av_xi= np.clip(np.floor((av_xi-params_inversion['xi_min'])*    dh/(params_inversion['xi_max']   -params_inversion['xi_min'])),0,dh-1)

            # print("_________________ANSWER: _________________")
            # print("av_vsv:     ",av_vsv[:]                    )
            # print("av_vsh:     ",av_vsh[:]                    )
            # print("av_vpv:     ",av_vpv[:]                    )
            # print("av_vph:     ",av_vph[:]                    )
            # print("av_xi:      ",av_xi[:]                     )
            # print("__________________________________________")
            # print("av_vsv:     ",i_av_vsv[:]                  )
            # print("av_vsh:     ",i_av_vsh[:]                  )
            # print("av_vpv:     ",i_av_vpv[:]                  )
            # print("av_vph:     ",i_av_vph[:]                  )
            # print("av_xi:      ",i_av_xi[:]                   )
            # print("__________________________________________")

            param_pairs=[('vph','xi'), ('vph','vsv'), ('vsv','xi')]

            for pp in param_pairs:
                for d1, deps_1 in enumerate(av_u_dep):
                    for d2, deps_2 in enumerate(av_u_dep):
                        deps_1=str(int(deps_1))
                        deps_2=str(int(deps_2))
                        
                        hist_name='h_'+pp[0]+'_'+deps_1+'_'+pp[1]+'_'+deps_2
                        # Need to make arrays size dh x dh 
                        outputs['stack'][str(av_int)][str(dh)][str(hist_name)]=np.zeros([dh,dh])

                        # print(pp, d1, deps_1, d2, deps_2, hist_name)       
                        if pp == ('vph','xi'):
                            outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vph[d1]),int(i_av_xi[d2])] = outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vph[d1]),int(i_av_xi[d2])]+1
                        if pp == ('vph','vsv'):
                            outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vph[d1]),int(i_av_vsv[d2])] = outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vph[d1]),int(i_av_vsv[d2])]+1
                        if pp == ('vsv','xi'):
                            outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vsv[d1]),int(i_av_xi[d2])] = outputs['stack'][str(av_int)][str(dh)][str(hist_name)][int(i_av_vsv[d1]),int(i_av_xi[d2])]+1
            # print(h)


    return outputs

def write_posterior_corr_dict(output_directory,input,params_inversion):
    '''
    Writes posterior correlations between vsv, xi and vph.
    
    Parameters
    ----------
    output_directory : str - location to save output

    input : dictionary of posterior correlations between vsv, xi, vph to be saved

    params_inversion : dict

    Returns
    -------
    outputs : saved output dictionary of correlations: json


    '''

    input['params_inversion']={}

    input['params_inversion']['burn-in'] = params_inversion['burn-in']
    input['params_inversion']['nsample'] = params_inversion['nsample']
    input['params_inversion']['thinning'] = params_inversion['thinning']
    input['params_inversion']['cores'] = params_inversion['cores']
    input['params_inversion']['d_min'] = params_inversion['d_min']
    input['params_inversion']['d_max'] = params_inversion['d_max']
    input['params_inversion']['vsref_min'] = params_inversion['vsref_min']
    input['params_inversion']['vsref_max'] = params_inversion['vsref_max']
    input['params_inversion']['xi_min'] = params_inversion['xi_min']
    input['params_inversion']['xi_max'] = params_inversion['xi_max']
    input['params_inversion']['vpref_min'] = params_inversion['vpref_min']
    input['params_inversion']['vpref_max'] = params_inversion['vpref_max']


    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)

    filename_out=output_directory+'/Posterior_corr_100_200.json'

    with io.open(filename_out, 'w', encoding='utf8') as outfile:
        str_ = json.dumps(input,indent=4, sort_keys=True,
                          separators=(',', ': '), ensure_ascii=False, cls=NumpyEncoder)
        outfile.write(to_unicode(str_))

    print('Written Posterior correlations to: '+str(filename_out))
    return()


############################################################################
#        main()
############################################################################


def main():
    print('start')
    model_ref=get_model_ref()
    print('got model ref')
    params_inversion=get_metadata(input_directory)
    print('got metadata')
    output=apply_stuff(input_directory,cores,[create_post_array_ref, posterior_correlations],params_inversion,model_ref)
    print('applied functions')

    # udep=200
    # ldep=300
    # lay_vsv_min, lay_vsv_max, lay_vph_min, lay_vph_max=get_layer_vel_bounds(model_ref,udep,ldep,-1.0*params_inversion['width_vsv'],params_inversion['width_vsv'],params_inversion['vp_min'],params_inversion['vp_max'])

    write_posterior_vsv_xi_vph(output_directory,input=output['create_post_array_ref'],params_inversion=params_inversion)

    write_posterior_corr_dict(output_directory,input=output['posterior_correlations'],params_inversion=params_inversion)




if __name__ == '__main__':
    main()
    