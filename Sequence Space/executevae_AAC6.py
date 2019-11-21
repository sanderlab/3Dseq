#!/usr/bin/env python

import os, sys, copy, time

import numpy as np
import pandas as pd

import model
import helper
import train

import theano
#theano.config.gcc.cxxflags = "-Wno-c++11-narrowing" #FIXES ISSUE WITH COMPILING THEANO LIVE COMPILATION

ALIGNMENT_DIR_AAC6 = './data.AAC6'

AAC6_ATTRIBUTES = {
    'alignment_directory': ALIGNMENT_DIR_AAC6,
    'alignment': ALIGNMENT_DIR_AAC6+'/AAC6_natural_>85%ungapped_5K_subsample.a2m',
    'vae_file_prefix': 'AAC6_2latent_5ksubsample',
}

#
#
# data munging functions
#
#
def getVAEDataHelper(ds_attributes):
    return helper.DataHelper(
        alignment_file=ds_attributes['alignment'],
        working_dir=ds_attributes['alignment_directory'],
        calc_weights=False,
    )

def loadVAEModel(ds_attributes):
    start = time.time()
    data_helper = getVAEDataHelper(ds_attributes)
    vae_model   = model.VariationalAutoencoder(
        data_helper,
        batch_size                = 100,
        encoder_architecture      = [1500,1500],
        decoder_architecture      = [100,500],
        n_latent                  = 2,
        n_patterns                = 4,
        warm_up                   = 0.0,
        convolve_patterns         = True,
        conv_decoder_size         = 40,
        logit_p                   = 0.001,
        sparsity                  = 'logit',
        encode_nonlinearity_type  = 'relu',
        decode_nonlinearity_type  = 'relu',
        final_decode_nonlinearity = 'sigmoid',
        output_bias               = True,
        final_pwm_scale           = True,
        working_dir               = data_helper.working_dir
    )
    print('Done loading model - took {0}s'.format(time.time()-start))
    return vae_model

def loadVAEParameters(ds_attributes, vae_model):
    start = time.time()
    vae_model.load_parameters(file_prefix=ds_attributes['vae_file_prefix'])
    print('Done Loading parameters - took {0}s'.format(time.time()-start))
    return vae_model

def trainVAEModel(ds_attributes, vae_model):
    '''
    very slow
    '''
    start = time.time()
    data_helper = getVAEDataHelper(ds_attributes)
    train.train(
        data_helper, vae_model,
        num_updates             =   300000,
        save_progress           =   True,
        verbose                 =   True,
        save_parameters         =   False,
        job_string              =   ds_attributes['vae_file_prefix']
    )
    vae_model.save_parameters(file_prefix=ds_attributes['vae_file_prefix'])
    print('Done training and saving parameters - took {0}s'.format(time.time()-start))


def loadOrCreateVAEModel(attributes):
    '''
    Shorthand to load model and parameters -- or create them if they
    don't exist -- and then return the loaded model.
    '''
    try:
        vae_model = loadVAEModel(attributes)
        vae_model = loadVAEParameters(attributes, vae_model)
    except IOError:
        trainVAEModel(attributes, vae_model)#takes days
        vae_model = loadVAEModel(attributes)
        vae_model = loadVAEParameters(attributes, vae_model)
    return vae_model


loadOrCreateVAEModel(AAC6_ATTRIBUTES)

