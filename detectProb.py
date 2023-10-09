#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:40:56 2020

@author: elisabethburnor
"""

# Function takes the simulated mean concentration of bacteria (simMean) at a 
# any branch and section.
# and three user-inputted fixed parameter values
# sampleVol: the volume to be sampled
# sensitivity: the probability of a positive result, given a certain number
# of bacteria in that sample
# bactNum: the number of bacteria in a sample that relates to the sensitivity
import math as math
import numpy as np
# import time as time

#simMean = 10E-06
#sampleVol = 1
#sensitivity = 0.99
#bactNum = 1
# Old version
def detectProbMC(simMean, sampleVol, sensitivity, bactNum):
   
    # Set beta parameter for sensitivity
    beta = -math.log(1-sensitivity) / bactNum
    
    # Determine mean number of bacteria in the sampled volume of water
    meanCFU = simMean * sampleVol
    
    # Create empty array to store monte carlo estimates
    detectProb = np.zeros(1000)
    
    # Monte Carlo: iterate 100 times to estimate detection probabilities
    for i in range(1000):
        # Estimate the number of bacteria in a sample volume with a poisson distr.
        k = np.random.poisson(lam = meanCFU)
    
        # Calculate the probability of detection in one sample, given k bacteria
        a = -beta * k
        prob = 1 - math.e ** a
        
        # Store the detection probability in an array
        detectProb[i] = prob
    
    
    detectProbMean = np.mean(detectProb)
    
    #print("Detect Prob ", detectProbMean)
    
    detectProbSD = np.std(detectProb)
    detectProbSE = detectProbSD / math.sqrt(1000)
    
    detectProbMedian = np.median(detectProb)
    
    detectUpper = detectProbMean + 1.96*detectProbSE
    detectLower = detectProbMean - 1.96*detectProbSE
    
    return detectProbMean, detectProbSD, detectLower, detectUpper, detectProbMedian

"""
def detectProbMC(simMean, sampleVol, sensitivity, bactNum):
   
    # Set beta parameter for sensitivity
    beta = -math.log(1-sensitivity) / bactNum
    
    
    # Determine mean number of bacteria in the sampled volume of water
    meanCFU = simMean * sampleVol
    
    # Estimate the number of bacteria in a sample volume with a poisson distr.
    k = np.random.poisson(lam = meanCFU)
    
    # Calculate the probability of detection in one sample, given k bacteria
    a = -beta * k
    prob = 1 - math.e ** a
        
    return prob

"""