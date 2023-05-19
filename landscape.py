#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 15:38:21 2023

@author: magdalena
"""

#import random
from itertools import product
#import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
#import Levenshtein
#from termcolor import colored
#import copy
#import warnings
#import yaml
#import os

def alignment_distance(word1, word2):
    #
    alignments = pairwise2.align.globalms(word1, word2, match=1, mismatch=0, open=-1.5, extend=-.1, penalize_end_gaps = (0,0))
    score_alignment = alignments[0].score
    if len(word1)<len(word2):
        lenght = len(word1)
    else:
        lenght = len(word2)
    #
    distance_alignment = lenght - float(score_alignment)
    #
    return distance_alignment

def create_initial_system(bs1,bs2,cs1,cs2,ss1,ss2):
    
    d_system = {}
    d_system['TFs'] = {}
    d_system['TF_signals'] = {}
    d_system['TF_concentration'] = {}
    d_system['promoters'] = {}

    d_system['promoters'][0] = bs1
    d_system['promoters'][1] = bs2
    d_system['TFs'][0] = cs1
    d_system['TFs'][1] = cs2
    d_system['TF_signals'][0] = ss1
    d_system['TF_signals'][1] = ss2
    d_system['TF_concentration'][0] = 25
    d_system['TF_concentration'][1] = 25

    return d_system

def mismatches_function(d_system, sequences):
    Mismatches = np.zeros((len(d_system['TFs']),len(sequences)))#each row is a TF, each column a sequence in sequence space

    for itf,tf in enumerate(d_system['TFs']):
        word1 = d_system['TFs'][tf]
        for isequence,sequence in enumerate(sequences):
            word2 = sequence
            Mismatches[itf,isequence] = alignment_distance(word1, word2)

    return(Mismatches)

def p_binding_function(d_system,sequences,environment,energy,Mismatches):
    #calculate probability that sequences are bound by any of the TFs that are expressed in the environment
    p_binding = [0]*len(sequences)
    for ipromoter in range(len(sequences)):
        numerator = 0
        for itf,tf in enumerate(d_system['TF_signals']):
            is_tf_in_env = [si for si,s in enumerate(environment) if environment[si]=='1' and d_system['TF_signals'][tf][si]=='1']
            if len(is_tf_in_env)>0:
                numerator = numerator + d_system['TF_concentration'][tf]*np.exp(-1*energy*Mismatches[itf,ipromoter])

        p_binding[ipromoter] = numerator/(1+numerator)
        
    return(p_binding)


energy = 3

nucleotides = ['a','c','g','t']
BS_length = 3
sequence_space_lenght3 = [''.join(comb) for comb in product(nucleotides, repeat=BS_length)]


save_folder = '/Users/magdalena/polybox/research/chance_contingency/'
save_file_name = 'landscape.txt'
# with open(save_folder+save_file_name, 'w') as file:

#     count = 0
#     for ss1 in ['00','01','10','11']:
#         for cs1 in sequence_space_lenght3:
#             for ss2 in ['00','01','10','11']:
#                 for cs2 in sequence_space_lenght3:
#                     for bs1 in sequence_space_lenght3:
#                         for bs2 in sequence_space_lenght3:
#                             count = count+1
#                             file.write(ss1+cs1+ss2+cs2+bs1+bs2)
                            
#                             if ss1=='00' and ss2=='00':
#                                 file.write(',[0,0],[0,0],[0,0]')
#                             else:
#                                 d_system = create_initial_system(bs1,bs2,cs1,cs2,ss1,ss2)
                                
#                                 sequences = [d_system['promoters'][p] for p in d_system['promoters']]
#                                 Mismatches = mismatches_function(d_system,sequences)
#                                 for environment in ['01','10','11']:
#                                     p_binding = p_binding_function(d_system,sequences,environment,energy,Mismatches)
#                                     file.write(','+str(p_binding))
                                
                            
with open(save_folder+save_file_name, 'w') as file:

    count = 0
    for ss1 in ['00','01','10','11']:
        for cs1 in sequence_space_lenght3:
            for ss2 in ['00','01','10','11']:
                for cs2 in sequence_space_lenght3:
                    for bs1 in sequence_space_lenght3:
                        for bs2 in sequence_space_lenght3:
                            count = count+1
                            
                            
                            if ss1=='00' and ss2=='00':
                                continue
                            else:
                                d_system = create_initial_system(bs1,bs2,cs1,cs2,ss1,ss2)
                                
                                sequences = [d_system['promoters'][p] for p in d_system['promoters']]
                                Mismatches = mismatches_function(d_system,sequences)
                                for environment in ['01','10','11']:
                                    p_binding = p_binding_function(d_system,sequences,environment,energy,Mismatches)
                                    if p_binding[0]==0 and p_binding[1]==0:
                                        continue
                                    else:
                                        file.write(ss1+cs1+ss2+cs2+bs1+bs2)
                                        file.write(','+str(round(p_binding),4))                               