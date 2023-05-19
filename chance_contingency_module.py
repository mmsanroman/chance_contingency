#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 21:11:14 2023

@author: magdalena
"""

import random
from itertools import product
#import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
#import Levenshtein
#from termcolor import colored
import copy
import warnings
import yaml
import os
import math
    

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

def sequence_space(BS_length):
    nucleotides = ['a','c','g','t']

    #create dictionary with the sequence space for different binding site lenghts
    d_sequence_space = {}
    d_sequence_space[BS_length] = [''.join(comb) for comb in product(nucleotides, repeat=BS_length)]
    
    return d_sequence_space

def create_initial_system(sequences):
    #description of the initial system
    
    d_system = {}
    d_system['TFs'] = {}
    d_system['TF_signals'] = {}
    d_system['TF_concentration'] = {}
    d_system['promoters'] = {}

    #duplicated TFs
    d_system['TFs'][0] = sequences[0]
    d_system['TFs'][1] = sequences[1]
    signals = '11'
    d_system['TF_signals'][0] = signals
    d_system['TF_signals'][1] = signals
    d_system['TF_concentration'][0] = 25
    d_system['TF_concentration'][1] = 25
    
    #genes activated by the duplicated TFs
    d_system['promoters'][0] = sequences[0]
    d_system['promoters'][1] = sequences[1]

    return d_system

def mismatches_function(d_system, sequences):
    Mismatches = np.zeros((len(d_system['TFs']),len(sequences)))#each row is a TF, each column a sequence in sequence space or binding site

    for itf,tf in enumerate(d_system['TFs']):
        word1 = d_system['TFs'][tf]
        for isequence,sequence in enumerate(sequences):
            word2 = sequence
            Mismatches[itf,isequence] = alignment_distance(word1, word2)

    return(Mismatches)

def p_binding_function(d_system,sequences,environment,Mismatches):
    #calculate probability that sequences are bound by any of the TFs that are expressed in the environment
    energy = 3
    p_binding = [0]*len(sequences)
    for ipromoter in range(len(sequences)):
        numerator = 0
        for itf,tf in enumerate(d_system['TF_signals']):
            is_tf_in_env = [si for si,s in enumerate(environment) if environment[si]=='1' and d_system['TF_signals'][tf][si]=='1']
            if len(is_tf_in_env)>0:
                numerator = numerator + d_system['TF_concentration'][tf]*np.exp(-1*energy*Mismatches[itf,ipromoter])

        p_binding[ipromoter] = numerator/(1+numerator)
        
    return(p_binding)

def p_binding_optimal_function():
    p_binding_optimal = [0]*2
    #
    TF_concentration = 25
    energy = 3
    #no binding
    p_binding_optimal[0] = 0
    #perfect binding
    k_tf_promoter = 0
    numerator = TF_concentration*np.exp(-1*energy*k_tf_promoter)
    p_binding_optimal[1] = numerator/(1+numerator)
    
    
    return(p_binding_optimal)
    

def mismatch_signal_geneExpression_function(d_system, environments):
    #calculate fitness (mismatch between observed and optimal gene expression)
    sequences = [d_system['promoters'][p] for p in d_system['promoters']]
    Mismatches = mismatches_function(d_system,sequences)
    #
    p_binding_optimal = p_binding_optimal_function()
    #
    mismatch_signal_genes = 0
    all_p_binding = []
    for environment in environments:
        #calculate probability of sequence (promoter) being bound by any TF
        sequences = d_system['promoters']
        p_binding = p_binding_function(d_system,sequences,environment,Mismatches)
        all_p_binding.append(p_binding) 
        ###
        #mismatch between the expected and the optimal
        genes = len(d_system['promoters'])
        
        #expression = p_binding
        expression = [math.floor(p*10)/10 for p in p_binding]
        # expression = [0]*genes
        # for ipromoter in range(genes):
        #     if p_binding[ipromoter]>2/3:
        #         expression[ipromoter] = 1
        #     elif p_binding[ipromoter]>1/3 and p_binding[ipromoter]<=2/3:
        #         expression[ipromoter] = 0.5
        #     elif p_binding[ipromoter]<=1/3:
        #         expression[ipromoter] = 0
        
        optimal_expression = [math.floor(p_binding_optimal[int(e)]*10)/10 for e in environment]
        mismatch_signal_genes = mismatch_signal_genes + sum([(optimal_expression[i]-expression[i])**2 for i in range(genes)]) 
        
    #
    return (all_p_binding, mismatch_signal_genes)

def mutation_TF(d_system, d_system_mutated, BS_length):
    nucleotides = ['a','c','g','t']
    #which TF?
    which_TF = random.sample(d_system['TFs'].keys(), 1)[0]
    #which nucleotide position?
    mutate_nucleotide_position = random.sample(range(BS_length), 1)[0]
    #and to which nucleatide?
    mutate_to_nucleaotide = random.sample(nucleotides, 1)[0]
    #mutated sequence
    to_list = list(d_system['TFs'][which_TF])
    to_list[mutate_nucleotide_position] = mutate_to_nucleaotide
    new_sequence = ''.join(to_list)
    #update system
    d_system_mutated['TFs'][which_TF] = new_sequence
    #
    return(d_system_mutated, which_TF)

def mutation_promoters(d_system, d_system_mutated, BS_length):
    nucleotides = ['a','c','g','t']
    #which promoter?
    which_promoter = random.sample(d_system['promoters'].keys(), 1)[0]
    #which nucleotide position?
    mutate_nucleotide_position = random.sample(range(BS_length), 1)[0]
    #and to which nucleatide?
    mutate_to_nucleaotide = random.sample(nucleotides, 1)[0]
    #mutated sequence
    to_list = list(d_system['promoters'][which_promoter])
    to_list[mutate_nucleotide_position] = mutate_to_nucleaotide
    new_sequence = ''.join(to_list)
    #update system
    d_system_mutated['promoters'][which_promoter] = new_sequence
    #
    return(d_system_mutated, which_promoter)

def mutation_signal(d_system, d_system_mutated, BS_length):
    #which TF?
    which_TF = random.sample(d_system['TF_signals'].keys(), 1)[0]
    #which signal?
    which_signal = random.sample(range(len(d_system['TF_signals'][which_TF])), 1)[0]
    #update system
    to_list = list(d_system['TF_signals'][which_TF])
    if d_system['TF_signals'][which_TF][which_signal] == '0':
        to_list[which_signal] = '1'
        new_sequence = ''.join(to_list)
    else:
        to_list[which_signal] = '0'
        new_sequence = ''.join(to_list)
    d_system_mutated['TF_signals'][which_TF] = new_sequence
    #
    return(d_system_mutated, which_TF)

def mutation_function(d_system):
    #copy system
    d_system_mutated = copy.deepcopy(d_system)
    #relative mutation rates
    r_tf = 1
    r_signal = 1
    r_promoters = 1
    #define probability of mutation happening
    p_mut_tf = len(d_system['TFs'])*r_tf #*BS_length
    p_mut_signal = len(d_system['TFs'])*r_signal
    p_mut_promoters = len(d_system['promoters'])*r_promoters #
    #
    possible_mutations = p_mut_tf + p_mut_signal + p_mut_promoters 
    #chose what will be mutated
    possible_mutations_name = ['TFs', 'signal' , 'promoters']#consensus sequence, signal specificity and binding sites
    probability_weights = [p_mut_tf/possible_mutations, p_mut_signal/possible_mutations, p_mut_promoters/possible_mutations]
    mutate = random.choices(possible_mutations_name,weights=probability_weights,k=1)[0]
    #
    #1: 'TFs', 2: 'TF_signals', 3: 'promoters'
    BS_length = 3
    if mutate=='TFs':
        (d_system, which_mutated) = mutation_TF(d_system, d_system_mutated, BS_length)
        mut_type = 1
    elif mutate=='signal':
        (d_system, which_mutated) = mutation_signal(d_system, d_system_mutated, BS_length) 
        mut_type = 2
    elif mutate=='promoters':
        (d_system, which_mutated) = mutation_promoters(d_system, d_system_mutated, BS_length) 
        mut_type = 3
    #
    return(d_system_mutated, mut_type, which_mutated)

def save_yaml(save_folder, save_file_name, d_initial_system, d_system, fitness_values, mutation_type, mutation_time, mutation_which_TF_or_promoter):
    d_save = {}
    d_save['initial_system'] = d_initial_system
    d_save['final_system'] = d_system
    d_save['fitness'] = [float(f) for f in fitness_values]
    d_save['mutation_type'] = mutation_type
    d_save['mutation_time'] = mutation_time
    d_save['mutation_which_TF_or_promoter'] = mutation_which_TF_or_promoter

    with open(save_folder+save_file_name, 'w') as file:
        yaml.dump(d_save, file)
        
def fix_mutation_function(delta_fitness,how_to_evaluate_delta_fitness, pop_size):
    #how_to_evaluate_delta_fitness = 'always_up' -> accept every mutation that increases fitness
    #how_to_evaluate_delta_fitness = 'probability_fixation' -> more realistic 
    fix_mutation = 0
    if how_to_evaluate_delta_fitness == 'always_up':
        if delta_fitness>0:
            fix_mutation = 1
    
    if how_to_evaluate_delta_fitness == 'probability_fixation':        
        if abs(delta_fitness)<0.000001:
            P_mutation_fixates = 1/(2*pop_size)
        elif -4*pop_size*delta_fitness>700: #I added this line because I was getting into numerical issues
            P_mutation_fixates = 0
        else:
            with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        P_mutation_fixates = (1-np.e**(-2*delta_fitness))/(1-np.e**(-4*pop_size*delta_fitness))
                    except Warning as e:
                        print('Houston, we have a warning:', e)
                        print('delta fitness: ', delta_fitness)
        #does mutation fixate?
        if random.random()<P_mutation_fixates:
            fix_mutation = 1
    
    return(fix_mutation)

def environments_function(n_genes):
    environments = [''.join(comb) for comb in product(['0','1'], repeat=n_genes)]
    return(environments)

def from_system_to_genotype(d_system):
    genotype = d_system['TF_signals'][0]+d_system['TFs'][0]+d_system['TF_signals'][1]+d_system['TFs'][1]+d_system['promoters'][0]+d_system['promoters'][1]
    return(genotype)

def evolution_simulation(seed_value,d_system,total_generations,general_save_folder,general_file_name,how_to_evaluate_delta_fitness, pop_size):
    ###
    # set seed
    random.seed(seed_value)
    
    ###
    os.makedirs(general_save_folder,exist_ok=True)
    file = open(general_save_folder+general_file_name, 'w')
    
    ###
    fitness_values = []

    ###
    # environments: selection in all possible environments
    n_genes = len(d_system['promoters'])
    environments = environments_function(n_genes)
    
    ###
    # generation 0
        
    d_initial_system = copy.deepcopy(d_system)
        
    #calculate mismatch between expected and observed gene expression
    all_p_binding, mismatch_signal_genes = mismatch_signal_geneExpression_function(d_system, environments)
    fitness = -1*mismatch_signal_genes
    fitness_values.append(fitness)
    
    #write to file
    genotype = from_system_to_genotype(d_system)
    file.write('0,'+genotype+','+str(all_p_binding)+','+str(fitness))

    ###
    # generation 1 and on
    for gen in range(1,total_generations):

        ###
        # create mutant system 
        d_system_mutated, mut_type, which_mutated = mutation_function(d_system) 

        ###
        #calculate mismatch between expected and observed gene expression in the mutant system
        all_p_binding, mismatch_signal_genes_mutant = mismatch_signal_geneExpression_function(d_system_mutated, environments)
        fitness_mutant = -1*mismatch_signal_genes_mutant
        
        ###
        #fix mutation?
        delta_fitness = fitness_mutant-fitness
        fix_mutation = fix_mutation_function(delta_fitness,how_to_evaluate_delta_fitness, pop_size)

        
        if fix_mutation==1:
            fitness = fitness_mutant
            #
            d_system = copy.deepcopy(d_system_mutated)
            genotype = from_system_to_genotype(d_system)
            #
            file.write('\n'+str(gen)+','+genotype+','+str(all_p_binding)+','+str(fitness))

        fitness_values.append(fitness)
    
    file.close()
     
    return(fitness_values)
        