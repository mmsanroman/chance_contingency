#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 07:09:30 2023

@author: magdalena
"""

#import packages
import random
from itertools import product
import numpy as np
import chance_contingency_module as mym
import os
import yaml
import sys
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"

BS_length = 3
d_sequence_space = mym.sequence_space(BS_length)

replicates = 100

total_generations = 5000
pop_size = 100
general_save_folder = '/Users/magdalena/polybox/research/chance_contingency/results/pop_size_'+str(pop_size)+'/'
how_to_evaluate_delta_fitness = 'probability_fixation' 

for word in d_sequence_space[BS_length][1:]:
    sequences = [word, word]
    d_system = mym.create_initial_system(sequences)
    #print(d_system)
    
    genotype = mym.from_system_to_genotype(d_system)
    #print(genotype)
    
    for r in range(replicates):
        seed_value = r
        general_file_name = genotype+'_'+str(seed_value)
        
        fitness_values = mym.evolution_simulation(seed_value,d_system,total_generations,general_save_folder,general_file_name,how_to_evaluate_delta_fitness, pop_size)
        
        #plt.plot(range(total_generations),fitness_values)
        #plt.xlabel('Generations')
        #plt.ylabel('Fitness')


