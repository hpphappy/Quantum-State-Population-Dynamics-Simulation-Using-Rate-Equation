#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 13:45:47 2019

@author: panpanhuang
"""

import numpy as np
import pandas as pd


PgopherLine = np.genfromtxt('AlH+_A-X_transition_transition_v1_Circular_JMax50_Threshold-6.csv', delimiter = ',' , skip_header = 2, skip_footer =1, usecols=(2,3,4,6,7,8,9,10,11,12,13,14,15))
JEx = PgopherLine[:,0]
RParity_Ex = PgopherLine[:,1]
NEx = PgopherLine[:,2]
JGnd = PgopherLine[:,3]
RParity_Gnd = PgopherLine[:,4]
NGnd = PgopherLine[:,5]
Transition = PgopherLine[:,6]
OStrength = PgopherLine[:,7]
Eupper = PgopherLine[:,8]
Elower = PgopherLine[:,9]
Spol = PgopherLine[:,10]
A = PgopherLine[:,11]
Width = PgopherLine[:,12]

ElowerSort = PgopherLine[np.argsort(PgopherLine[:,9])]
EupperSort = PgopherLine[np.argsort(PgopherLine[:,8])]

#Generate an array of state energy including the energy of GND and the EX electronic states
StateEnergy = []
StateEnergy_temp = np.concatenate((ElowerSort[:,9],EupperSort[:,8]), axis = None)

J_sort_by_Energy = []
J_sort_by_Energy_temp = np.concatenate((ElowerSort[:,3],EupperSort[:,0]), axis = None)

J_degeneracy_sort_by_Energy=[]

#Generate the transition energy matrix and write to a csv file
for i in np.arange(len(StateEnergy_temp)):
    if StateEnergy_temp[i] not in StateEnergy:
        StateEnergy.append(StateEnergy_temp[i])
        J_sort_by_Energy.append(J_sort_by_Energy_temp[i])
        J_degeneracy_sort_by_Energy.append(2*J_sort_by_Energy_temp[i]+1)

TransitionEnergy = np.zeros((len(StateEnergy), len(StateEnergy)))
for i in np.arange(len(StateEnergy)):
    for j in np.arange(len(StateEnergy)):
        TransitionEnergy[i,j] = StateEnergy[j]-StateEnergy[i]

np.savetxt("v1_J_degeneracy_sort_by_Energy_Threshold-6.csv", J_degeneracy_sort_by_Energy, delimiter=',')        
np.savetxt("v1_TransitionEnergy_Threshold-6.csv.csv", TransitionEnergy, delimiter=',')

#Gerate the oscillator strength matrix and write to a csv file
OscillatorStrength = np.zeros((len(StateEnergy), len(StateEnergy)))
for i in np.arange(PgopherLine.shape[0]):
    OscillatorStrength[StateEnergy.index(Elower[i]), StateEnergy.index(Eupper[i])] = OStrength[i]
    OscillatorStrength[StateEnergy.index(Eupper[i]), StateEnergy.index(Elower[i])] = OStrength[i]
    
np.savetxt("v1_OscillatorStrength_HC_Threshold-6.csv.csv", OscillatorStrength, delimiter=',')
