#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 10:55:34 2019

@author: panpanhuang
"""

import os
import numpy as np
import csv
from scipy.integrate import odeint
import matplotlib
matplotlib.use('pdf')
from matplotlib import gridspec
import matplotlib.pyplot as plt
import time

t1=time.time()
Molecule = "AlH"


#Import folder and import file names
ImportDir = r""
TransitionE = r"AlH+_v1_TransitionEnergy_Threshold-6.csv"
Degeneracies = r"AlH+_v1_J_degeneracy_sort_by_Energy_Threshold-6.csv"
OStrength = r"AlH+_v1_OscillatorStrength_Threshold-6.csv"
TranLabel = r"AlH+_A-X_transition_v1_random_polarization"

#Output path and output file names
outputpath = r""
if not os.path.exists(outputpath):
    os.makedirs(outputpath)
PopDyn = str(TranLabel) + r"_popdyn.csv"


# Define physical and environmental constants (c is in cm/s; k_B is in cm-1/K; T is in K)
q_e = 1.60E-19
m_e = 9.11E-31
c = 3.00E+10
h = 6.626E-34
eps_0 = 8.85E-12
k_B = 0.695
T = 295

# Define simulation granularity:
# Time sampling (Default: 100,000 step/1 ms)
t_i = 0.00E-00
t_f = 400E-03
t_steps = 1E+04

# Energy level bin width
Delta_E = 0.00

# Define Rotational cooling laser parameters 
## Define frequency-domain parameters - MaiTai laser parameters: 80-MHz doubled-NIR, laser with ~200-cm-1 FWHM)
## Calculate the averaged laser intensity of the pulsed laser
## I_0 = 85.52463783*(the real intensity). The unit is W/cm^2
Off_On_tRatio = 12.5E-9/1.6E-13 
I_0 = 2.00E-9/Off_On_tRatio

E_nu_tilde_0 = 27620.00
## FWHM bandwidth ~ 2.355*sigma (200 cm-1 FWHM ~ 85 cm-1 sigma)
E_nu_tilde_sigma = 85.00
## The spacing (in Hz and in cm^-1) between the comb teeth of the MaiTai laser output
E_nu_rep = 80.00E06
E_nu_tilde_rep = E_nu_rep/c
print("Comb teeth spacing = " + str(E_nu_tilde_rep) + " cm-1")
## Distribute the comb teeth within the 2*sigma range from the center wavenumber
num_combteeth = np.int((4*E_nu_tilde_sigma)/E_nu_tilde_rep)+1 
E_nu_tilde_rep_mask = np.linspace(E_nu_tilde_0 - 2*E_nu_tilde_sigma, E_nu_tilde_0 + 2*E_nu_tilde_sigma, num=num_combteeth)
## Define the laser intensities at the positions of combteeth
I_w_0 = I_0*np.exp(-np.square((E_nu_tilde_rep_mask - E_nu_tilde_0)/(np.sqrt(2)*E_nu_tilde_sigma)))

#Define the pulse-shaping parameters
##Define removed frequency range from the MaiTai output through pulse-shaping
notchfilters = [(27655.00,100000.00)]

## Apply pulse-shaping mask to the laser intensities using the specified notchfilter(s)
I_w = I_w_0
PSmask = np.ones_like(I_w_0)
for n in notchfilters:
    cuton_indices = np.nonzero(E_nu_tilde_rep_mask >= n[0])
    cutoff_indices = np.nonzero(E_nu_tilde_rep_mask <= n[1])
    nnum = cutoff_indices[0][-1] - cuton_indices[0][0]+1
    notch = np.linspace(cuton_indices[0][0],cutoff_indices[0][-1],num=nnum,dtype=int)
    for i in range(notch.shape[0]):
        PSmask[notch[i]] = 0.00
        I_w = I_w_0*PSmask


#Import files
TransitionEnergies = np.loadtxt(ImportDir + TransitionE, delimiter=",")
LevelEnergies = TransitionEnergies[0,:]
LevelDegeneracies = np.loadtxt(ImportDir + Degeneracies, delimiter=",")
OscillatorStrengths = np.loadtxt(ImportDir + OStrength, delimiter=",")
K = TransitionEnergies.shape[0]

"""
Oscillator strengths imported from pgopher include degeneracies!
"""
#Generate Einstein AB coefficiens
LevelDegeneracyRatios = np.ones((K,K))
A_if_coeffs = ((2*np.pi*1.00E+06*(TransitionEnergies**2)*(q_e**2))/(eps_0*m_e*c))*LevelDegeneracyRatios*OscillatorStrengths
np.fill_diagonal(TransitionEnergies, 1.00)
B_if_coeffs = ((q_e**2)/(4*eps_0*m_e*h*c*np.abs(TransitionEnergies)))*LevelDegeneracyRatios*OscillatorStrengths
B_fi_coeffs = ((q_e**2)/(4*eps_0*m_e*h*c*np.abs(TransitionEnergies)))*OscillatorStrengths        


#MaiTai Laser comb teeth resonance intensity test station
##Generate D in which each element is the summation of the total intensity the comb teeth that satify
##the resonance condition for each transition
## Presumed transition linewidths (cm^-1, 0.0053 cm^-1 ~ 160 MHz)
transition_nu_tilde_sigma = E_nu_tilde_rep;
DrivenLaser = np.zeros((K,K))
for i in range(K):
    for j in range(i,K):
        print("i=" + str(i) + " --> " + "f=" + str(j))
        cuton_indices = np.nonzero(E_nu_tilde_rep_mask >= (np.abs(TransitionEnergies[i,j]) - transition_nu_tilde_sigma))
        cutoff_indices = np.nonzero(E_nu_tilde_rep_mask <= (np.abs(TransitionEnergies[i,j]) + transition_nu_tilde_sigma))
        if cutoff_indices[0].shape[0] != 0 and cuton_indices[0].shape[0] != 0:
            inum = cutoff_indices[0][-1] - cuton_indices[0][0]+1
            inclusion = np.linspace(cuton_indices[0][0],cutoff_indices[0][-1],num=inum,dtype=int)
            for index in inclusion:
                DrivenLaser[i,j] += I_w[index]
                DrivenLaser[j,i] = DrivenLaser[i,j]

# Determine FWHM bounds (TransitonEnergy -/+ 2.414*sigma)
# Determine if and how many comb frequencies fall within this bandwidth
# Sum together intensities of those comb frequencies (taking into account their
# respective detunings) to arrive at an effective driving strength
# Assume Q of axial modes is infinite such that the comb lines are delta functions?
# (So that each comb line has an "area" that is approximated by its height)
# DOES NOT TAKE INTO ACCOUNT NON-RESONANT CONTRIBUTIONS!!!

# Use this array as a mask to add driven terms to the K matrix
# Calculate BBR spectral intensity (https://nanohub.org/wiki/DerivationofPlancksLaw)
# The 'delta_nu_tilde' parameter is the bandwidth of the transition (for now, assumed to be 1)
                
##BBR field driven intensity
delta_nu_tilde = 1
np.fill_diagonal(TransitionEnergies, 1.00)
BBR = ((8*h*np.pi*1.00E+06*np.abs(TransitionEnergies)**3)/(np.exp(np.abs(TransitionEnergies)/(k_B*T)) - 1))*delta_nu_tilde
np.fill_diagonal(TransitionEnergies, 0.00)

# Define the field intensities matrix (Laser field + BBR field)
FIntensity = BBR + DrivenLaser


# From the field intensities matrix, calculate the rate matrix contributed by the field
FieldStimulated_lower_feedMatrix = np.tril(B_if_coeffs, k=-1)*np.tril(FIntensity, k=-1)    
FieldStimulated_diagonal_depleteMatrix = np.zeros((K,K), dtype = float)
for i in range(K):
    FieldStimulated_diagonal_depleteMatrix[i,i] = np.sum((np.triu(B_fi_coeffs,k=1)*np.triu(FIntensity,k=1) + np.tril(B_fi_coeffs,k=-1)*np.tril(FIntensity,k=-1))[i,:])
FieldStimulated_upper_feedMatrix = np.triu(B_fi_coeffs, k=1)*np.triu(FIntensity, k=1)
FieldStimulated_tot = FieldStimulated_lower_feedMatrix - FieldStimulated_diagonal_depleteMatrix + FieldStimulated_upper_feedMatrix


#Spontanieous decay rate
## Spontaneous emission TO states of lower energy DEPLETE the given state
## The rows of the K matrix are defined in ascending order of energy
## For a given K-row index, I SUBTRACT a term that corresponds to contributions
## from ALL SMALLER indices
## For state (column) i, just need to sum all rows of 'A_if_coeffs' up to the (i-1)th element!
A_diagonal_deplete = np.zeros((K,K), dtype = float)
for i in np.linspace(1,K-1,num=K-1,dtype=int):
    A_diagonal_deplete[i,i] = np.cumsum(A_if_coeffs[:][i])[i-1]

## Spontaneous emission FROM states of higher energy FEED the state
## The rows of the K matrix are defined in ascending order of energy
## For a given K-row index, I define positive terms that correspond to contributions
## from each of the LARGER indices
A_offdiagonal_feed = np.triu(A_if_coeffs, k=1)
A_tot = A_offdiagonal_feed - A_diagonal_deplete

#Combine Field driven terms with spontaneous decay rate
DrivenTot = FieldStimulated_tot + A_tot

# Define initial values to evaluate constants of the general solution
# (i.e. define initial state populations as Boltzmann distribution)
#"""
Q = np.sum(LevelDegeneracies*np.exp(-LevelEnergies/(k_B*T)))
y0 = LevelDegeneracies*np.exp(-LevelEnergies/(k_B*T))/Q

# Numerical Integration
# Define system of ordinary differential equations
def f(y, t, K):
    return np.dot(DrivenTot,y)

t = np.linspace(t_i,t_f,num=int(t_steps + 1),dtype=float)

# Solve system of equations via numerical integration
print("Numerically-integrating the system...")
yt = np.transpose(odeint(f, y0, t, args=(K,), hmin=0, mxstep=0))
    
print(yt[:,0]-yt[:,-1])
t2=time.time()
print(t2-t1,'s')


print("Initial total probability=" + str(np.sum(yt[:,0])))
print("Midpoint total probability=" + str(np.sum(yt[:,int(yt.shape[1]/2)])))
print("Final total probability=" + str(np.sum(yt[:,-1])))


# Determine states with population above specified threshold
pop_threshold = 0.03
nonzerosteadystates = np.nonzero(yt[:,-1] >= pop_threshold)[0]
if len(nonzerosteadystates) > 0:
    print("At the end of the simulation time, the following state populations are above the specified threshold of " + str(pop_threshold) + ":\n")
    for n in nonzerosteadystates:
        enrg = LevelEnergies[n]
        pop = yt[n,-1]
        print("P_" + str(n) + "=" + str(pop) + "; level energy=" + str(enrg) + " cm-1" + "\n")
else:
    print("At the end of the simulation time, no state populations are above the specified threshold of " + str(pop_threshold) + ":\n")

# Export system dynamics to a CSV file
with open(outputpath + PopDyn,"w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(yt)

# Plot the above-threshold populations as a function of time
# Global plot format settings
Fontsize = 10
LW = 0.500
LS = "--"
GridMajorLW = 0.25
GridMinorLW = 0.15

if notchfilters == []:
    subtitle_state = "OFF"
else:
    subtitle_state = "ON"

# Plot the quantum state population dynamics
Figure2 = plt.figure(figsize=(10,5))
gs = gridspec.GridSpec(1, 2, width_ratios=[7, 1]) 
ax0 = Figure2.add_subplot(gs[0])
for n in nonzerosteadystates:
    LW = 1.00
    label = str(n)
#    j = r"|\," + label + r"\left.\right\rangle|\,masked\,\left.\right\rangle"
    j = r"|\," + label + r"\left.\right\rangle"
    formula =  "${{{}}}$".format(j)
    ax0.plot(t*1.00E+06,yt[n,:],label=formula, markersize=3, alpha=0.5, linewidth=LW)
ax0.legend(loc=3, bbox_to_anchor=(1.0,0.0,0.5,0.5),fontsize=Fontsize,ncol=1)

plt.title(r"Population Dynamics of the Driven %s$^+\,$X(v''=0,1)$\,\leftrightarrow\,$A(v'=0) Manifold" %Molecule,fontsize=Fontsize+2)
plt.xlabel("Time ($\mu s$)",fontsize=Fontsize)
plt.ylabel("Probability",fontsize=Fontsize)
plt.ylim(-0.10,1.10)
Figure2.savefig(outputpath + "%splus_XA_pop_dyn_"%Molecule + "filter_" + str(subtitle_state) + ".pdf")
print('%.6g' %np.max(FieldStimulated_tot))