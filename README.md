# QSPD
QSPD is used to solve the rate equations to obtain the **q**uantum **s**tates **p**opulation **d**ynamics. 


## Introduction
QSPD uses Einstein rate equation model to construct the rate equations, and solve it to get the population dynamics. <br>
The simplest case is a 2-level system namely there're two states with different energies in the system. According to the model, three types of transition can happen. <br>
  1. Absorption process - absorbing a photon to go to the higher-energy level
  2. Spontaneous emission process 
      * 
  3. Stimulated emission process (the transition rate per unit density is described by B12)
      * ![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/75/AtomicLineAb.svg/340px-AtomicLineAb.svg.png)
The population in those 2 states

## Usage
* Prepare necessary files for constructing the rate equation
To construct the rate equations, you need the following excel or txt files.

1. the transition energy matrix, TransE[i,j] = StateEnergy[i]-StateEnergy[j]. 
    The example file is "AlH+_v1_TransitionEnergy_Threshold-6.csv" in the repository.
2. the oscillator strength matrix
    The example file is "AlH+_v1_J_degeneracy_sort_by_Energy_Threshold-6" in the repository
3. the array of state degeneracies.
    The example file is "AlH+_v1_J_degeneracy_sort_by_Energy_Threshold-6.csv"
    
Note: the states are sort by state energy to generate the oscillator strength matrix and the state degenracies array.

If you have a Pgopher file constructed, you can change the unit of the intensity in Pgopher to oscillator strength and export a transition line file from Pgopher. Use the file "TransItion_Info_from_Pgopher.ipynb" to generate the three required files above using the exported Pgopher file (example pgopher output file: AlH+_A-X_transition_v1_0Gauss_Random_JMax50_Threshold-6.csv).
