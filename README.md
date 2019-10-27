# Quantum-State-Population-Dynamics-Simulation-Using-Rate-Equation
Solving rate equations to obtain the quantum states population dynamics. 

To construct the rate equations, you need the following excel or txt files.

1. the transition energy matrix, TransE[i,j] = StateEnergy[i]-StateEnergy[j]. 
    The example file is "AlH+_v1_TransitionEnergy_Threshold-6.csv" in the repository.
2. the oscillator strength matrix
    The example file is "AlH+_v1_J_degeneracy_sort_by_Energy_Threshold-6" in the repository
3. the array of state degeneracies.
    The example file is "AlH+_v1_J_degeneracy_sort_by_Energy_Threshold-6.csv"
    
Note: the states are sort by state energy to generate the oscillator strength matrix and the state degenracies array.

If you have a Pgopher file constructed, you can change the unit of the intensity in Pgopher to oscillator strength and export a transition line file from Pgopher. Use the file "TransItion_Info_from_Pgopher.ipynb" to generate the three required files above using the exported Pgopher file (example pgopher output file: AlH+_A-X_transition_v1_0Gauss_Random_JMax50_Threshold-6.csv).
