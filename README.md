# QSPD
QSPD is used to solve the rate equations to obtain the **q**uantum **s**tates **p**opulation **d**ynamics. 


## Introduction
QSPD uses Einstein rate equation model to construct the rate equations, and solve it to get the population dynamics. <br>
The simplest case is a 2-level system namely there're two states with different energies in the system. According to the model, three types of transition can happen. <br>
  1. Absorption process - absorbing a photon to go to the higher-energy level
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/75/AtomicLineAb.svg/340px-AtomicLineAb.svg.png)
      
  2. Spontaneous emission process 
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/7a/AtomicLineSpEm.svg/339px-AtomicLineSpEm.svg.png)
      
  3. Stimulated emission process (the transition rate per unit density is described by B12)
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b9/AtomicLineInEm.svg/340px-AtomicLineInEm.svg.png)

![](http://latex.codecogs.com/png.latex?\dpi{110}&space;\begin{align}\begin{aligned}\frac{\partial\rho_i}{\partial&space;t}=-\sum_{j\neq&space;i}B_{ij}(I_{\mathrm{BBR}}&plus;I_{\mathrm{laser}})\rho_i&space;-&space;\sum_{j<i}A_{ij}\rho_i&space;\\&plus;\sum_{j\neq&space;i}B_{ji}(I_{\mathrm{BBR}}&plus;I_{\mathrm{laser}})\rho_i&space;&plus;&space;\sum_{j>i}A_{ji}\rho_i&space;\end{aligned}\end{align}" title="http://latex.codecogs.com/png.latex?\dpi{110} \begin{align}\begin{aligned}\frac{\partial\rho_i}{\partial t}=-\sum_{j\neq i}B_{ij}(I_{\mathrm{BBR}}+I_{\mathrm{laser}})\rho_i - \sum_{j<i}A_{ij}\rho_i \\+\sum_{j\neq i}B_{ji}(I_{\mathrm{BBR}}+I_{\mathrm{laser}})\rho_i + \sum_{j>i}A_{ji}\rho_i \end{aligned}\end{align})>

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
