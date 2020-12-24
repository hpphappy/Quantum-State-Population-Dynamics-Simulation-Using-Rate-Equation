# QSPD
QSPD is used to solve the rate equations to obtain the **q**uantum **s**tates **p**opulation **d**ynamics. 


## Introduction
QSPD uses Einstein rate equation model to construct the rate equations, and solve it to get the population dynamics. <br>
The simplest case is a 2-level system namely there're two states with different energies in the system. According to the model, three types of transition can happen. <br>
  1. Absorption process - absorbing a photon to go to the higher-energy level
  
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/75/AtomicLineAb.svg/340px-AtomicLineAb.svg.png)
      
  2. Spontaneous emission process 
  
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/7/7a/AtomicLineSpEm.svg/339px-AtomicLineSpEm.svg.png)
      
  3. Stimulated emission process
  
      ![](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b9/AtomicLineInEm.svg/340px-AtomicLineInEm.svg.png)

![](http://latex.codecogs.com/png.latex?%5Cdpi%7B110%7D%20%5Cbegin%7Balign%7D%5Cbegin%7Baligned%7D%5Cfrac%7B%5Cpartial%5Crho_i%7D%7B%5Cpartial%20t%7D=-%5Csum_%7Bj%5Cneq%20i%7DB_%7Bij%7D(I_%7B%5Cmathrm%7BBBR%7D%7D&plus;I_%7B%5Cmathrm%7Blaser%7D%7D)%5Crho_i%20-%20%5Csum_%7Bj%3Ci%7DA_%7Bij%7D%5Crho_i%20%5C%5C&plus;%5Csum_%7Bj%5Cneq%20i%7DB_%7Bji%7D(I_%7B%5Cmathrm%7BBBR%7D%7D&plus;I_%7B%5Cmathrm%7Blaser%7D%7D)%5Crho_i%20&plus;%20%5Csum_%7Bj%3Ei%7DA_%7Bji%7D%5Crho_i%20%5Cend%7Baligned%7D%5Cend%7Balign%7D)

where ![](http://latex.codecogs.com/png.latex?\dpi{110}&space;\rho_i) is the population fraction in state i. The system of equations includes the rovibronic and hyperfine states of interest. The initial population was assumed to be thermal with a temperature of 300 K. ![](http://latex.codecogs.com/png.latex?\dpi{110}&space;I_{BBR}) and ![](http://latex.codecogs.com/png.latex?\dpi{110}&space;I_{laser}) are the energy densities of the blackbody radiation and laser. ![](http://latex.codecogs.com/png.latex?\dpi{110}&space;A) and ![](http://latex.codecogs.com/png.latex?\dpi{110}&space;B) are the spontaneous emission and stimulated emission Einstein coefficients, respectively.


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
