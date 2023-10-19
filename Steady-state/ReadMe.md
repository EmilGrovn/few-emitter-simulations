# Steady-State
Contains 4 folders:
## Low_pump_limit
In this folder, we consider the limit of low pumping.
## data
In this folder, data generated from the current folder is saved.
## plots
In this folder, plots generated from the current folder is saved.
## deprecated
In this folder, files, which are no longer of use are put.
## current folder
Contains scripts for steady state simulations, which are used to generate data for various pump values. The scripts for generating data are excactly analogous to the data generating scripts in the low pump folder. INSTEAD OF HAVING ONE SCRIPT FOR EACH PARAMETER IN OUR PARAMETER SPACE, CONSIDER CHANGING THE FUNCTIONALITY OF THE DATA GENERATING SCRIPTS, SO YOU CAN CHOOSE, WHICH OF THE 5 PARAMETERS (kappa, gammaA, gammaD, Pump. g) TO SWEEP OVER.

We go through the structure of the scripts below.

# Generating data
## perform_pump_sweep: 
in this script, you choose your system: 

physical input parameters (N_em,g,kappa,pump_logmin,pump_logmax,gamma,gamma2) (where logmin/max gives min/max value of pump/g, which we want to sweep over)

numerical parameters: cannot choose these here, only N_Hilbert. Consider changing this. Ties into our considerations

Then run one or more of the following functions:

pump_sweep_fixed_NH

pump_sweep_variable_NH

pump_sweep_spec

Let us take a look at each function

## pump_sweep_fixed_NH. Function.
input: physical parameters (pump is given as boundaries) and N_Hilbert

1. constructs pumplist from pump_logmin/max and an internal step parameter, dLogP_over_g

2. for each pump: call getSteady state. Calculate nP and g2

out: no output, but saves an numpyzip-file in the datafolder with nP (mean photon population) and g^(2)(0), and an np-array of corresponding values of pump/g

Further comments: this function is limited in the sense that for most pump-sweeps (when we reach lasing) a much bigger Hilbert space is needed than for the small pump
Preferably, we should alter it so it takes in a list of N_Hilbert

## pump_sweep_variable_NH. Function.
input: physical parameters (pump is given as boundaries)
1. constructs pumplist from pump_logmin/max and an internal step parameter, dLogP_over_g
2. for each pump: call getSteady state. Calculate nP and g2. 
      Do this until N_Hilbert is chosen such that a certain convergence criteria is met. We want to just meet the criteria.
      Consider which convergence criteria to use!! Right now, we have not thought enough about it
      Update: we have moved the old version of conv crit to the deprecated folder.
      We have now implemented and tested a convergence criterion where we consider the population in the higest photon state included
      
out: no output, but saves an numpyzip-file in the datafolder with nP (mean photon population), g^(2)(0), 
N_Hilbert used, Occupation propability for N_Hilbert photon population, and emitter occupation, and an np-array of corresponding values of pump/g

## getSteadyState. Function
input: physical parameters and N_Hilbert
1. Constructs Hamiltonian and collapse operators. Uses this to calculate steady state density operator, rho_ss
2. Constructs c_occ, c_occ^2 and emitter ocupation operators
return: rho_ss and c_occ, c_occ^2 and emitter ocupation operators

## pump_sweep_spec. Function.
input: physical parameters (pump is given as boundaries) and N_Hilbert
1. constructs pumplist from pump_logmin/max and an internal step parameter, dLogP_over_g
2. Use N_Hilbert
3. for each pump: call getSteadyStateSpectrum. Save spectrum
out: no outout, but saves an numpyzip-file in the datafolder with spectrum-sweep and frequencylist


## getSteadyStateSpectrum. Function
input: physical parameters and N_Hilbert
1. Constructs Hamiltonian and collapse operators. Uses this to calculate steady state density operator, rho_ss
2. Construct frequencylist. In this context, we choose number of points and range of frequencies
3. Use qutip function spectrum to calculate the spectrum. THIS IS BY FAR THE MOST COMPUTATIONALLY DEMANDING TASK
      Should consider: how do we choose frequencylist best
return: frequencylist and spectrum

## getLinewidth. Function
input: physical parameters

imports spectrum data and calculates and saves linewidth

# analysing data
## plot_linewidths: 
from generated data; makes a single plot of the linewidth with a lot of comparisons

## superplot_regimes
plots sweeps over photon number, the g2-function, the spectrum, and the linewidth with shared 1st axis of pumping strengths.

## comparison of N emitter to effective 1-emitter
Makes a plot of the same type as in Figure 5 in special course report.
