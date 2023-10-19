The structure is the same as in the Steady state-folder. In this case, we just have scripts for each Lindblad decay rate sweep: kappa-sweep, gamma-sweep, gammaD-sweep. Let us give ourselves the full overview of the folder.

# Low pump limit
Contains 4 folders:
## Low_pump_limit_report
In this folder, we consider the limit of low pumping.
## data
In this folder, data generated from the current folder is saved.
## plots
In this folder, plots generated from the current folder is saved.
## Dynamical
This folder is not relevant to us anymore. Here we consider the transient case of zero pumping. We used this to check that it agrees with the spectrum for our steady state solutions in the limit of very low pumping.
## current folder
In this folder, we generate data through the exact same structure, as we did in the steady state folder. We have scripts for pump, gammaA, gammaD, and kappa.

We also have scripts for plotting the spectrum and the linewidth. Finally, a script-type called multiplot for linewidth and spectrum, where we plot 9 pump sweeps in one plot for 9 different values of either gammaA, gammaD, or kappa. These multiplots have been of limited use. In terms of plotting, the scripts in this folder are no longer relevant.
