# Dynamical simulations

This folder is concerned with simulations of the time-evolution of our system. It contains the following codes:

## Population dynamics:
solveMasterEq_Nemitter: en funktion der tager fysiske (Nem, g, kappa, pump, gamma, gamma2) og numeriske (tlist, N_Hilbert) parametre som input – derpå bygger vores Hamiltonian og collapse operatorer, og udregner time-evolution result med mesolve - returner result

time_evolution: en funktion der tager fysiske (Nem, g, kappa, pump, gamma, gamma2) og numeriske (Tw, Nt, N_Hilbert) parametre som input – laver derpå tlist (ud fra Tw og Nt) bruger solveMasterEq_Nemitter – genererer en data-fil med cavity- og emitter-occuptation samt tlist – kan bruges til population dynamik

peform_time_evolution: tager fysike og numeriske parametre som input – kalder time_evolution¬ – indlæser genereret data-fil og laver time-evolution plots (nu også for N_em>1)

## spectrum:
compute2op2t_spectrum_trapz: input: fysiske og numeriske parametre. - Action: laver H og collapse operatorer (samme som solveMasterEq). Skal vælge passende omega-liste. Udregner 2-time correlation function med correlation_2op_2t, som giver matrix. Udregner spectrum ved at integrere numerisk to gange vha. np.trapz. – return: wlist, S (spectrum)

perform_spectrum: kalder compute2op2t_spectrum_trapz: udregner linewidth. Gemmer wlist, spectrum, linewidth 

plot_spectrum: indlæser wlist, spectrum, linewidth. Plotter spektrum