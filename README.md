# few-emitter-simulations
In this project, we model a system of N emitters coupled to a cavity.

To do this, we use the master equation approach with the Born-Markov approkximation, where interaction with the environment are modelled using Lindblad terms.
Currently, we include cavity leakage through (with rate $\kappa$), background and non-radiative decay ($\gamma$), pure dephasing ($\gamma_D=\sqrt{2}\gamma_2$), and incoherent pumping ($P$).

The limitations of the Born-Markov approach are summarized briefly in e.g. the Nanophotonics lecture note.
