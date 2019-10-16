
#############################
####  CompleteModel 
#############################

In this folder the files for modeling the developing of the RRMS disease considering the complete stochastic model. In the "NoDAC" folder are grouped the simulations without the therapy, and in "YesDAC" the same experiments (i.e.  using the same conditions) but with the  Daclizumab  therapy.

#############################
####  HybridModel 
#############################

Similar to above but considering the hybrid model.
Opening the ESSN model stored in the file .PNPRO it is possoble to observe which transitions are considered stochastic (labeled with FN:D: "transition's name") and which deterministic (labeled with just the "transition's name").

###################################
## How to generate the .solver:

ESSN -> name of the file ESSN.PNPRO

    - `unfolding2 ESSN -long-names`
    
    - `PN2ODE.sh ESSN -M -C transitions.cpp`


###################################
## How to run the .solver:

Usually:

    -  `./ESSN.solver NameOutput -type SolverType ... `
    
Since the injecions are modeled by 1) stopping the simulation, 2) modifyin the end marking of the ESSN adding the EBV cells or DAC cells, and 3) using the last marking modified as initial marking to continue the simulation.
The code is available in the following bash scripts:

- Resolution.sh
- ParallRes.sh

