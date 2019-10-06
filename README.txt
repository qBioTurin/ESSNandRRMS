
The ESSN model is stored into the file ModelRRMS.PNPRO (the first net is the one with the general functions, and the second without them)

By using these commands is possible to generate the system of ODEs stored in a R file or in a C++ solver

1) Unfolding of the ESSN in order to obtain the underline PN:

- unfolding2 ModelRRMS -long-names

2) Generation of the ODEs system:

- PN2ODE.sh ModelRRMS -M -R   (R file, in the example ModelRRMS.R we added the general functions)
- PN2ODE.sh ModelRRMS -M -C transitions (C++ solver, the file transitions.txt stores the general transitions)

 
----------------------------------------------------

To reproduce the results presented in the paper, the following files are provided:

1) Resolution.R -> for the Healthy and Unhealthy cases
2) PregnatResolution.R -> for the pregnancy 
