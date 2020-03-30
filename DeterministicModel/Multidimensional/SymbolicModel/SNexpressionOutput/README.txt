INSTRUCTIONS:

Download the SNexpression tool from http://www.di.unito.it/~depierro/SNexpression/#Download

- unzip the file SNex.zip;
- double-click the jar file SNexpression.jar or in a OS shell change the working directory to the folder where the archive has been decompressed and type the following command:
java -jar SNexpression.jar

A command line interface (CLI) will open. 
- In the CLI command window type
load "pathname/model-filename.sn" 
where pathname is the path to reach the file model-filename.sn
(e.g. load "./sn/SMmodelCIBB19.sn"
- type print_ode

New files are generated in the same directory where the model file is, among which:
model-filename.sn.gen.ode  (containing the system of SODE in R syntax, to be completed)
model-filename.sn.gen.map  (to be filtered as explained hereafter)

execute the bash script

script_prepare model-filename

(Note the script needs the files awk_filter_func and awk_NoDup)


The script generates a new file 
model-filename.sn.genNODUP.map obtained by filtering the contents of  model-filename.sn.gen.map

The new file contains the definition of the functions used in the SODE formulae, assuming all transitions are of type MASS-ACTION. The user needs to redefine the functions for the transitions of type GENERAL.

