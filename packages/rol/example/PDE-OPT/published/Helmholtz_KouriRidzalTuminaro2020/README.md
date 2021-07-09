This subdirectory of the *ROL* package, i.e.,

/rol/example/PDE-OPT/published/Helmholtz_KouriRidzal2019

accompanies the paper

D. P. Kouri, D. Ridzal, R. Tuminaro, "KKT Preconditioners for
PDE-Constrained Optimization with the Helmholtz Equation"

To reproduce the numerical results, follow these steps:

1. Build Trilinos, with the ROL package and all its dependencies
   enabled.  For help building Trilinos, consult the Trilinos
   community or ROL developers.  In the main Trilinos build
   directory, you should see the subdirectory

   /packages/rol/example/PDE-OPT/published/Helmholtz_KouriRidzal2019

2. Obtain Matlab.  The code has been tested with Matlab R2019a.
   Start Matlab in command-line mode in the above subdirectory.
   Run
   >> study

3. To change study parameters, edit the study.m file.  The file
   contains descriptions of all parameters, and references to
   running the tests on finer meshes.

