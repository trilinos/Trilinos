Muemex is the MATLAB interface for MueLu.

Basic Instructions:
-Run "matlab" script in this directory to launch matlab with proper shared libraries
-Run "help muelu" from matlab to see detailed help for "muelu" function
-Basic setup for muelu is "problemID = muelu('setup', A);"
-Basic solve for muelu is "x = muelu(problemID, b);"
-Run "ctest" in this directory to run the experimental matlab tests for MueLu.

Factory Instructions:
-MatlabSmoother, SingleLevelMatlabFactory and TwoLevelMatlabFactory are
implementations of those types of factories that use matlab functions
to generate level data.

Parameters:
MatlabSmoother: "Needs", "Setup Function", "Solve Function", "Number of Solver
Args"
-"Needs" is a comma-separated list of hierarchy data that the smoother will
request, and pass into setup function in order that they are listed.
-"Setup Function" is the matlab function to run to set up the smoothing. Must
take a sparse matrix (A) followed by the matlab types corresponding to "Needs"
-"Solve Function" is the matlab function to run the solve phase of the
smoother. Must take sparse matrix A, double array X, double array B, followed by
the outputs of "Setup Function" as arguments. Shouldn't return anything.
-"Number of Solver Args" (int) is the number of expected outputs of "Setup
Function" that will be stored until the solve phase, when they will be
passed in after A, x, b.

SingleLevelMatlabFactory: "Needs", "Provides", "Function"
-"Needs" is list of inputs to the matlab function that will be pulled from
level.
-"Provides" is what will be returned by the matlab function and added to the
level.
-"Function" is the name of the matlab function to run. Parameters/return
values must match Needs/Provides.

TwoLevelMatlabFactory: "Needs Fine", "Needs Coarse", "Provides", "Function"
-Just like SingleLevelMatlabFactory, but inputs come from both fine and coarse
levels.
