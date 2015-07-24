Muemex is the MATLAB interface for MueLu.

<><><><><><><><><><>
Basic Instructions:
<><><><><><><><><><>

-Run "matlab" script in this directory to launch matlab with proper shared libraries
-Run "help muelu" from matlab to see detailed help for "muelu" function
-Basic setup for muelu is "problemID = muelu('setup', A);"
-Basic solve for muelu is "x = muelu(problemID, b);"
-Run "ctest" in this directory to run the experimental matlab tests for MueLu

<><><><><><><>
Muemex Usage:
<><><><><><><>

Go to muelu/matlab/bin.
Run MATLAB through the "matlab" script.

With a sparse matrix A, set up a MueLu hierarchy:

    problemID = muelu('setup', A);

Note: any number of problems can be set up at a time.
Optionally, pass A and fine level coordinates:

    problemID = muelu('setup', A, coords);

Parameters for the "easy parameter list" are listed after A and coords:

    problemID = muelu('setup', A, coords, 'coarse: max size', 50);

Sublists can be passed using MATLAB cell arrays.

    problemID = muelu('setup', A, 'level 0', {'aggregation: drop tol', 0.03}, 'level 1', {'aggregation: drop tol', 0.01});

Parameters can be strings, integers, booleans or arrays (full or sparse, real
or complex). 

Solve a problem with b as the right-hand side vector(s):

    x = muelu(problemID, b);
    x = muelu(problemID, A, b);

The first version uses the A that was passed when setting up the problem.
The second version uses a new A with the same hierarchy.

Remove a specific problem and free memory associated with it:

    muelu('cleanup', problemID);

or all problems:

    muelu('cleanup');

List basic information about problem(s):

    muelu('status', problemID);
    muelu('status');

Get data from a level:

    data = muelu('get', problemID, levelID, 

<><><><><><><><><><><>
Factory Instructions:
<><><><><><><><><><><>

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
level. "Level" is a special key for the Needs list, and will be passed to MATLAB with the
current level ID.
-"Provides" is what will be returned by the matlab function and added to the
level.
-"Function" is the name of the matlab function to run. Parameters/return
values must match Needs/Provides.

TwoLevelMatlabFactory: "Needs Fine", "Needs Coarse", "Provides", "Function"
-Just like SingleLevelMatlabFactory, but inputs come from both fine and coarse
levels.

<><><><><><><><><>
Custom Variables:
<><><><><><><><><>

Muemex also supports setting custom data in the hierarchy. To use this
feature, set the data you want as a parameter in a level sublist when
setting up the problem. For example, if you want a matrix "MyMatrix" to
be available to a matlab factory, set up with this command:

muelu('setup', A, 'level 0', {'Matrix MyMatrix', MyMatrix}, 'xml parameter
file', 'myXMLParams.xml');

Notice the type specifier "Matrix" before the variable name. This is part
of the name, and has to be included in every mention of the variable. The
type is not case sensitive; 'MultiVector' and 'multivector' are equivalent.
The name itself can be any string without whitespace. Names have to unique
within their level - there can't be "Double d1" and "Int d1".

If an aggregates factory needs that custom matrix, the XML parameter list for
the problem might look like this:
...
<Parameter name="factory" type="string" value="SingleLevelMatlabFactory"/>
<Parameter name="Provides" type="string" value="Aggregates"/>
<Parameter name="Needs" type="string" value="A, Matrix MyMatrix"/>
<Parameter name="Function" type="string" value="simpleAggregation"/>
...

This will send MyMatrix to simpleAggregation as the second argument.
There must always be exactly one space between the
type and the name. Leading and trailing spaces are always ignored.

The following custom variable types are supported:
-Matrix (sparse, MxM, real or complex depending on MueLu context)
-MultiVector
-OrdinalVector (must be Mx1 column vector of int32)
-Int (32 bit signed)
-Scalar (double or complex depending on context)
-Double
-Complex

Note: Custom variables added to the hierarchy either through muelu('setup') or
by a matlab factory are never removed from their level - a UserData keep
flag is set for them.
