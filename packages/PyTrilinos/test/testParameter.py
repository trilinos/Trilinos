#! /usr/bin/env python

import setpath
from PyTrilinos.NOX import Parameter

# Create a test parameter list
plTest = Parameter.List()
plTest.setParameter("Integer Parameter", 2003    )
plTest.setParameter("Float Parameter",   3.14    )
plTest.setParameter("String Parameter",  "Sandia")

# Check parameters' existence
assert(    plTest.isParameter("Integer Parameter"))
assert(    plTest.isParameter("Float Parameter"  ))
assert(    plTest.isParameter("String Parameter" ))
assert(not plTest.isParameter("Fake Parameter"   ))

# Check parameters' values
assert(plTest.getParameter("Integer Parameter") == 2003    )
assert(plTest.getParameter("Float Parameter"  ) == 3.14    )
assert(plTest.getParameter("String Parameter" ) == "Sandia")

# Emulate NOX example parameter list structure
nlParams     = Parameter.List()
nlParams.setParameter("Nonlinear Solver", "Line Search Based")
printParams  = nlParams.sublist("Printing")
printParams.setParameter("MyPID",            0); 
printParams.setParameter("Output Precision", 3);
printParams.setParameter("Output Processor", 0);
printParams.setParameter("Output Information", 
                         Parameter.Utils.Warning        +
                         Parameter.Utils.OuterIteration +
                         Parameter.Utils.InnerIteration +
                         Parameter.Utils.Parameters     +
                         Parameter.Utils.Details        +
                         Parameter.Utils.OuterIterationStatusTest);
searchParams = nlParams.sublist("Line Search")
searchParams.setParameter("Method", "Full Step")
dirParams    = nlParams.sublist("Direction")
dirParams.setParameter("Method", "Newton")
newtonParams = dirParams.sublist("Newton")
newtonParams.setParameter("Forcing Term Method", "Constant")
lsParams     = newtonParams.sublist("Linear Solver")
lsParams.setParameter("Aztec Solver",     "GMRES"                   )
lsParams.setParameter("Max Iterations",   800                       )
lsParams.setParameter("Tolerance",        1e-4                      )
lsParams.setParameter("Output Frequency", 50                        )
lsParams.setParameter("Scaling",          "None"                    )
lsParams.setParameter("Preconditioning",  "AztecOO: Jacobian Matrix")
print nlParams
