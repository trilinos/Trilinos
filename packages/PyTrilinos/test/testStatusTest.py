#! /usr/bin/env python

import setpath
from PyTrilinos.NOX import StatusTest

# Test NormF status test
absresid = StatusTest.NormF(1.0e-6)

# Test NormF status test (unscaled)
relresid = StatusTest.NormF(1.0e-6,StatusTest.NormF.Unscaled)

# Test NormUpdate status test
update = StatusTest.NormUpdate(1.0e-5)

# Test NormWRMS status test
wrms = StatusTest.NormWRMS(1.0e-2, 1.0e-8)

# Test MaxIters status test
maxiters = StatusTest.MaxIters(20)

# Stagnation status test
stagnate = StatusTest.Stagnation()

# Finite value status test
finiteVal = StatusTest.FiniteValue()

# Test the Combo StatusTest constructors and addStatusTest member function
converged = StatusTest.Combo(StatusTest.Combo.AND)
converged.addStatusTest(absresid)
converged.addStatusTest(relresid)
converged.addStatusTest(wrms)
converged.addStatusTest(update)
combo = StatusTest.Combo(StatusTest.Combo.OR)
combo.addStatusTest(converged)
combo.addStatusTest(maxiters)
print combo
