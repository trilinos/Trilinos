#  NOX: An Object-Oriented Nonlinear Solver Package and LOCA: Library of Continuation Algorithms Package

## Introduction

The User's Guide and Developer's Guide for NOX and LOCA is combined on
a single website found at:

  http://trilinos.sandia.gov/packages/nox

Alternatively, you can build the documentation as html pages yourself
from the source code using the doxygen program.  To build the
documentation type the following:

1. Go to the top level trilinos build directory and configure a bare
   bones build of trilinos by typing:

   configure --enable-teuchos --enable-nox

2. Change into the top of the nox build directory:

   cd packages/nox

3. Create the documentation by typing the following at the top level nox
   build directory:

   make dox

4. Finally, bring up the following web page in your web browser for the
   remaining instructions.  The documentation is NOT generated in the
   build tree but rather in the source tree:

   Trilinos/packages/nox/doc/html/index.html


## Copyright and License
See nox/COPYRIGHT\_NOX, nox/LICENSE\_NOX, nox/COPYRIGHT\_LOCA, nox/LICENSE\_LOCA, https://trilinos.github.io/license.html and individual file headers for additional information.


## Questions? 
Contact lead developers:

* NOX team        (GitHub handle: @trilinos/nox)
* Roger Pawlowski (GitHub handle: [rppawlo](https://github.com/rppawlo) or rppawlo@sandia.gov)
* Eric Phipps     (GitHub handle: [etphipp](https://github.com/etphipp) or etphipp@sandia.gov)

