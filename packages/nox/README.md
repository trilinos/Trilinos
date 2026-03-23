#  NOX: An Object-Oriented Nonlinear Solver Package and LOCA: Library of Continuation Algorithms Package

NOX and LOCA are a combined package for robustly solving and analyzing large-scale systems of nonlinear equations.

NOX is short for Nonlinear Object-Oriented Solutions, and its objective is to enable the robust and efficient solution of the equation: $ F(x)=0 $, where $ F:\Re^n \rightarrow \Re^n $ using globalized Newton methods such as line search and trust region methods. NOX is designed to work with any linear algebra package and to be easily customized. NOX is part of Sandia’s Trilinos project.

LOCA, distributed as part of NOX, is short for Library of Continuation Algorithms, and its objective is to compute families of solutions to \f$ F(x,p)=0 $ and their bifurcations, where $ F:\Re^n\times\Re^m\rightarrow\Re^n $.

The User’s Guide and Developer’s Manual are built using doxygen and can be accessed via the web. Please choose the documentation corresponding to your software version on the left hand toolbar

## Overview

**NOX**  
NOX is short for Object-Oriented Nonlinear Solutions. NOX is an object-oriented C++ library designed to solve large-scale systems of nonlinear equations. It implements Newton-based globalization techniques including line search and trust region algorithms. NOX defines interfaces to users codes through the abstract group and vector pure virtual classes. In this way, NOX allows the use of any linear algebra format (such as lapack, epetra, or petsc) and any linear solver compatible with the linear algebra format.

Models alternative to Newon-based techniques are also supplied in the form of lower order Broyden approximations and a higher-order Tensor-based model.

The native linear algebra format for Trilinos (and thus for NOX) is the epetra package. If using this format, NOX requires only the residual equaitons (optionally the user can provide the Jacobian, and preconditioner too). NOX’s epetra support can evaluate the Jacobian so users don’t have to code up derivatives of their nonlinear equation system. This includes brute force finite differencing, a much more efficient colored finite differencing, and a Jacobian-Free Newton-Krylov (JFNK) algorithm.

**LOCA**  
LOCA is generic continuation and bifurcation analysis package that is designed for large-scale applications. The algorithms are designed with minimal additional interface requirements over that needed for a Newton method to reach an equilibrium solution. LOCA is built upon the NOX nonlinear solver package. The algorithms in NOX are generic and written to the NOX::Abstract::Group and NOX::Abstract::Vector, which provide abstract interfaces to the linear algebra, data structures, and nonlinear equations to be solved. LOCA uses the NOX interface and extends it via additional abstract groups that provide the interface needed for continuation and bifurcation tracking, such as setting parameters and computing derivatives with respect to parameters.

LOCA provides several generic groups that take the NOX group representing the equilibrium equations and implement the extended sets of nonlinear equations representing various forms of continuation and bifurcations (such as the additional equation for arc-length continuation). These extended groups also include generic algorithms for computing the Newton step for the extended system based on the Newton step for the equilibrium equations (e.g. Sherman-Morrison-Woodbury formula for arclength continuation). They are then sent to NOX for nonlinear solution.

Finally, LOCA provides a stepper class that repeatedly calls the NOX nonlinear solver to compute points along a continuation curve. The design allows for continuation of bifurcations so two-parameter bifurcation sets can be generated. The stepper class relies on several support classes that compute predictors, step sizes, etc. Each of these are discussed in more detail in the LOCA Class Overview.

Unlike NOX which can provide a range of nonlinear solvers using a single abstract interface to the nonlinear equations and linear algebra, LOCA provides several different levels of functionality, each requiring additional information from the underlying problem. Therefore, the interface to LOCA is split among several abstract classes each encapsulating a different level of functionality. To interface to LOCA, the user need only provide implementations of those abstract classes for the functionality the user is interested in. LOCA provides two complete interfaces:

```
LOCA::LAPACK::Group  
LOCA::Epetra::Group    
```

both of which implement the interface required for all levels of functionality provided by LOCA.

## Documentation

NOX and LOCA is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [NOX and LOCA's Doxygen webpages](https://trilinos.github.io/docs/nox/index.html).


### Building NOX and LOCA Doxygen

The User's Guide and Developer's Guide for NOX and LOCA is combined on
a single website found at:

  https://trilinos.github.io/docs/nox/index.html

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


## Questions? 

Contact lead developers:

* **NOX team**        (GitHub handle: @trilinos/nox)
* **Roger Pawlowski** (GitHub handle: [rppawlo](https://github.com/rppawlo) or rppawlo@sandia.gov)
* **Eric Phipps**     (GitHub handle: [etphipp](https://github.com/etphipp) or etphipp@sandia.gov)

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For NOX-and-LOCA-specific copyright and license details, refer to the [nox/COPYRIGHT_NOX](COPYRIGHT_NOX) - [nox/COPYRIGHT_LOCA](COPYRIGHT_LOCA) and [nox/LICENSE_NOX](LICENSE_NOX) - [nox/LICENSE_LOCA](LICENSE_LOCA) files located in the `nox` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.
