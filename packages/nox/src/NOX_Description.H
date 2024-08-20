// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
  \mainpage NOX & LOCA: Object-Oriented Nonlinear Solver and Continuation Packages

  \section introduction Introduction

  %NOX is short for <em>Nonlinear Object-Oriented Solutions</em>, and
  its objective is to enable the robust and efficient solution of the
  equation: \f$ F(x)=0 \f$, where \f$F:\Re^n \rightarrow \Re^n\f$.
  %NOX implements a variety of Newton-based globalization techniques
  including %Line Search and %Trust Region algorithms.  In additon it
  provides higher and lower order models using Broyden and Tensor
  algorithms.  Special algorithms have been developed for use with
  inexact linear solvers such as Krylov subspace techniques. %NOX is
  designed to work with any linear algebra package and to be easily
  customized.  %NOX is part of Sandia's <A
  HREF="http://trilinos.sandia.gov/">Trilinos project</a>.

  %LOCA, distributed as part of %NOX,
  is short for <em>Library of Continuation Algorithms</em>, and its
  objective is to compute families of solutions to
  \f$
  F(x,p)=0
  \f$
  and their bifurcations, where
  \f$F:\Re^n\times\Re^m\rightarrow\Re^n\f$.  For %LOCA specific information,
  see \ref loca_overview.

  \section announcements Announcements

  NOX/LOCA version 10.2 will be released as part of Trilinos 10.2 in March 2010.  Please read the \ref nox_release_information page for more details on these changes.

  \section technical_information Technical Information

  <ul>
  <li> \ref nox_release_information - Latest news and notes on released software.
  <li> \ref nox_user_information - Instructions for building and linking %NOX and %LOCA to your appplication.
  <li> \ref nox_developer_information - Developer reference pages.
  <li> \ref portability_issues - Platform specific issues.
  </ul>

  \section license License

  %NOX is available for download under the terms
  of the <a href="http://www.gnu.org/copyleft/lesser.html">GNU Lesser
  General Public License</a>.

  \section downloads Downloads

  %NOX and %LOCA are packaged and released as part of the <a href="http://trilinos.sandia.gov/">Trilinos Solver Framework</a>.  You may download %NOX (which includes %LOCA) by downloading a copy of Trilinos.  Please see the <a href="http://trilinos.sandia.gov/download">Trilinos downloads page</a> for the latest releases. <br />

  \section bugs Reporting Bugs and Making Enhancement Requests

  To reports bugs or make enhancement requests, visit the <A HREF="http://software.sandia.gov/bugzilla/">Trilinos Bugzilla (Bug Tracking) Database</A>, and use the following instructions.
      <UL>
      <LI>Click on "Enter a new bug report"
      <LI>Choose "NOX"
      <LI>Either login or create a new account
      <LI>Submit your bug report
      </UL>

  \section email Mailing Lists

  <p>We recommend that users and developers subscribe to the following mailing lists as appropriate.
  <ul>
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-announce">NOX-Announce Mailing List</a> -
  Low-volume, primarily for release announcements. (Highly Recommended)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-users">NOX-Users Mailing List</a> -
  Discussion forum for users of %NOX. (Recommended)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-developers">NOX-Developers Mailing List</a> -
  Discussion forum for %NOX developers. (Developers)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-checkins">NOX-Checkins Mailing List</a> -
  CVS Commit Messages (Active  Developers)
  </ul>

  \section contributors NOX & LOCA Contributors

  The following people have contributed code to %NOX and/or %LOCA:

  <ul>
  <li> Brett Bader, Sandia National Labs, bwbader@sandia.gov
  <li> Russ Hooper, Sandia National Labs, rhoope@sandia.gov
  <li> <a href="http://csmr.ca.sandia.gov/~tgkolda/">Tammy Kolda</a>,
  Sandia National Labs,
  tgkolda@sandia.gov
  <li> Roger Pawlowski, Sandia National Labs, rppawlo@sandia.gov (%NOX Project Lead)
  <li> Eric Phipps, Sandia National Labs, etphipp@sandia.gov (Primary %LOCA developer)
  <li> Andy Salinger, Sandia National Labs, agsalin@sandia.gov (%LOCA Project Lead)
  <li> Alex Toth, NC State University
  </ul>

  \section questions For All Other Questions and Comments...

  Please contact Roger Pawlowski (rppawlo@sandia.gov) or Eric Phipps (etphipp@sandia.gov).

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_user_information NOX User Information

\section overview Overview

To get started with %NOX, you'll want to follow these four steps.
<ul>
<li> \ref step1
<li> \ref step2
<li> \ref step3
<li> \ref step4
</ul>

\section documentation  Generating Documentation

The documentation is available online at
<a href="http://trilinos.sandia.gov/packages/nox/">http://trilinos.sandia.gov/packages/nox/</a>.

To generate a local copy of the documentation, you'll need to have doxygen installed.
Then follow the instructions in \c nox/doc/README.

The documentation is generated in the source directory tree.  The main page is \c nox/doc/html/index.html.

\section Quicklinks

- \ref nox_configuration_options
- \ref parameters
- \ref downloads
- \ref bugs
- \ref epetra_interface
- \ref nox_epetra_tutorial
- \ref thyra_interface
- \ref petsc_interface
- \ref nox_class_overview

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step1 Step 1: Download, Configure, Compile, and Install NOX

\section step1_download Downloading and Unpacking NOX

%NOX is a package in the <a href="http://trilinos.sandia.gov">Trilinos Solver Framework</a>.  You will have to <a href="http://trilinos.sandia.gov/download">download</a> the Trilinos software, configure for your machine, build the libraries, and install them into your working directories.

The <a href="http://trilinos.sandia.gov/getting_started.html">Getting Started With Trilinos</a> page is an excellent resource for configuring and building Trilinos. Important links include:
<ul>
<li> <a href="http://trilinos.sandia.gov/Trilinos10CMakeQuickstart.txt">CMake Quickstart Guide</a>
</ul>

\section step1_problems If you have problems...

If you have problems configuring and compiling %Trilinos, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step1_moving_on Moving on...

Go on to \ref step2.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step2 Step 2: Create concrete implementations of the NOX::Abstract classes

\section step2_abstract NOX abstract classes

  %NOX's flexibility is based on the fact that it does not depend on
  any particular linear algebra package. Interfacing to %NOX requires
  concrete implementations or instatiations that derive from the
  following abstract classes:

  - NOX::Abstract::Vector
  - NOX::Abstract::Group

  The \b Vector supports basic vector operations such as dot products
  and so on. The \b Group supports the linear algebra functionality as
  well as the interface to evaluate the function (nonlinear equations)
  and, optionally, the Jacobian. Complete details are provided in the
  class descriptions.

  Users have two choices on how to supply the concrete implementations
  of the above classes: (1) write your own implementation of the
  NOX::Abstract::Vector and NOX::Abstract::Group classes (see \ref
  step2_writing) or (2) use one of the predefined implementations (see
  \ref step2_instantiations).

\section step2_writing Writing your own instantiation

  We recommend using the NOX::LAPACK concrete implementation as a
  guide.  In this case, the underlying vectors are C++ STL
  std::vector<double> objects, and the matrices are stored in our own
  NOX::LAPACK::Matrix class.  Both the vectors and matrices are
  manipulated using LAPACK.  The NOX::LAPACK::Vector is
  straightforward; see NOX_LAPACK_Vector.H and NOX_LAPACK_Vector.C.
  The NOX::LAPACK::Group uses the NOX::LAPACK::Vector and
  NOX::LAPACK::Matrix objects. The interface with the application is
  handled in a separate class (NOX::LAPACK::Interface) that is passed
  to the NOX::LAPACK::Group when it is constructed.

\section step2_instantiations Instantiations provided with NOX

  %NOX includes four ready-made instantiations.

  <ul>
  <li>
  The NOX::LAPACK instantiation is an interface to the BLAS/LAPACK
  library.  It is not intended for large-scale computations, but to
  serve as an easy-to-understand example of how one might interface to
  %NOX.  To compile the NOX::LAPACK library, use the CMake flag
  "NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON" option to enable.
  The NOX::LAPACK instatiation is serial only.
  <li>
  The NOX::Epetra instantiation is an interface to the Trilinos/Epetra
  linear algebra library.  It is Sandia's workhorse code for parallel
  applications.  For more information see the <a
  href="http://trilinos.sandia.gov/packages/epetra">Epetra Home
  Page</a>. The epetra support interface actually requires four
  Trilinos libraries: %Epetra, %EpetraExt, %AztecOO, and %Ifpack.  The
  NOX::Epetra library is automatically enabled if you configure
  Trilinos with the four required packages enabled.  Alternatively,
  you can force epetra support to be enabled with the CMake flag
  "NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA=ON".  The configure will
  fail if you do not have the required packages (%Epetra, %EpetraExt,
  %AztecOO, and %Ifpack) enabled.  For instructions on interfacing
  your code to the epetra implementation see \ref epetra_interface.
  Optionally, the epetra support can also use the %ML library.
  <li>
  The NOX::Thyra instantiation is an interface to the Trilinos/Thyra
  package.  The Thyra package contains a set of interfaces and
  supporting code that defines basic interoperability mechanisms
  between different types of numerical software. The foundation of all
  the interfaces related to abstract numerical algorithms (ANAs) are
  the mathematical concepts of vectors, vector spaces, and linear
  operators. All other ANA interfaces and support software are built
  on these fundamental operator/vector interfaces.  We expect this
  interface to supersceed the Epetra interface above when all
  functionality is moved into the Trilinos/Stratimikos layer.  For
  more information see the <a
  href="http://trilinos.sandia.gov/packages/thyra">Thyra Home
  Page</a>. Thyra will automatically be enabled if the Trilinos
  package Thyra is enabled. To explicitly force the NOX::Thyra library
  to be enabled, use the CMake option
  "NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA=ON" during configure.  For
  instructions on interfacing your code to the thyra implementation
  see \ref thyra_interface.
  <li>
  The NOX::PETSc instantiation is an interface with the PETSc
  library. PETSc is developed at <a href="http://www.anl.gov">Argonne
  National Labs</a>; for more information see the <a
  href="http://www-fp.mcs.anl.gov/petsc">PETSc Home Page</a> To
  compile the NOX::Petsc library and examples, use the CMake flag
  "NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC=ON".  For instructions on
  interfacing your code to the petsc implementation see \ref
  petsc_interface.
  </ul>

\section step2_problems If you have problems...

Please report any issues regarding %NOX and %LOCA to <A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step2_moving_on Moving on...

Go on to \ref step3.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step3 Step 3: Call NOX from your code

A small example based on the NOX::LAPACK Vector and Group is included
with %NOX; see the file
Trilinos/packages/nox/examples/lapack/NOX_SimpleExamples/Rosenbrock.C.

There are four steps to calling the solver.

\section step3_stepa Step A: Instantiate a Group object containing the initial guess.

The first step is to somehow construct an object deriving from the
Abstract::Group which contains the initial guess as its x-vector. See
\ref step2 for more details.  For example:

\code

// Construct an object that derives from the abstract Group object. Here
// we call it the ExampleGroup. It should contain the initial guess.
Teuchos::RCP<UserConcreteGroup> grp =
  Teuchos::rcp(new UserConcreteGroup);
\endcode

Note that the group is not handled as a raw pointer, but is wrapped in a reference counter smart pointer (RCP).  RCPs enforce strict memory management in %NOX and %LOCA and can prevent many errors that users have encountered with objects going out of scope.  We use the RCP implementation in the %Trilinos package <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a>.  For information on the use of smart pointers, see <a href="http://trilinos.sandia.gov/RefCountPtrBeginnersGuideSAND.pdf">Teuchos Reference Counted Pointer Beginners Guide</a>.

\section step3_stepb Step B: Construct the status tests

The NOX::StatusTest objects are used to determine when %NOX should
terminate the iterative process due to either convergence or
failure. For example...

- The NOX::StatusTest::MaxIters test is used to control the maximum
  number of iterations that the nonlinear solver is allowed to take.

- The NOX::StatusTest::NormF test checks the size of \f$\| F(x)
  \|\f$. This check can be done based on relative or absolute norms,
  scaled by \f$ \sqrt{n} \f$, and more. See the description for more
  details.

These and other status tests can be combined using the
NOX::StatusTest::Combo object. It takes an arbitrary number of status
tests and either AND's or OR's them together.

Alternatively, users can build a set of status tests from a Teuchos::ParameterList by using the NOX::StatusTest::Factory object.  See NOX::StatusTest::Factory for parameter list arguments and examples.

Finally, the user can create their own status tests, so long as they
derive from NOX::StatusTest::Generic.

Here is some sample code.

\code
// Set up the status tests
Teuchos::RCP<NOX::StatusTest::NormF> statusTestNormF =
  Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestMaxIters =
  Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo =
  Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                          statusTestNormF,
                                          statusTestMaxIters));
\endcode

This test uses an OR combination, meaning that either test registered with the Combo test can trigger termination.  The first test requires that the norm of the residual is less that 1.0e-4.  If met, %NOX will terminate returning a "Converged" result.  The second test will trigger %NOX to terminate if a maximum number of iterations is reached and return a "Failed" result.

\section step3_stepc Step C: Set the solver parameters

Next, we set up a parameter list containing the information on what
types of solvers and so forth should be used.

\code
// Create the list of solver parameters
Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
  Teuchos::rcp(new Teuchos::ParameterList);

// Select the solver (this is the default)
solverParametersPtr->set("Nonlinear Solver", "Line Search Based");

// Create the line search parameters sublist
Teuchos::ParameterList& lineSearchParameters = solverParametersPtr->sublist("Line Search");

// Set the line search method
lineSearchParameters.set("Method","More'-Thuente");
\endcode

For a full list of parameters; see \ref parameters.

\section step3_stepd Step D: Construct the solver and solve

The last step is to create the solver, passing in the group with the
initial guess, the status tests, and the solver parameters.

\code
// Create the solver
Teuchos::RCP<NOX::Solver::Generic> solver =
  NOX::Solver::buildSolver(grp, statusTestsCombo, solverParametersPtr);

// Solve the nonlinear system
NOX::StatusTest::StatusType status = solver->solve();

// Print the parameter list
cout << "\n" << "-- Parameter List From Solver --" << "\n";
solver->getList().print(cout);

// Get the answer
grp = solver->getSolutionGroup();

// Print the answer
cout << "\n" << "-- Final Solution From Solver --" << "\n";
grp.printSolution();
\endcode

\section step3_problems If you have problems...

Please report any issues regarding %NOX and %LOCA to <A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step3_moving_on Moving on...

Go on to \ref step4.

**/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step4 Step 4: Link your code to NOX

\section step4_export_makefiles Using the Export Makefile System

%NOX and %LOCA support the export makefile system to streamline issues dealing with linking against user applications.  When %NOX and %LOCA are installed using the "make install" command, it will install two Makefiles that contain all flags needed to link an application against the solvers.  The first Makefile is Makefile.export.NOX.  This contains the flags needed to build and link your code against nox.  The second file is Makefile.client.NOX.  This file can be used to link directly against your code.  See the files for more details.

\section step4_basics Compiling Without using Export Makefiles

Although we don't recommend it, users can link %NOX and %LOCA to their application without using the export makefile system.  To compile your code in this manner, your Makefile should look something like the following, assuming %NOX was installed in /usr/local.

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib
LIBS = -lnox
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

Compiling with LAPACK support:

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib
LIBS = -lnox -lnoxlapack -llapack -lblas
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

Compiling with Epetra support:

We assume that Epetra, EpetraExt, Aztec, and Ifpack are installed in the same directory as %NOX.

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib
LIBS = -lnox -lnoxepetra -laztecoo -lifpack -lepetraext -lepetra  -llapack -lblas
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

\note If compiling without using the export makefile system, you have to be very careful about which libraries are needed to link against.  Enabling configure options typically brings in new library requirements.  The export makefile system catches this automatically.  Not using export makefiles means YOU will have to catch the additional dependencies.

\section step4_problems If you have problems...

If you have problems linking %NOX and/or %LOCA, please report them using <A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_developer_information NOX Developer Information

To become a %NOX developer, contact Roger Pawlowski (rppawlo@sandia.gov).

Bug Tracking: <A HREF="http://software.sandia.gov/bugzilla/">Trilinos Bugzilla Database</A>.

CVS Respository View/Query Tool: <a
href="http://software.sandia.gov/bonsai/cvsqueryform.cgi?cvsroot=/space/CVS&module=nox">Trilinos Bonsai</a>

The following pages are relevant to developers.
<ul>
<li>\ref cvsrepos
<li>\ref coding
<li>\ref prerelease
</ul>


*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_release_information NOX and LOCA Release Information (Last updated March 1st 2010)


 %NOX/%LOCA version 10.2 will be release as part of Trilinos 10.2 in March 2010.
 <ul>
 <li> We have fixed a long standing issue with nox regarding overloading of the std::ostream operator.  If you are compiling nox and get a compile/link error that fails to find the "operator<<" method on one of your classes, then you need to move your overloaded std::ostream function for that object into the correct namespace (the same namespace as the object it is operating on).

The following artice by Herb Sutter explains this issue:

<a href="http://www.gotw.ca/publications/mill08.htm">Herb Sutter Article</a>

 </ul>


 %NOX/%LOCA version 10.0 will be release as part of Trilinos 10.0 in September 2009.
 <ul>
 <li> NOX has been changed over to use the cmake build system.  The autotools have been deprecated and will be completely removed in the next release.
 </ul>
 %NOX/%LOCA version 9.0 will be release as part of Trilinos 9.0 in September 2008.
 <ul>
 <li> No major changes to nox/loca.
 </ul>
 %NOX/%LOCA version 8.0 will be release as part of Trilinos 8.0 in early September 2007.  Current changes/additions:
<ul>
  <li> All %Manager classes have been replaced with Factories.  This change clarifies what occurs during reset calls and simplifies object management.  No more handles.  The factories return "Generic" base class objects.
  <li> A new concrete implementation of the group/vector has been added to support %Thyra objects.  %Thyra is a set of abstract interfaces designed to provide developers of abstract numerical algorithms a single point of contact with codes.  It plays the same role as the NOX::Abstract::Group and NOX::Abstract::Vector play in NOX. More information on %Thyra can be found at <A HREF="http://trilinos.sandia.gov/packages/thyra">Thyra</A>.  The group leverages the <A HREF="http://trilinos.sandia.gov/packages/stratimikos">Stratimikos library</A> to build the linear solver and preconditioner.  In the future, we plan to move away from the direct Epetra support class as the %Thyra wrappers gain more functionality.  Examples can be found in nox/test/epetra/Thyra_2Dsim directory.
  <li> A wrapper called Thyra::NOXNonlinearSolver has been added that derives from the Thyra::NonlinearSolverBase class.  This allows abstract numerical agorithms to use NOX as the nonlinear solver through a Thyra abstract interface.  In particular, NOX can now be used by Rythmos.
  <li> A NOX::StatusTest::Factory has been added to allow for the construction of an entire status test set via a parameter list.  Please use the nonmember functions to build parameter lists instead of using the Factory object directly.
  <li> The NOX::Utils object can now take a parameter sublist to enable the print options specified with the "Output Information" parameter.  Previously the user had to sum the print options into an integer that was set in the parameter list.  Now the user can just specify a sublist of print options for the "Output Information" parameter.  Setting the integer still works.  Here is an example:
\code
Teuchos::ParameterList p;
p.set("MyPID", myPID);
p.set("Output Processor", printProc);
p.set("Output Precision", 3);
p.set("Output Stream", outputstream);
p.set("Error Stream", outputstream);
Teuchos::ParameterList& info = p.sublist("Output Information");
info.set("Error", true);
info.set("Details", true);
info.set("Outer Iteration StatusTest", true);
info.set("Outer Iteration", true);
NOX::Utils utils3(p);
\endcode
  <li> Removed all instances of enumerated types from the parameter list.  This was causing problems with the python interfaces.  Enumerated type are now std::string types.
  <li> The RefCountPtr class in Teuchos was changed to RCP to shorten the amount of code.  All instances of RefCountPtr should be changed to RCP.  There is a script available in Teuchos to automatically make this change to any package.
</ul>

 %NOX/%LOCA version 7.0 will be release as part of Trilinos 7.0 in early September 2006.  %NOX and %LOCA have undergone many important changes between the 6.0 and 7.0 release:

<ul>
  <li> %NOX is now using reference counted smart pointers (RCPs) for memory management.  Passing an RCP wrapped object through a method call implies a persisting relationship.  We are using the Teuchos::RCP object for the RCP.
  <li> The NOX::Parameter::List has been removed.  We are now using the Teuchos::ParameterList object.  This has the same functionality with minor changes to the syntax of method calls.
  <li> %NOX now has a non-optional dependency on the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> package (also part of Trilinos).  While it can be built stand-alone, we recomend building %NOX as part of Trilinos with the teuchos library enabled to simplify dependencies.
  <li> %NOX and %LOCA now generate Makefiles that users can include in their own makefiles to get all dependencies correct.  For example if %NOX is built with Epetra support enabled, nox depends on epetra, aztecoo and ifpack headers and libraries.  To keep the user from having to learn about those packages so they can link nox into their codes, just include the %NOX file Makefile.export.nox file to get variables that have all sub package dependencies correct. In the installation directory, look at the files include/Makefiles.export.nox and include/Makefiles.export.nox.macros.  The important flags in Makefile.export.nox are NOX_INCLUDES, NOX_LIBS, LOCA_INCLUDES, and LOCA_LIBS flags.  The Makefile.export.nmox.macros file contains the compiler flags that were used to build the %NOX and %LOCA packages.  It is very important to use consistent flags when linking your application against Trilinos.
  <li> Multivector support is no longer a configure option, but is automatically built as part of the library.
  <li> The epetra interfaces (NOX::Epetra::Interface namespace) have been changed.  They now take an Epetra_Operator as an argument when evaluating the Jacobian and preconditioner objects.  While this is redundant for our base epetra group (since it shares the same operator between groups), this allows users to override the eptra group and use separate operators for each group if desired.
  <li> There are 2 options for creating an Ifpack preconditioners.  The original code where the "Preconditioner" option in the "Linear Solver" sublist was set to "Ifpack" and the new ifpack code "New Ifpack".  The original option explicitly constructs an ifpack ilu object that has since been deprecated and will most likely be removed in the next release (Trilinos 8.0).  The option "New Ifpack" uses the new ifpack factory object.  You can now create a parameter list that directly uses the ifpack factory arguments and can create any ifpack preconditioner (not limited to ILU).  We recommend switching to the "New Ifpack" preconditionerr interface asap.
  <li> In the epetra support, the preconditioner options "AztecOO", "New Ifpack", and "ML" allow users to recompute the preconditioner using the same graph so that preconditioner memory is not allocated and destroyed each linear solve.  The new "Linear Solver" sublist key is "Preconditioner Reuse Policy" with the possible values "Rebuild", "Reuse", and "Recompute".  If you choose "Rebuild", this is what nox has done in previous releases - destroy, reallocate and recompute the preconditioner for each linear solve.   "Reuse" wil use the same preconditioner for every linear solve.  "Recompute" will recompute the preconditioner, but will try to efficiently reuse any objects that don't need to be destroyed.  How efficient the "Recompute" option is depends on the type of preconditioner.  For example if we are using ILU in Ifpack, the graph that is allocated does not need to be destroyed and recomputed even though the computation will recompute the factorization when the matrix changes.  These changes support native Aztec, Ifpack and ML preconditioners.
  <li> Build options have changed.  please run "configure --help" to see the new options.
  <li> Many of the internals of %LOCA have been rewritten to more fully support multi-parameter continuation, constraint tracking, higher-order predictors and advanced bifurcation algorithms.  The old framework has been completely removed, and therefore all users must upgrade their interface to the new framework.  For a brief description of the new framework, see the \ref new_loca_framework page.  In addition, some of the parameter list entries that control %LOCA have changed.  See the \ref loca_parameters for a description of the new parameter list structure.
  <li> A minimally augmented turning point formulation has been added to %LOCA.  This method appears to be significantly faster and more robust than the turning point methods based on the Moore-Spence formulation.
  <li> A modified pitchfork bordering method has been added to %LOCA, which is similar to the modified turning point method already available.  This method avoids the singular matrix solves in the original pitchfork bordering method.
  <li> Parallel and serial Hopf support has been added to %LOCA.
</ul>

NOX/LOCA 4.1 is being released as part of Trilinos 6.0 on September 1st 2005.

General release 4.1 notes:
<ul>
<li> %NOX and %LOCA now use the new export makefile system.  When compiled and installed into the prefix directory, two new makefiles will appear: Makefile.export.nox and Makefile.export.nox.macros.  These files contain all variables needed my a simple Makefile to build and link against nox.  All major packages in Trilinos now support this feature.
<li> %NOX and %LOCA's examples now work in a consistent manner with the tests.  The examples and tests can be enabled or disabled with the flags --enable-nox-tests and --enable-nox-examples.  The specific examples that will be compiled depend on which support packages are enabled (lapack, epetra, petsc, belos, tsf). You no longer need to specify specific example support like --enable-nox-epetra-examples.
</ul>

%LOCA release 4.1 notes:
<ul>
<li> A constraint tracking feature has been added to %LOCA where the user can
specify additional algebraic constraints that must be satisfied along with
the residual.  Examples include boundary conditions, phase conditions in
periodic orbits, etc.  The resulting Newton system can be efficiently solved
using the QR method discussed below.  This capability is only available using
the LOCA::NewStepper also discussed below.
<li> A QR method for solving augmented systems of equations in parallel has
been added, and is an extension of the Householder method already available
to systems with more than one augmented row/column.  This method allows
bordered systems of equations to be solved in roughly the same time as
original system and is well-conditioned even when the original system is
nearly singular.  This capability is only available using the LOCA::NewStepper
discussed below.
<li> Many of the internals of %LOCA have been rewritten to more fully support
multi-parameter continuation, constraint tracking, higher-order predictors
and advanced bifurcation algorithms.  This process is mostly, but not fully
complete.  This functionality is available by choosing the LOCA::NewStepper
for continuation instead of the original LOCA::Stepper.  However pitchfork
and Hopf continuation are currently not available in the new framework.
Documentation for this new framework can be found on the
\ref new_loca_framework page.
<li> LOCA's interface to eigensolvers has changed slightly and now uses
a strategy interface approach instead of relying on the group to provide
computeEigenvalues() and saveEigenvalues() methods.  This provides a more
natural decoupling of the parallel linear algebra data structures (i.e., the
group) from the eigensolver method/package (eigensolver strategy object).
The use of strategy interfaces instead of group methods is a major theme in
the new %LOCA framework, and is discussed more fully on the new framework page.
<li> Teuchos::RCP support has been added to improve software quality.
For the most part, ref-count pointers are only used in the new framework.
<li> Many regression tests have been added to improve software quality.
<li> Libtool support has been removed from %LOCA due to portability issues.
<li> Python support is temporarily disabled in this release due to a major
refactor of the build system for python wrappers in Trilinos.
</ul>

%LOCA portability notes:
<ul>
<li>
There currently is a portability issue with %LOCA to the Solaris platform using the Sun Forte 7 compilers.  %LOCA builds on this platform but does not run correctly.  It appears the problem is compiler related as the virtual table for some %LOCA objects is not setup correctly.  All of the %LOCA examples and any %LOCA tests involving continuation runs fail with a segmentation fault.  It is not known whether %LOCA runs correctly when built with version 10 of the Sun compilers.  See bugzilla bug <a href="http://software.sandia.gov/bugzilla/show_bug.cgi?id=1238"> 1238 </a> for more information on this issue.
</ul>

NOX release 4.0 notes:
<ul>
<li> A new Tensor based solver that is independent of the linear solver is now implemented (in the prerelease section).
<li> Inexact Trust Region Techniques are fully supported (also in the prerelease section).
<li> Improved support for Epetra objects - The NOX::EpetraNew::Group should be used in place of the NOX::Epetra::Group.  In the future we will be phasing out the NOX::Epetra::Group and replace it with NOX::EpetraNew::Group.
<li> A multi-physics technique for code coupling has been demonstrated in the test/epetra/Multiphysics directory.  This will be expanded to a generic capability in the next year.
</ul>

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_configuration_options NOX Configuration Options

Here is a list and explanation of %NOX specific configuration options.  See the <a href="http://trilinos.sandia.gov/TrilinosUserGuide.pdf">Trilinos User Guide</a> for general options supported by all packages in %Trilinos.

To get a complete list of all configure options supported by %NOX and %LOCA do the following:
\code
    $ cd $BUILD_DIR
    $ rm CMakeCache.txt
    $ cmake -LAH -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=ON \
      $TRILINOS_HOME

    NOTE: You can also just look at the text file CMakeCache.txt after
    configure which gets created in the build directory and has all of the
    cache variables and documentation.
\endcode

\section basic_config Basic Options

<ul>
  <li> \c NOX_ENABLE_DEBUG=ON
    <ul>
      <li> This turns on special debugging checks and output in %NOX.
    </ul>
  <li> \c NOX_ENABLE_EXAMPLES=ON
    <ul>
      <li> Builds the example code in %NOX and %LOCA. Overrides the Trilinos_ENABLE_EXAMPLES configure flag for the nox package only.
    </ul>
  <li> \c NOX_ENABLE_TESTS=ON
    <ul>
      <li> Builds the tests in %NOX and %LOCA. Overrides the Trilinos_ENABLE_TESTS configure flag for the nox package only.
    </ul>
</ul>

\section lapack_options LAPACK-Related Options

<ul>
  <li> \c  NOX_ENABLE_ABSTRACT_IMPLEMENTATION_LAPACK=ON
    <ul>
      <li> Compile noxlapack library (libnoxlapack.a). Default is not to compile the library.
    </ul>
</ul>

\section epetra_options Epetra-Related Options

<ul>
  <li> \c NOX_ENABLE_ABSTRACT_IMPLEMENTATION_EPETRA=ON
    <ul>
      <li> Compile noxepetra library (libnoxepetra.a). Default is not to compile the library unless Epetra, EpetraExt, Ifpack, and AztecOO are enabled.
    </ul>
</ul>

\section thyra_options Thyra-Related Options

<ul>
  <li> \c NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA=ON
    <ul>
      <li> Compile noxthyra library (libnoxthyra.a). Default is not to compile the library unless the Thyra package is enabled.
    </ul>
</ul>

\section petsc_options PETSc-Related Options

<ul>
  <li> \c NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC=ON
    <ul>
      <li> Compile noxpetsc library (libnoxpetsc.a). Default is not to compile the library.
      <li> Requires existing PETSc source specified by environment variables,
           \c PETSC_DIR and \c PETSC_ARCH
    </ul>

  See \ref petsc_interface for more information regarding building
       %NOX with PETSc
</ul>

\section loca_config LOCA-Related Options

<ul>
  <li> \c NOX_ENABLE_LOCA=ON
    <ul>
      <li> Compile the loca library (libloca.a). Default is to always compile the library.
      <li> If %LAPACK support is enabled, Trilinos must be configured with Teuchos_ENABLE_COMPLEX for required complex support.
    </ul>
  <li> \c NOX_ENABLE_Anasazi=ON
    <ul>
      <li> Enable <a href="http://trilinos.sandia.gov/packages/anasazi/index.html"> Anasazi </a> eigensolver support in %LOCA.
      <li> %LOCA has a generic interface to Anasazi allowing it to compute eigenvalues/eigenvectors for any group/vector implementation.
      <li> Automatically enabled if the Anasazi library is enabled.
    </ul>
  <li> \c NOX_ENABLE_MF=ON
    <ul>
      <li> Enable Multifario support in %LOCA.
      <li> <a href="http://www.coin-or.org/multifario/"> Multifario </a> (MF) is a multi-parameter continuation library written by Mike Henderson at IBM.  Including this option allows %LOCA to use MF for multi-parameter continuation.
      <li> The version of MF that is known to work with %LOCA is available in the Trilinos3PL.
      <li> The multifario libraries can be specified using the TPL support.  See the Trilinos CMake quickstart guide for more details on enabling third party libraries.
    </ul>
</ul>

\section developer_config Developer-Related Options

<ul>
<li>  \c NOX_BUILD_PRERELEASE=ON
  <ul>
    <li> Cannot be used with the distribution version of the code - only works with the development version.
         Compiles code that is not included in the releases but is part of the CVS repository.
         See \ref prerelease.
  </ul>
</ul>

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page prerelease NOX Developer's Guide to Prerelease Code


We define <em>prerelease code</em> to be code that should not be
distributed, documented, or tested as part of the automatic nightly
scripts.

Prerelease code will only be used if configure is executed with the \c
NOX_BUILD_PRERELEASE=ON option. In that case,
<ul>
<li>all prerelease code will be compiled (\c -DWITH_PRERELEASE is added to \c CPPFLAGS),
<li>doxygen will document prerelease code (\c WITH_PRERELEASE is added to the \c PREDEFINED variable in Doxyfile),
<li>\c BUILD_PRERELEASE is true for automake/autoconf Makefile generation.
</ul>

There are two steps to adding prerelease code.

<ol>

<li> Add ifdef's to the source and header files.

In other words, the code should be surrounded by a \c WITH_PRERELEASE
ifdef as follows.

\verbatim
#ifdef WITH_PRERELEASE
// Insert All Code Here
#endif
\endverbatim

<li> Modify the CMakelists.txt file.

For example, the following modifications were made to \c
nox/src/Makefile.am to add NOX_Direction_QuasiNewton as prerelease code
in the libnox.a library.

\verbatim
# Prerelease Code
IF (${PACKAGE_NAME}_BUILD_PRERELEASE)

  APPEND_SET(HEADERS
    NOX_Solver_TensorBasedTest.H
    NOX_LineSearch_Tensor.H
    NOX_Direction_Tensor.H
    NOX_Direction_QuasiNewton.H
    NOX_Direction_ModifiedNewton.H
    )

  APPEND_SET(SOURCES
    NOX_Solver_TensorBasedTest.C
    NOX_LineSearch_Tensor.C
    NOX_Direction_Tensor.C
    NOX_Direction_QuasiNewton.C
    NOX_Direction_ModifiedNewton.C
    )

ENDIF()
\endverbatim

</ol>

\note Prerelease code must still be able to compile.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page cvsrepos NOX Developer's Guide to the Git Repository for NOX

The git repository is located on software.sandia.gov. See Bill Goldman (bgoldma@sandia.gov) for an account.

You should not do any development on software.sandia.gov. Instead,
check out a copy to your local machine.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page coding NOX Developer's Coding Guidelines

\htmlonly
Adapted largely from <a
href="http://www.doc.ic.ac.uk/lab/cplus/c++.rules/">Programming in
C++, Rules and Recommendations, by Mats Henricson and Erik
Nyquist</a>.


<H3>Structure of the Files</H3>

<H5>Naming Conventions</H5>
<ul>
<li>C++ header files end in <code>.H</code> and source files end in
<code>.C</code>
<br>&nbsp;
<li>The name of the files should correspond to the name of the class
they define, with double-colons replaced by underscores. For example,
the definition of the class <code>NOX::Abstract::Group</code> is in
the file <code>NOX_Abstract_Group.H</code>.
<br>&nbsp;
</ul>

<H5>General File Structure</H5>
<ul>
<li>Each file should begin as follows:

<pre>
// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
</pre>

Once the file is committed to the CVS repository, the CVS Information lines
will look something like the following:

<pre>
</pre>

The header information is automatically filled in between the two
<code>//@HEADER</code> keys when we run the <code>nox/maintenance/update_nox_headers.sh</code> command.

</ul>


<H5>Include File Structure</H5>
<ul>
<li>No include file should define more than one class.
<ul>
<li> The exception is for utility classes that are <em>only</em> used by the main class defined in the header file.
</ul>
<br>&nbsp;
<li>Every include file must contain a mechanism that prevents multiple
inclusions of the file. For example, the following should follow the
header information for the <code>NOX_Abstract_Vector.H</code> header
file.

<pre>
#ifndef NOX_ABSTRACT_VECTOR_H
#define NOX_ABSTRACT_VECTOR_H

<em>...body of include file goes here...</em>

#endif
</pre>

<li>Do not include system files (e.g., <CODE>iostream</CODE>) directly
in your files. Instead, include <CODE>NOX_Common.H</CODE>. The goal is
to better enable system portability since some machines have
<CODE>math.h</CODE> and others have <CODE>cmath</CODE> and so
on. Currently, we have the following system headers:
<br>&nbsp;
<ul>
<li><code>cstdlib</code>
<li><code>cstdio</code>
<li><code>cassert</code>
<li><code>iostream</code>
<li><code>iomanip</code>
<li><code>string</code>
<li><code>cmath</code>
<li><code>vector</code>
<li><code>map</code>
<li><code>deque</code>
<li><code>algorithm</code>
<li><code>sstream</code>
<li><code>fstream</code>
<br>&nbsp;
</ul>
<li>Do not use relative or absolute paths when including header files.  The autotools will set the paths correctly for you.  This is especially important since installing the headers and libraries moves the headers into a flat directory structure.
<ul>
<li> DO NOT do this:<br>
     <code>#include "../test/utils/NOX_TestCompare.H"</code>
<li> Do this: <br>
     <code>#include "NOX_TestCompare.H"</code>
</ul><br>&nbsp;
<li>Definitions of classes that are only accessed via pointers
(<CODE>*</CODE>) or references (<CODE>&amp;</CODE>) should be declared
using forward declarations, and <EM>not</EM> by including the header files.
<br>&nbsp;
<li>These are the cases when header files should be included in the
header file:
<br>&nbsp;
<UL>
<LI>classes that are used as <EM>base classes</EM>,
<LI>classes that are used as <EM>member variables</EM>,
<LI>classes that appear as <EM>return types</EM> or as <EM>argument
types</EM> in function/member function prototypes.
</UL>
</ul>

<H3>Naming Conventions</H3>
<ul>
<li>Everything must be declared within the <code>NOX</code>
namespace. No exceptions!
<br>&nbsp;
<li>Furthermore, each class should be within the appropriate
sub-namespace. The choices are:
<br>&nbsp;
<ul>
<li><code>Abstract</code>
<li><code>Parameter</code>
<li><code>Solver</code>
<li><code>StatusTest</code>
<li><code>LineSearch</code>
<li><code>Direction</code>
<li><code>Epetra</code>
<br>&nbsp;
</ul>
<li>Class names should begin with an uppercase letter. Variable and
function names should begin with a lowercase letter.
<br>&nbsp;
<li>In names (function, class, variable, etc)
    which consist of more than one word, the words are
    written together and each word that follows the first is begun
    with an uppercase letter. (e.g., <code>NOX::Linesearch::MoreThuente</code>).
<br>&nbsp;
<li>Do not use identifiers which begin
with one or two underscores (`<CODE>_</CODE>' or `<CODE>__</CODE>').
<br>&nbsp;
<li>Names should not include abbreviations that are not generally
    accepted.
<br>&nbsp;
<li>Choose variable names that suggest the usage.
</ul>

<H3>Style</H3>

<h5>Classes</h5>
<ul>
<li>The public, protected, and private sections of a class are to be
    declared in that order (the public section is declared before the
    protected section which is declared before the private section).
<br>&nbsp;
<li>No inline functions, except the empty <code>{}</code>
function. The reason for this is that if the function does not inline
correctly, it can actually lead to slower code rather than faster
code.
</ul>
<h5>Functions</h5>
<ul>
<li>Always provide the
<EM>return type</EM> of a function explicitly.
<br>&nbsp;


<li>When declaring functions, the leading parenthesis and the first
argument (if any) are to be written on the <EM>same line</EM> as the
function name. If space permits, other arguments and the closing
parenthesis may also be written on the same line as the function
name. Otherwise, each additional argument is to be written on a
separate line (with the closing parenthesis directly after the last
argument).
<br>&nbsp;

<li>Always write the left parenthesis
directly after a function name.
<PRE>
   void foo ();    // No!!
   void foo();     // Better
</PRE>
</ul>

<h5>Variable declarations</h5>

<ul>
<li>Only one variable per line.
<PRE>
   int i,j;   // No!!

   int i;     // Yes
   int j;
</PRE>

<li>
The characters `<CODE>*</CODE>' and `<CODE>&amp;</CODE>' should
be written together with the types of variables instead of with the
names of variables in order to emphasize that they are part of the type
definition. Instead of saying that <CODE>*i</CODE> is an <CODE>int</CODE>,
say that <CODE>i</CODE> is an <CODE>int*</CODE>.
<PRE>
   int *i;   // No!!
   int* i;   // Yes
</PRE>
</ul>

<h5>Loops and conditionals: <code>if</code>, <code>for</code>,
<code>while</code>, etc.</h5>

<ul>

<li>User parens to make code readable.
<pre>
  if (a == b && c < d || e == f)  // No!
  {
    /* Stuff */
  }

  if (((a == b) && (c < d)) || (e == f))  // Yes
  {
    /* Stuff */
  }
</pre>
<li>The block of any <code>if</code> statement should always follow on
a separate line.
<PRE>
   if ( /*Something*/ ) i++; // No!!

   if ( /*Something*/ )      // Yes!
     i++;
</PRE>
<li>Braces ("{}") which enclose a block should be aligned as follows:
<PRE>
   if ( /*Something*/ ) // Yes!
   {
     i++;
     j++;
   }

   if ( /*Something*/ ) { // Okay
     i++;
     j++;
   }

   if ( /*Something*/ ) // No!
     {
       i++;
       j++;
     }
</PRE>

Adding the following line to your \c .emacs file will help:

<PRE>
(c-set-offset 'substatement-open 0)
</PRE>

</ul>


<h5>Miscellaneous</h5>
<UL>
<li>Always provide a space on both sides of <code>=</code> signs and all logical
operators.
<br>&nbsp;
<li>Each statement shall always be in a separate line, however small it may appear.
<br>&nbsp;

<Li>Do not use <EM>spaces</EM> around
`<CODE>.</CODE>' or `<CODE>-&gt;</CODE>', nor between unary operators
and operands.
<br>&nbsp;
<li>Use the c++ mode in GNU Emacs to
format code.
</UL>

<H3>Coding Rules</H3>


<ul>
<li>A public member function must never return a non-const reference
or pointer to member data.
<br>&nbsp;
<li>Constants are to be defined using <code>const</code> or
<code>enum</code>; never using <code>#define</code>.
<br>&nbsp;
<li>A switch statement must always contain a default branch which
handles unexpected cases.
</ul>


<H3>Output</H3>

The <code>NOX::Utils</code> class has utility functions related
to printing. To use it, include <code>NOX_Utils.H</code>.

<ul>
<li>For <b>any</b> non-error print statement, call the
<code>NOX::Utils::out</code> or <code>NOX::Utils::pout</code> functions with the appropriate
MsgType flag. The flags are:
<br>&nbsp;
<ul>
<li><code>NOX::Utils::Error  </code>
<li><code>NOX::Utils::Warning</code>
<li><code>NOX::Utils::OuterIteration</code>
<li><code>NOX::Utils::InnerIteration</code>
<li><code>NOX::Utils::Parameters</code>
<li><code>NOX::Utils::Details</code>
<li><code>NOX::Utils::OuterIterationStatusTest</code>
<li><code>NOX::Utils::LinearSolverDetails</code>
<li><code>NOX::Utils::TestDetails</code>
<li><code>NOX::Utils::StepperIteration</code>
<li><code>NOX::Utils::StepperDetails</code>
<li><code>NOX::Utils::StepperParameters</code>
<li><code>NOX::Utils::Debug</code>
</ul>

</ul>


<H3>Error Handling</H3>

<ul>
<li>Always check return values of functions for errors.
<li>In general, try to recover from errors.
<li>If you must throw an exception, always print an explanation with
the function name to <code>NOX::Utils::err</code> or <code>NOX::Utils::perr</code> and then throw an exception with the std::string
<code>"NOX Error"</code>. For example if utils is a NOX::Utils object,

<pre>
  if (/* Error Condition */)
  {
    utils.err() << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << std::endl;
    throw std::runtime_error("NOX Error");
  }
</pre>

</ul>

<H3>Comments</H3>

We use Doxygen for the comments. To generate the documentation.  Instructions can be found in the README file in the top-level nox directory.

<ul>
<li> Document each class, function, and enum in the header
files.
<br>&nbsp;
<ul>
<li> The one exception is that functions in derived objects do not
need to be documented <b>if</b> the documentation is inherited from
the base class. This should be tested in Doxygen to be sure that it
works correctly.
<br>&nbsp;
</ul>
<li> Here's an example of documented a class. Note the formatting of
the comments. It's a C-style comment. The open comment marker
(<code>/*</code>) is followed by an exclamation mark to indicate that
it's a Doxygen comment. The open and close comment markers are on
lines by themselves, and the text of the comment is indented two
spaces. Always include a <code>\brief</code> description. The long
description follows. Observe the use of the formatting tags
<code>\c</code> and <code>\e</code>. The <code>\note</code> tag is
used for any special notes. The <code>\author</code> tag is
recommended.

<pre>
/*!
  \brief Arbitrary combination of status tests.

  In the \c AND (see NOX::Status::Combo::ComboType) combination, the
  result is \c Unconverged (see NOX::Status::StatusType) if \e any of
  the tests is \c Unconverged. Otherwise, the result is equal to the
  result of the \e first test in the list that is either \c Converged
  or \c Failed. It is not recommended to mix \c Converged and \c
  Failed tests in an \c AND combination.

  In the \c OR combination, the result is \c Unconverged if \e all of
  the tests are \c Unconverged. Otherwise, it is the result of the \e
  first test in the list that is either \c Converged or \c
  Failed. Therefore, it will generally make sense to put the \c Failed
  -type tests at the end of the \c OR list.

  \note We always runs through all tests, even if we don't need
  to. This is useful so that the user knows which tests have and have
  not be satisfied.

  \author Tammy Kolda (SNL 8950)
*/
class Combo : public Test {
...
}; // class Combo
</pre>

<li>Any parameters that are used within the class <b>must</b> be
documented in the class description <b>and</b> in the file NOX_Description.H on the parameters "page".
Note that the name, a brief description, and the default value for
each parameter is listed.
<pre>
/*!
  \brief %Newton-like solver with a line search.

  The following parameters are valid for this solver:

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::Linesearch::Factory. Defaults to an empty list.

  - "Linear Solver" - Sublist of the linear solver parameters, passed
    to Abstract::Group::computeNewton(). Furthermore, the "Tolerance"
    within this list may be modified by the
    resetForcingTerm(). Defaults to an empty list.

  - "Forcing Term Method" - Method to compute the forcing term, i.e.,
    the tolerance for the linear solver. Defaults to ""
    (nothing). Choices are "Type 1" and "Type 2".

  - "Forcing Term Minimum Tolerance" - Minimum acceptable linear
    solver tolerance. Defaults to 1.0e-6.

  - "Forcing Term Maximum Tolerance" = Maximum acceptable linear
    solver tolerance. Default to 0.01.

  - "Forcing Term Alpha" - Used for the "Type 2" forcing term
    calculation. Defaults to 1.5.

  - "Forcing Term Gamma" - Used for the "Type 2" forcing term
    calculation. Defaults to 0.9.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/
</pre>
Here's a more complicated example to produce a two-tiered list.
<pre>
/*!
   The parameters must specify the type of line search as well as all
   the corresponding parameters for that line search.

   &lt;ul&gt;
   &lt;li&gt; "Method" - Name of the line search. Valid choices are
     &lt;ul&gt;
     &lt;li&gt; "Full Step" (NOX::Linesearch::FullStep)
     &lt;li&gt; "Interval %Halving" (NOX::Linesearch::Halving)
     &lt;li&gt; "Polynomial" (NOX::Linesearch::Polynomial)
     &lt;li&gt; "More'-Thuente" (NOX::Linesearch::MoreThuente)
     &lt;/ul&gt;
   &lt;/ul&gt;
 */
</pre>


<li>Constants and enums can generally be described with simple
<code>\brief</code> comments. Those can be formatted in either of two
ways, as follows.

<pre>
  /*!
    \brief The test can be either the AND of all the component
    tests, or the OR of all the component tests.
  */
  enum ComboType {AND, OR};

  //! Constructor
  Combo(ComboType t = OR);
</pre>

<li>Doxygen does automatically cross-linking, which is very
convenient. However, sometimes it cross-links when you don't intend
for it to. For example, the following line would automatically
generate a link from the word <code>Newton</code> to the
NOX::Solver::Newton class.

<pre>
//! Newton-like solver with a line search.
</pre>

To prevent that automatic link, insert a percent sign (<code>%</code>)
immediately before the word that is causing the link. For example,

<pre>
//! %Newton-like solver with a line search.
</pre>
</ul>
\endhtmlonly

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_class_overview NOX Class Overview

  %NOX's flexibility is based on the fact that it does not depend on
  any particular linear algebra package. In order to interface to %NOX,
  the user needs to supply methods that derive from the following
  abstract classes:

  - NOX::Abstract::Vector
  - NOX::Abstract::Group

  The \b Vector supports basic vector operations such as dot products
  and so on. The \b Group supports the linear algebra functionality as
  well as the interface to evaluate the function and,
  optionally, the Jacobian. Complete details are provided in
  the class description.

  Included with %NOX are four such interfaces:

  - NOX::LAPACK
  - NOX::Epetra
  - NOX::Thyra
  - NOX::Petsc

  The \b LAPACK interface is an interface to the BLAS/LAPACK library.
  It is not intended for large-scale computations, but to serve as an
  easy-to-understand example of how one might interface to %NOX.

  The \b Epetra interface is an interface to the
  Trilinos/Epetra library developed by
  <a href="http://www.sandia.gov">Sandia National Labs</a>;
  for more information see the
  <a href="http://trilinos.sandia.gov/trilinos">Trilinos Home Page</a>.

  The \b Thyra interface is an interface to the
  Trilinos/Thyra library developed by
  <a href="http://www.sandia.gov">Sandia National Labs</a>;
  for more information see the
  <a href="http://trilinos.sandia.gov/trilinos">Trilinos Home Page</a>.

  The \b PETSc interface is an interface with the PETSc library. PETSc
  was developed at
  <a href="http://www.anl.gov">Argonne National Labs</a>;
  for more information see the
  <a href="http://www-fp.mcs.anl.gov/petsc">PETSc Home Page</a>

  All solvers are in the NOX::Solver namespace. The solvers are
  constructed via the NOX::Solver::Factory. The recommended solver is
  NOX::Solver::LineSearchBased, which is a basic nonlinear solver
  based on a line search. Each solver has a number of options that can
  be specified, as documented in each class or on the \ref parameters.

  The search directions are in the NOX::Direction namespace and
  accessed via the NOX::Direction::Factory. The default search
  direction for a line-search based method in NOX::Direction::Newton.

  Several line searches are available, as defined in the
  NOX::LineSearch namespace, and accessed via the NOX::LineSearch::Factory
  class. Examples include

  - NOX::LineSearch::FullStep
  - NOX::LineSearch::Backtrack
  - NOX::LineSearch::Polynomial
  - NOX::LineSearch::MoreThuente

  Convergence or failure of a given solver method is determined by
  the status tests defined in the NOX::StatusTest namespace. Various
  status tests may be combined via the NOX::StatusTest::Combo
  object. Users are free to create additional status tests that derive
  from the NOX::StatusTest::Generic class.


    */

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page parameters NOX Parameter Reference Page

  All parameters are specified via parameter lists as defined by the
  Teuchos::ParameterList class.

  Here is the basic hierarchical structure of the top level parameter
  list that is passed to a NOX::Solver::Factory object or the builder
  method.

  \note These parameters are valid for the general %NOX library.  The
  LAPACK, Epetra, Thyra, and PETSc Group/Vector support objects have
  their own set of parameters that are documented on each object's
  doxygen page.

<ul>
<li> "Nonlinear Solver" [std::string] - see NOX::Solver::Factory for complete list of options. Valid choices are:
   <ul>
     <li> "Line Search Based" or "Newton" (NOX::Solver::LineSearchBased) [Default]
     <li> "Trust Region Based" (NOX::Solver::TrustRegionBased)
     <li> "Inexact Trust Region Based" (NOX::Solver::InexactTrustRegionBased)
     <li> "Tensor Based" (NOX::Solver::TensorBased)
   </ul>


<li> "Direction" [sublist] - used by NOX::Solver derived objects to compute the primary direction associated with the method:

     <ul>

     <li> "Method" [std::string] - see NOX::Direction::Factory for complete list of options.  Valid options are:
           <ul>
       <li>
           <li> "Newton" NOX::Direction::Newton
       <li> "Steepest Descent" NOX::Direction::SteepestDescent
       <li> "NonlinearCG" NOX::Direction::NonlinearCG
       <li> "Broyden" NOX::Direction::Broyden
       <li> "Tensor" PRERELEASE ONLY! NOX::Direction::Tensor
       <li> "Modified-Newton" PRERELEASE ONLY! NOX::Direction::ModifiedNewton
       <li> "Quasi-Newton" PRERELEASE ONLY! NOX::Direction::QuasiNewton
       <li> "User Defined" - see NOX::Direction::Factory for details
       </ul>

     <li> "User Defined Constructor" [Teuchos::RCP<NOX::Direction::UserDefinedFactory>] - see NOX::Direction::Factory for details.

     <li> "Newton" [sublist] - used by NOX::Direction::Newton,
          see class documentation for details.

      <ul>

      <li> "Forcing Term Method" [std::string]
      <li> "Forcing Term Initial Tolerance" [double]
      <li> "Forcing Term Minimum Tolerance" [double]
      <li> "Forcing Term Maximum Tolerance" [double]
      <li> "Forcing Term Alpha" [double]
      <li> "Forcing Term Gamma" [double]
      <li> "Rescue Bad Newton Solver" [boolean]

      <li> "Linear Solver" [sublist] - passed to
           NOX::Abstract::Group::computeNewton() by
           NOX::Direction::Newton.

           <ul>
           <li> "Tolerance" [double]
           </ul>

           NOTE: This list will carry any number of linear solver
           options depending upon the concrete implementation of
           the linear solver used in the applyJacobianInverse()
           method of the NOX::Abstract::Group class.  For Epetra
           support, this option list can be found in the
           NOX::Epetra::LinearSystemAztecOO object.

      </ul>

     <li> "Steepest Descent" [sublist] - used by
          NOX::Direction::SteepestDescent, see class documentation for
          details

      <ul>
      <li> "Scaling Type" [std::string]  Type of scaling to use.  Options are:
               <ul>
           <li> "2-Norm"
           <li> "Quadratic Model Min"
           <li> "F 2-Norm"
           <li> "None"
           </ul>
      </ul>

     </ul>

<li> "Cauchy Direction" [sublist] - used to specify teh second direction in dogleg trusregions methods.  See NOX::Solver::TrustRegionBased and NOX::Solver::InexactTrustRegionBased solvers for mor einformation on how this direction is used.  This sublist uses the same options as "Direction" sublist.

<li> "Line Search" [sublist] - used by NOX::Solver::LineSearchBased solvers to specify the line search algorithm to use.
     <ul>
     <li> "Method" [std::string] - see NOX::LineSearch::Factory for a detailed list of options.  Valid options are:
          <ul>
      <li> "Full Step" (NOX::LineSearch::FullStep)
      <li> "Backtrack" (NOX::LineSearch::Backtrack)
      <li> "Polynomial" (NOX::LineSearch::Polynomial)
      <li> "More'-Thuente" (NOX::LineSearch::MoreThuente)
      <li> "User Defined" - see NOX::LineSearch::Factory.
      </ul>

     <li> "User Defined Constructor" [Teuchos::RCP<NOX::LineSearch::UserDefinedFactory>] - see NOX::LineSearch::Factory for details

     <li>"Full Step" [sublist] - used by NOX::LineSearch::FullStep,
         see class documentation for all details
     <ul>
     <li> "Full Step" [double]
     </ul>

     <li> "Backtrack" [sublist] - used by NOX::LineSearch::Backtrack,
          see class documentation for all details.
      <ul>
      <li> "Minimum Step" [double]
      <li> "Default Step" [double]
      <li> "Recovery Step" [double]
      <li> "Maximum Step" [double]
      <li> "Max Iters" [int]
      <li> "Reduction Factor" [double]
      </ul>

     <li> "Polynomial" [sublist] - used by NOX::LineSearch::MoreThuente,
          see class documentation for details.
      <ul>
      <li> "Default Step" [double]
      <li> "Max Iters" [int]
      <li> "Minimum Step" [double]
      <li> "Recovery Step" [double]
      <li> "Interpolation Type" [std::string] Valid options are:
           <ul>
           <li> "Quadratic"
           <li> "Quadratic3"
           <li> "Cubit"
           </ul>
      <li> "Min Bounds Factor" [double]
      <li> "Max Bounds Factor" [double]
      <li> "Sufficient Decrease Condition" [std::string] Valid options are:
           <ul>
           <li> "Armijo-Goldstein"
           <li> "Ared/Pred"
           </ul>
      <li> "Alpha Factor" [double]
      <li> "Force Interpolation" [boolean]
      <li> "Use Counters" [boolean]
      <li> "Maximum Iteration for Increase" [int]
      <li> "Allowed Relative Increase" [boolean]
      </ul>

     <li> "More'-Thuente" [sublist] - used by NOX::LineSearch::MoreThuente,
     see class documentation for all details.
      <ul>
      <li> "Sufficient Decrease" [double]
      <li> "Curvature Condition" [double]
      <li> "Interval Width" [double]
      <li> "Minimum Step" [double]
      <li> "Maximum Step" [double]
      <li> "Max Iters" [int]
      <li> "Default Step" [double]
      <li> "Recovery Step Type" [double]
      <li> "Sufficient Decrease Condition" [std::string] Valid options are:
           <ul>
           <li> "Armijo-Goldstein"
           <li> "Ared/Pred"
           </ul>
      <li> "Optimize Slope Calculation" [boolean]
      </ul>
     </ul>


<li> "Trust Region" [sublist] - used by NOX::Solver::TrustRegionBased,
     see class documentation for details
     <ul>
     <li> "Minimum Trust Region Radius" [double]
     <li> "Maximum Trust Region Radius" [double]
     <li> "Minimum Improvement Ratio" [double]
     <li> "Contraction Trigger Ratio" [double]
     <li> "Contraction Factor" [double]
     <li> "Expansion Trigger Ratio" [double]
     <li> "Expansion Factor" [double]
     <li> "Recovery Step" [double]
     </ul>

<li> "Solver Options" [sublist] - used by all solvers for basic options that all solvers should support (see individual solvers to make sure the option is supported).
     <ul>
     <li> "User Defined Pre/Post Operator" [Teuchos::RCP<NOX::Abstract::PrePostOperator>] - an observer object that
     allows users to insert pre and post method calls into nox's
     solvers (before and after the NOX::Solver::Generic::step() and
     before and after the NOX::Solver::Generic::solve() methods).
     Users write a class that derives from
     NOX::Abstract::PrePostOperator.  They then wrap the construction
     in a Teuchos::RCP object.  These objects can be set in
     the parameter list with the following code.

     \code
  Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
    Teuchos::rcp(new UserPrePostOperator(printing));
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator",
                     ppo);
     \endcode

     NOTE: The type of object set in the parameter list MUST be an RCP of the type NOX::Abstract::PrePostOperator.  It can NOT be a derived type!

     <li> "Status Test Check Type" [std::string] - References a NOX::StatusTest::CheckType that determines if all ("Complete") or just a minimal ("Minimal") set of status tests should be run at each iteration.  The status tests can be run efficiently if set to NOX::StatusTest::Minimal and all test run every iteration if it is set to NOX::StatusTest::Complete.  Valid options are:
     <ul>
     <li> "Complete" - All status tests are checked at each iteration.  Use this if you want to monitor progess in certain values computed in a status test.
     <li> "Minimal" - Compute only enough tests as necessary.  This is useful if some status tests are expensive.  For example, if you have 3 convergence tests in an "AND" Combo, a minimal check will only check enough tests until one is not converged.  It will then not check the others since it can't return a convecged status.
     <li> "None" - Use internally to turn off tests that should not be checked when using the "Minimal" option.  Users should not call this!
     </ul>

     <li> "User Defined Merit Function" [Teuchos::RCP<NOX::MeritFunction::Generic>] - An object that allows the user to define their own merit function for use in the globalization algorithms.  The user must implement a derived class from NOX::MeritFunction::Generic and then pass it into %NOX through the parameter list.

     \code
  Teuchos::RCP<UserMeritFunction> mf =
    Teuchos::rcp(new UserMeritFunction);
  nlParams.sublist("Solver Options").set
    <NOX::MeritFunction::Generic>("User Defined Merit Function", mf);
     \endcode

     </ul>

<li> "Printing" [sublist] - used to construct a NOX::Utils object.  See class documentation for details.  The basic options are shown below:

     <ul>

     <li> "Output Information" [sublist|int] - Either a sublist of output information booleans listed below or the integer value of the sum of the enumerated types requested.  See NOX::Utils for details on how to use the integer value.  The following sublist parameters will print information if the flag is set to true.
        <ul>
    <li> "Error" [boolean]
    <li> "Warning" [boolean]
    <li> "Outer Iteration" [boolean]
    <li> "Inner Iteration" [boolean]
    <li> "Parameters" [boolean]
    <li> "Details" [boolean]
    <li> "Outer Iteration StatusTest" [boolean]
    <li> "Linear Solver Details" [boolean]
    <li> "Test Details" [boolean]
    <li> "Stepper Iteration" [boolean]
    <li> "Stepper Details" [boolean]
    <li> "Stepper Parameters" [boolean]
    <li> "Debug" [boolean]
    </ul>
     <li> "MyPID" [int]
     <li> "Output Processor" [int]
     <li> "Output Precision" [int]

     </ul>

</ul>

A user can build a set of status tests from an independent parameter list using the NOX::StatusTest::Factory.  See this objects documentation for details on how to build that parameter list.

*/


/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page epetra_interface The Epetra Interface

 \section epetra_interface_derived Implementing the Interfaces

 The Epetra support code in %NOX implements a concrete Group and
 Vector class using Epetra objects.  The user only needs to write an
 interface that supplies the nonlinear equations in residual form, \f$
 F(x)=0 \f$.  Additionally, the user may provide details such as the
 Jacobian and Preconditioner if using an iterative linear solver.
 Below, we describe each of these interfaces.

 To link your code to %NOX using the epetra interface, you must write a
 concrete class derived from the NOX::Epetra::Interface::Required
 object.  This object defines the method for evaluating the set of
 nonlinear equations to be solved:

 <li> computeF(const Epetra_Vector& x, Epetra_Vector& RHS, const FillType flag) -
      Computes the set of nonlinear equations, \f$ F \f$, to be solved
      by %NOX. This method must be supplied by the user.
      - x - solution vector specified from %NOX.
      - RHS - the Epetra_Vector to be filled with the \f$ F \f$ values that
        correspond to the input solution vector x.
      - flag - enumerated type (see NOX::Epetra::FillType) that
        tells a users interface why computeF() was called.  %NOX
        has the ability to generate Jacobians based on numerical
        differencing using calls to computeF().  In this case, the
    user may want to compute an inexact (and hopefully cheaper)
        \f$ F \f$ since it is only used in the Jacobian (or preconditioner).

 Additonally if the user wishes to provide a Jacobian operator, they
 will need to implement a concrete class derived from the
 NOX::Epetra::Interface::Jacobian object.  This object defined the
 method for evaluating the Jacobian:

 <li> computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) - this
      is an optional method that the user can implement if they wish to
      supply their own evaluation of the Jacobian. If the user does not
      wish to supply their own Jacobian, they should implement this method
      so that it throws an error if it is called.  This method should update
      the Jac operator so that subsequent Epetra_Operator::Apply() calls
      on that operator correspond to the Jacobian at the current solution
      vector x.
      - x - solution vector specified from %NOX.
      - Jac - a reference to the Jacobian operator that the user supplied
        in the NOX::Epetra::Group constructor.

 Additonally if the user wishes to provide a preconditioner operator,
 they will need to implement a concrete class derived from the
 NOX::Epetra::Interface::Preconditioner object.  This object defined
 the method for evaluating the Preconditioner or a Matrix to be used
 with an interal preconditioner:

 <li> computePreconditioner(const Epetra_Vector& x, Epetra_Operator&
      M, Teuchos::ParameterList* precParams) - This method allows a user to supply their own preconditioner.
      The method should compute a preconditioner based upon the solution
      vector x and store it in the Epetra_Operator M.  Subsequent calls to
      the Epetra_Operator::Apply method will apply this user supplied
      preconditioner to epetra vectors. The Epetra_Operator M can also be
      an Epetra_RowMatrix that can be used by internally constructed
      preconditioner objects such as AztecOO, Ipfack, and ML.

      - x - solution vector specified from %NOX.
      - M - a reference to the operator that is to be filled.  This operator
    should be updated to correspond to the current solution vector x.
      - precParams - a pointer to a parameter list that can be used by
        the preconditioner.

 \section epetra_interface_operators Optional Operators

  The user can write their own operators for the Jacobian or
  preconditioner as discussed above, or they can use one of the %NOX
  pre-supplied operators.  Available Operators include:

  - Jacobian-Free Newton-Krylov Operator (see NOX::Epetra::MatrixFree)
  - Finite Difference Operator (see NOX::Epetra::FiniteDifference)
  - Colored Finite Difference Operator (see NOX::Epetra::FiniteDifferenceColoring)
  - Broyden Operator (see NOX::Epetra::BroydenOperator)

  These operators are built automatically by the linear system object
  depending on the constuctor used.

 \section epetra_interface_notes Notes on NOX/Epetra:

  - The concrete Epetra implementations for the %NOX Group and Vector
    are in a separate library from the nox solver algorithms.  To build
    the library for nox epetra support use the flag --enable-nox-epetra
    in NOX's configure script.  This will generate a separate library called
    libnoxepetra.a

  - In addition to the %NOX headers, be sure to include the %NOX %Epetra
    support specific headers in your interface:
    \verbatim
    // NOX headers for an epetra interface
    #include "NOX.H"
    #include "NOX_Epetra.H"
    \endverbatim

  \section epetra_features Special Epetra Features

  The %Epetra support in %NOX provides a number of tools that allow users to perform complex algorithms.  This section lists some of the important features.

  <ul>
  <li>  Multivector Support.  %NOX has an implementation of multivector support that is typically used for solving multiple RHS vectors in an efficeint manner.  See NOX::Abstract::Multivector for more information.
  <li>  Scaling.  %NOX supports a scaling object that allows users to specify diagonal scaling matrices as Epetra_Vectors that will scale the linear system using either right or left scaling.  See NOX::Epetra::Scaling for more information.
  <li> Arbitrary vector space.  The %NOX algorithms were written to support arbitrary vectors spaces.  This allows users to define scaled norms and dot products and have them used consistently throughout the code.  See NOX::Epetra::VectorSpace for more details.
  <li> Model evaluator interface.  %NOX supports a concrete implementation of the NOX::Epetra::Interface objects that wraps an EpetraExt::ModelEvaluator.  If a code implements their residual evaluation, Jacobian, or preconditioner via the EpetraExt::ModelEvaluator class, then it can be automatically used by the %NOX's %Epetra Group. See NOX::Epetra::ModelEvaluatorInterface for more details.
  <li> The NOX::Epetra::Group implementation allows users to define their own linear solver via the NOX::Epetra::LinearSystem pure virtual class.  %NOX comes with the concrete NOX::Epetra::LinearSystemAztecOO class the uses iterative Krylov subspace methods. See NOX::Epetra::LinearSystem for more details.
  </ul>

  \section nox_tutorial_link NOX Epetra Tutorial

  See the page \ref nox_epetra_tutorial for more information on how to
  use the epetra support group.

  The directory Trilinos/packages/nox/test/epetra/1Dfem contains many
  examples and is the best resource on how to use %NOX with %Epetra.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page thyra_interface The Thyra Interface

 \section thyra_interface_derived Implementing the Interfaces

 The Thyra support code in %NOX implements a concrete Group and Vector
 class using thyra objects.  The user must implement a
 Thyra::ModelEvaluator that supplies the nonlinear equations in
 residual form, \f$ F(x)=0 \f$ and the Jacobian (optional but highly
 recommended).

 To link your code to %NOX using the thyra interface, you must write a
 concrete class derived from the Thyra::ModelEvaluator
 object.

  The directory Trilinos/packages/nox/test/epetra/Thyra_2Dsim contains
  examples and is the best resource on how to use %NOX with %Thyra.
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page petsc_interface The PETSc Interface

At present, %NOX supports a concrete implementation of the
<b>NOX::Abstract::Group</b> and <b>NOX::Abstract::Vector</b> using
PETSc data structures in the namespace NOX::Petsc.  The current
implementation has been tested against PETSc 3.0.0-p5.

\section build_nox_petsc_library Building the Library

To compile the %NOX/PETSc library \c libnoxpetsc.a, use the following
configuration sequence:

\verbatim
setenv PETSC_DIR ${HOME}/petsc-3.0.0-p5
setenv PETSC_ARCH linux
configure --enable-nox-petsc
\endverbatim
This will build \c libnoxpetsc.a using a pre-built PETSc
version 2.1.6 installed in the user's home directory on a linux machine:

If the environment variables \c PETSC_DIR and \c PETSC_ARCH have already been
defined, simply use:
\verbatim
configure --enable-nox-petsc
\endverbatim

A third option is to specify \c PETSC_DIR and \c PETSC_ARCH as arguments
to the configure script, eg
\verbatim
configure --enable-nox-petsc PETSC_DIR=${HOME}/petsc-3.0.0-p5 PETSC_ARCH=linux
\endverbatim


\section build_nox_petsc_examples Building the Examples

To compile the examples for %NOX/PETSc, an additional environment variable,
\c PETSC_BOPT, must be defined and set to the value of BOPT used to build the
user's PETSc installation.  For example, for PETSc previously built
debug-enabled for C++ using BOPT=g_c++, simply set PETSC_BOPT=g_c++.
The configure script then allows the examples to be built using the following:
\verbatim
configure --enable-nox-petsc --enable-nox-petsc-examples.
\endverbatim


\section additional_nox_petsc_build_information Additional Information

Depending on the user's PETSc build, dynamically linked libraries may be used.
If an error occurs when trying to run the example executable, make sure the
environment variable \c LD_LIBRARY_PATH points to the PETSc library directory.
From the example above, setting the path would likely require:
\verbatim
setenv LD_LIBRARY_PATH ${PETSC_DIR}/lib/lib_${PETSC_BOPT}/${PETSC_ARCH}
\endverbatim


*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_epetra_tutorial NOX Epetra Tutorial

 This page describes how to write the main file for using %NOX's
 Epetra implementation of the Group and Vector.

 This example can be found in Trilinos/packages/nox/test/epetra/1Dfem/1Dfem.C

 Begin by including the %NOX library header files and any other headers
 needed by your application:

\code
// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

// Epetra Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"

// User's application specific files
#include "1DfemInterface.H"
#include "1DfemPrePostOperator.H"

#include "Teuchos_ParameterList.hpp"
\endcode


 For convenience define the following:

\code
using namespace std;
\endcode

 Begin by initializing MPI (if building a parallel version) and create
 the Epetra communicator for MPI.  In this case autoconf defines the
 flag HAVE_MPI if we are building a parallel version.

\code
int main(int argc, char *argv[])
{

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( NOX::get_global_comm() );
#else
  Epetra_SerialComm Comm;
#endif
\endcode

 Setup some initial values.

\code
  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Get the number of elements from the command line
  int NumGlobalElements = 0;
  if ((argc > 2) && (verbose))
    NumGlobalElements = atoi(argv[2]) + 1;
  else if ((argc > 1) && (!verbose))
    NumGlobalElements = atoi(argv[1]) + 1;
  else
    NumGlobalElements = 101;

  // The number of unknowns must be at least equal to the
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements
     << " cannot be < number of processors = " << NumProc << std::endl;
    cout << "Test failed!" << std::endl;
    throw std::runtime_error("NOX Error");
  }
\endcode

 Create your object that derives from NOX::Epetra::Interface::Required so that nox has access to the set of nonlinear equations to solve.

\code
  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Teuchos::RCP<Interface> interface =
    Teuchos::rcp(new Interface(NumGlobalElements, Comm));
\endcode

 Grab the initial guess vector and initialize it and the interface object.

\code
  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = interface->getSolution();
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
    Teuchos::rcp(new NOX::Epetra::Vector(soln,
                     NOX::Epetra::Vector::CreateView));

  // Set the PDE factor (for nonlinear forcing term).  This could be specified
  // via user input.
  interface->setPDEfactor(1000.0);

  // Set the initial guess
  soln->PutScalar(1.0);
\endcode

 Create a parameter list and choose the options for the %NOX solver.  Note that the linear solver parameters in teh "Linear Solver" sublist are dependent upon what NOX::Epetra::LinearSystem object you are using.  Currently %NOX only comes with one concrete linear system implementation: NOX::EpetraLinearSystemAztecOO.

\code
  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID);
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  if (verbose)
    printParams.set("Output Information",
                 NOX::Utils::OuterIteration +
                 NOX::Utils::OuterIterationStatusTest +
                 NOX::Utils::InnerIteration +
                 NOX::Utils::LinearSolverDetails +
                 NOX::Utils::Parameters +
                 NOX::Utils::Details +
                 NOX::Utils::Warning +
                             NOX::Utils::Debug +
                 NOX::Utils::TestDetails +
                 NOX::Utils::Error);
  else
    printParams.set("Output Information", NOX::Utils::Error +
                 NOX::Utils::TestDetails);

  // Sublist for line search
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 800);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Preconditioner", "New Ifpack");
  lsParams.set("Preconditioner Reuse Policy", "Reuse");
  lsParams.set("Max Age Of Prec", 5);
\endcode

 Optionally the user can define methods that will be called before and after each nonlineaer iteration and before and after each nonlinear solve.  This is done by creating an object derived from NOX::Abstract::PrePostOperator.

\code
  // Add a user defined pre/post operator object
  Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
    Teuchos::rcp(new UserPrePostOperator(printing));
  nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator",
                     ppo);
\endcode

 Set the status test check option.

\code
  // Let's force all status tests to do a full check
  nlParams.sublist("Solver Options").
    set("Status Test Check Type", NOX::StatusTest::Complete);
\endcode

 Create a Jacobian-Free Newton-Krylov method by using a NOX::Epetra::MatrixFree operator for the Jacobian.

\code
  // 2. Matrix-Free (Epetra_Operator)
  Teuchos::RCP<NOX::Epetra::MatrixFree> MF =
    Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface,
                         *noxSoln));
\endcode

 Create a Finite Difference operator to estimate the Jacobian matrix for preconditioning.

\code
  // 3. Finite Difference (Epetra_RowMatrix)
  Teuchos::RCP<NOX::Epetra::FiniteDifference> FD =
    Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface,
                           *soln));
\endcode

 Create a linear system derived from NOX::Epetra::LinearSystem.  NOX comes with a concrete implementation of the linear system class that uses AztecOO as the linear solver called NOX::Epetra::LinearSystemAztecOO.  In this case we will use the constructor that specifies the Jacobian operator and a preconditioner matrix to be used with an internally constructed AztecOO preconditioner.

\code
  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = MF;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = FD;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                              iJac, MF,
                              iPrec, FD,
                              *soln));
\endcode

 Create the group using the linear system.

\code
  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
                    iReq,
                    initialGuess,
                    linSys));
  NOX::Epetra::Group& grp = *grpPtr;
\endcode

 Create the status tests to determine how to terminate the nonlinear solve.

\code
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);
\endcode

 Create the solver and solve the problem

\code
  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  // End Nonlinear Solver **************************************

\endcode

 Get the underlying solution vector in terms of an Epetra_Vector and do whatever you want with it.

\code
  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).
    getEpetraVector();
\endcode

 Exit simulation.

\code
  return 0;
}
\endcode
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */


/*! \page portability_issues Portability Issues

8/29/2005:  There currently is a portability issue with %LOCA to the Solaris platform using the Sun Forte 7 compilers.  %LOCA builds on this platform but does not run correctly.  It appears the problem is compiler related as the virtual table for some %LOCA objects is not setup correctly.  All of the %LOCA examples and any %LOCA tests involving continuation runs fail with a segmentation fault.  It is not known whether %LOCA runs correctly when built with version 10 of the Sun compilers.  See bugzilla bug <a href="http://software.sandia.gov/bugzilla/show_bug.cgi?id=1238"> 1238 </a> for more information on this issue.

*/
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */
