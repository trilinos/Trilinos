//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

/*! 
\page loca_overview LOCA Overview

\section loca_intro Introduction
  
%LOCA is generic continuation and bifurcation analysis package that is designed for large-scale applications. The algorithms are designed with minimal additional interface requirements over that needed for a Newton method to reach an equilibrium solution. %LOCA is built upon the <a href="../../doc/html/index.html">NOX</a> nonlinear solver package.  The algorithms in %NOX are generic and written to the NOX::Abstract::Group and NOX::Abstract::Vector, which provide abstract interfaces to the linear algebra, data structures, and nonlinear equations to be solved. %LOCA uses the %NOX interface and extends it via additional abstract groups that provide the interface needed for continuation and bifurcation tracking, such as setting parameters and computing derivatives with respect to parameters. 

%LOCA provides several generic groups that take the %NOX group representing the equilibrium equations and implement the extended sets of nonlinear equations representing various forms of continuation and bifurcations (such as the additional equation for arc-length continuation). These extended groups also include generic algorithms for computing the Newton step for the extended system based on the Newton step for the equilibrium equations (e.g. Sherman-Morrison-Woodbury formula for arclength continuation). They are then sent to %NOX for nonlinear solution.  

Finally, %LOCA provides a stepper class that repeatedly calls the %NOX nonlinear solver to compute points along a continuation curve. The design allows for continuation of bifurcations so two-parameter bifurcation sets can be generated. The stepper class relies on several support classes that compute predictors, step sizes, etc.  Each of these are discussed in more detail n the \ref loca_class_overview.

Unlike %NOX which can provide a range of nonlinear solvers using a single abstract interface to the nonlinear equations and linear algebra, %LOCA provides several different levels of functionality, each requiring additional information from the underlying problem.  Therefore, the interface to %LOCA is split among several abstract classes each encapsulating a different level of functionality.  To interface to %LOCA, the user need only provide implementations of those abstract classes for the functionality the user is interested in.  %LOCA provides two complete interfaces:

 - LOCA::LAPACK::Group
 - LOCA::Epetra::Group

both of which implement the interface required for all levels of functionality provided by %LOCA.

\section loca_user LOCA User Information

See \ref loca_user_info.

\section loca_classes Class Overview

See \ref loca_class_overview.

\section loca_params Parameter Specifications

See \ref loca_parameters.
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_class_overview LOCA Class Overview

\section continuation Continuation Classes

The lowest level of functionality provided by %LOCA is that for continuation.  The LOCA::Continuation::AbstractGroup provides the interface required for single parameter continuation and is derived from the NOX::Abstract::Group.  The additional methods provided primarily consist of setting/retrieving parameter values and computing the derivative \f$\partial F/\partial p\f$ for a given parameter \f$p\f$.  For application codes that cannot provide this derivative directly, the LOCA::Continuation::FiniteDifferenceGroup provides a generic implementation using finite differencing.

Multiple continuation algorithms have been implemented, including:

- Natural Continuation (LOCA::Continuation::NaturalGroup)
- Arc-length Continuation (LOCA::Continuation::ArcLengthGroup)

The continuation groups implement the extended set of nonlinear equations to be solved at each continuation step.  

\section bifurcation Bifurcation Classes

%LOCA can locate and track several codimension 1 bifurcations:

- Turning point (a.k.a. fold) (LOCA::Bifurcation::TPBord::ExtendedGroup)
- Pitchfork (LOCA::Bifurcation::PitchforkBord::ExtendedGroup)
- Hopf (LOCA::Bifurcation::HopfBord::ExtendedGroup).

For each bifurcation, the ExtendedGroup for that bifurcation implements a set of nonlinear equations representing that bifurcation.  Each of these groups provide a complete implementation of the LOCA::Continuation::AbstractGroup and therefore can be used in continuation to track the bifurcation in a second parameter.  The interfaces for these bifurcation groups are encapsulated in the following abstract groups:

- LOCA::Bifurcation::TPBord::AbstractGroup (turning points and pitchforks)
- LOCA::Bifurcation::HopfBord::AbstractGroup (Hopfs)

See the relevant class documentation for more details on the additional interfaces provided by these abstract classes.

\section stepping Parameter Stepping Classes

The final component of %LOCA consists of several parameter stepping classes.  These consist of the LOCA::Stepper which drives a %NOX nonlinear solver at each parameter step, and several support classes that compute predictors (classes belonging to the LOCA::Predictor namespace) and encapsulate various step size control strategies (classes in the LOCA::StepSize namespace).

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_parameters LOCA Parameter Reference Page

As in %NOX, all parameters are passed to %LOCA through parameter lists defined by the NOX::Parameter::List class.  Nearly all parameters are passed to the LOCA::Stepper class during construction, although a few can be passed to the LOCA::Abstract::Group.  Below is a listing of all of the parameters that can be passed to the stepper with links to relevant classes that use those parameters, followed by parameters that can be passed to the LOCA::Abstract::Group constructor.

\section stepper_params Stepper Parameters

The parameter list passed to the Stepper has two sublists, "LOCA" and "NOX".  The "LOCA" sublist contains parameters and sublists for the Stepper and various %LOCA support classes, while the "NOX" sublist contains parameters and sublists specific to %NOX solvers.  The %NOX sublist is passed directly to an instance of the NOX::Solver::Manager class in the Stepper and therefore has the same structure as any %NOX solver parameter list hierarchy.  Please see the <a href="../../doc/html/parameters.html">%NOX Parameter Reference</a> page for more details on the structure of this list.  

<ul> 
<li> "NOX" sublist -- See <a href="../../doc/html/parameters.html">%NOX Parameter Reference</a> for more details
<li> "LOCA: sublist
     <ul>
     <li> "Stepper" sublist -- used by LOCA::Stepper to set parameters relevant for continuation run
          <ul>
          <li> "Continuation Method" -- [string] (default: "Natural") Type of continuation to use.  Choices are:
	          <ul>
                  <li> "Natural" -- Natural continuation (LOCA::Continuation::NaturalGroup)
                  <li> "Arc Length" -- Pseudo arc-length continuation using bordering (LOCA::Continuation::ArcLengthGroup)
                  <li> "Householder Arc Length" -- Pseudo arc-length continuation using Householder transformations (LOCA::Continuation::HouseholderGroup)
                  </ul>
	  <li> "Continuation Parameter" -- [string] (Must be supplied) Name of continuation parameter
	  <li> "Initial Value" -- [double] (Must be supplied) Initial value of continuation parameter
	  <li> "Max Value" -- [double] (Must be supplied) Maximum value of continuation parameter
	  <li> "Min Value" -- [double] (Must be supplied) Minimum value of continuation parameter
	  <li> "Max Steps" -- [int] (default: 100) Maximum number of continuation steps (including failed steps)
	  <li> "Max Nonlinear Iterations" -- [int] (default: 15) Maximum number of nonlinear iterations per continuation step
          <li> "Enable Arc Length Scaling" -- [bool] (default: true) Enable arc-length scaling to equilibrate solution and parameter components to arc-length equations (see LOCA::Continuation::ArcLengthGroup)
	  <li> "Goal Arc Length Parameter Contribution" -- [double] (default: 0.5) Goal for parameter contribution to arc-length equation
	  <li> "Max Arc Length Parameter Contribution" -- [double] (default: 0.0) Max for parameter contribution to arc-length equation, triggering rescaling
	  <li> "Initial Scale Factor" -- [double] (default: 1.0) Initial scale factor for parameter term of arc-length equation
	  <li> "Min Scale Factor" -- [double] (default: 1.0e-3) Minimum scale factor for scaling parameter term of arc-length equation
	  <li> "Enable Tangent Factor Step Size Scaling" -- [bool] (default: false) Enable step size scaling by cosine between two consective tangent vectors \f$v_0\f$ and \f$v_1\f$ to continuation curve:  \f$|v_0\cdot v_1|^{\alpha}\f$ where \f$\alpha\f$ is the tangent factor exponent.
	  <li> "Min Tangent Factor" -- [double] (default: 0.1) Minumum cosine between two consecutive tangent vectors, below which the continuation step is failed
	  <li> "Tangent Factor Exponent" -- [double] (default: 1.0) Exponent on the cosine between two consecutive tangent vectors, which then modifies the step size
	  <li> "Compute Eigenvalues" -- [bool] (default: false) Flag for requesting eigenvalue calculation after each continuation step
	  </ul>
     <li> "Bifurcation" sublist -- used by LOCA::Bifurcation::Manager to determine what type of bifurcation calculation, if any, to use
          <ul>
          <li> "Method" -- [string] (default: "None") Bifurcation method to use.  Choices are
               <ul>
               <li> "None" -- No bifurcation and do equilibrium continuation instead
               <li> "Turning Point" -- Turning point tracking using bordering (LOCA::Bifurcation::TPBord::ExtendedGroup) 
               <li> "Modified Turning Point" -- Turning point tracking using modified bordering to avoid singular Jacobian solves (LOCA::Bifurcation::TPBord::ModifiedBorderingGroup) 
               <li> "Nic-Day Modified Turning Point" -- Turning point tracking using a combination of Nic-Day deflation ideas with the modified turning point approach (LOCA::Bifurcation::TPBord::NicDayModifiedBorderingGroup)
               <li> "Pitchfork" -- Pitchfork tracking using bordering (LOCA::Bifurcation::PitchforkBord::ExtendedGroup) 
               <li> "Hopf" -- Hopf point tracking using bordering (LOCA::Bifurcation::HopfBord::ExtendedGroup)
               </ul>
          <li> "Bifurcation Parameter" -- [string] (Must be supplied if "Method" is not "None") Name of bifurcation parameter
          <li> "Length Normalization Vector" -- [NOX::Abstract::Vector*] (Must be supplied if "Method" is not "None) Pointer to NOX::Abstract::Vector holding length normilization vector \f$\phi\f$ used in turning point, pitchfork, and Hopf bifurcation algorithms
          <li> "Initial Null Vector" -- [NOX::Abstract::Vector*] (Must be supplied for tuning point and pitchfork problems) Initial guess for null vector in turning point and pitchfork algorithms
          <li> "Asymmetric Vector" -- [NOX::Abstract::Vector*] (Must be supplied for pitchfork problems) Pointer to NOX::Abstract::Vector holding asymmetric vector \f$\psi\f$ used in pitchfork algorithm
          <li> "Initial Real  Eigenvector" -- [NOX::Abstract::Vector*] (Must be supplied for Hopf problems) Initial guess for real component of generalized eigenvector in Hopf algorithm
          <li> "Initial Imaginary  Eigenvector" -- [NOX::Abstract::Vector*] (Must be supplied for Hopf problems) Initial guess for imaginary component of generalized eigenvector in Hopf algorithm
          <li> "Initial Frequency" -- [double] (Must be supplied for Hopf problems) Initial guess for Hopf frequency \f$\omega\f$ in Hopf algorithm
          <li> "Perturb Initial Solution" -- [bool] (default: false) Flag indicating whether to apply an initial perturbation to the initial guess for the solution vector before starting bifurcation algorithm
          <li> "Relative Perturbation Size" -- [double] (default: 1.0e-3) Size of relative perturbation of initial guess for solution vector
          </ul>
     <li> "Predictor" sublist -- used by LOCA::Predictor::Manager to determine what type of predictor to use for each continuation step
          <ul>
	  <li> "Method" -- [string] (default: "Constant") Predictor method to use for computing the initial guess for each continuation step.  Choices are 
                 <ul>
                 <li> "Constant" -- Use previous solution as initial guess for next step (LOCA::Predictor::Constant)
                 <li> "Secant" -- Use secant vector to continuation curve to compute initial guess (LOCA::Predictor::Secant)
                 <li> "Tangent" -- Use tangent vector to continuation curve to compute initial guess (LOCA::Predictor::Tangent)
                 <li> "Random" -- Use a random perturbation of previous solution to compute initial guess (LOCA::Predictor::Random)
                 </ul> 
	  <li> "Epsilon" -- [double] (default: 1.0e-3) Relative size of perturbation for random predictor
	  <li> "First Step Predictor" sublist -- used by the secant predictor to determine which predictor to use for the first continuation step -- Replicates "Predictor" sublist structure
	  <li> "Last Step Predictor" sublist -- Predictor to use for last step of arc-length continuation to hit target (max or min) value exactly (usually "Constant" or "Random") -- Replicates "Predictor" sublist structure
	  </ul>
     <li> "Step Size" sublist -- used by LOCA::StepSize::Manager to determine step size constrol strategies
          <ul>
	  <li> "Method" -- [string] (default: "Constant") Step size control strategy to use.  Choices are
                <ul>
                <li> "Constant" -- Use a constant step size in general, reducing the step size after a failure and increasing step size back up to original value after subsequent successes (see LOCA::StepSize::Constant)
                <li> "Adaptive" -- Use an adaptive step size control strategy that adjusts step size according to the number of Newton iterations per step (see LOCA::StepSize::Adaptive)
	  <li> "Initial Step Size" -- [double] (default: 1.0) Initial parameter step size
	  <li> "Min Step Size" -- [double] (default: 1.0e-12) Minimum parameter step size
	  <li> "Max Step Size" -- [double] (default: 1.0e+12) Maximum parameter step size
	  <li> "Failed Step Reduction Factor" -- [double] (default: 0.5) Factor by which step size is reduced after a failed step
	  <li> "Successful Step Increase Factor" - [double] (default: 1.26) Factor by which step size is increased after a successful step when the step size is smaller than the initial step size (Constant step size method only). 
	  <li> "Aggressiveness" -- [double] (default: 0.0) Aggressiveness factor in adaptive step size adjustment
	  </ul>
     <li> "Utilities" sublist -- used by LOCA::Utils to control screen output
          <ul>
	  <li> "Output Information" -- [int] (default: 0xf) Sum of message types (see LOCA::Utils)
	  <li> "Output Processor" -- [int] (default: 0) ID of printing processor
	  <li> "MyPID" -- [int] (default: 0) Processor ID
	  <li> "Output Precision" -- [int] (default: 3) Precision of output
	  </ul>
     <li> "Anasazi" sublist -- used by LOCA::Continuation::AnasaziGroup to compute eigenvalues using the Anasazi package
          <ul>
	  <li> "Block Size" -- [int] (default: 1) block size
	  <li> "Arnoldi Size" -- [int] (default: 30) maximum length of the Arnoldi factorization
	  <li> "NEV" -- [int] (default: 4) number of requested eigenvalues
	  <li> "Tol" -- [double] (default: 1.0e-7) tolerance for the converged eigenvalues
	  <li> "Convergence Check" -- [int] (default: 1) checks convergence every so many steps
	  <li> "Restarts" -- [int] (default: 1) number of restarts allowed
	  <li> "Frequency" -- [int] (default: 1) how often to recalculate eigenvalues
	  <li> "Debug Level" -- [int] (default: 1) Anasazi Debug level
	  <li> "Sorting Order" -- [string" (default: "LM") Sorting order of printed eigenvalues
	  </ul>
     </ul>
</ul>

\section group_params LOCA::Abstract::Group parameters

Parameters can also be passed to the LOCA::Abstract::Group constructor.  Currently, the only parameters used by the group are passed to the LOCA::Bifurcation::TPBord::SingularSolveGroup parent class, which provides techniques for approximating the solution to a (nearly) singular set of linear equations.

<ul> 
<li> "Singular Solve" sublist -- Used by LOCA::Bifurcation::TPBord::SingularSolveGroup to determine which method to use for computing solutions to nearly singular systems
     <ul> 
     <li> "Method" -- [string] (default: "Default") Method to use.  Choices are
           <ul>
           <li> "Default" -- Use default solution method
           <li> "Nic" -- Use "Nic" approach (LOCA::SingularJacobianSolve::Nic)
           <li> "Nic-Day" -- Use "Nic-Day" approach (LOCA::SingularJacobianSolve::NicDay)
           <li> "Iterative Refinement" -- Use iterative refinement (LOCA::SingularJacobianSolve::ItRef)
           </ul>
     </ul>
</ul>

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_user_info LOCA User Information

\section loca_user_info_overview Overview

Since %LOCA is built upon %NOX, using %LOCA is quite similar to using %NOX and
an understanding of how to use %NOX is a prerequisite for using %LOCA.  See
\ref nox_user_information for a description of how to use %NOX.  Here we will
discuss the additional steps required to use %LOCA.

\section loca_user_info_step_1 Step 1: Download, Configure, Compile, and Install LOCA

%LOCA is part of %NOX, so by downloading %NOX you have %LOCA.  To enable %LOCA 
during configuration and compilation, supply the --enable-loca configuration 
option.  %LOCA %LAPACK and %Epetra support and examples are automatically 
enabled if they are enabled in %NOX (see \ref nox_configuration_options).  
Typing "make" and "make install" in the top-level %NOX directory then
compile and install %LOCA as well as %NOX.

\section loca_user_info_step_2 Step 2: Create concrete implementations of the LOCA AbstractGroup classes

See \ref loca_class_overview for a description of the %LOCA AbstractGroup class
hierarchy used to supply interfaces for various parameter continuation and
bifurcation tracking algorithms.  These AbstractGroups extended the interface
supplied by the NOX::Abstract::Group to provided the functionality needed for
continuation and bifurcation tracking, and to use %LOCA, concrete 
implementations of these interfaces must be provided.  

\subsection loca_user_info_step_2_writing_implementation Writing your own implementation

As in %NOX, we recommend using the LOCA::LAPACK::Group as a guide for writing
your own implementation.  The simplest approach is to derive your group
from the LOCA::Abstract::Group and provide implementations of the appropriate
AbstractGroup methods for the level of functionality you desire.  At a 
minimum for parameter continuation, you must provide implementations of the
parameter setting/retrieval methods.  Note there is no additional %LOCA 
abstract vector interface to implement to use %LOCA.

\subsection loca_user_info_step_2_loca_implementations Implementations provided with LOCA

%LOCA provides two ready-made implementations:

<ul>
<li>
The LOCA::LAPACK implementation extends the NOX::LAPACK implementation to
%LOCA.  %LAPACK support and examples are automatically enabled if they are
enabled in %NOX.  The LOCA::LAPACK implementation supports all levels of
functionality provided by %LOCA and can be used in any continuation or
bifurcation calculation.

<li>
The LOCA::Epetra implementation extends the NOX::Epetra implementation to
%LOCA.  %Epetra support and examples are automatically enabled if they are
enabled in %NOX.  All levels of functionality provided by LOCA are supported
except Hopf calculations and generalized eigenvalue calculations, which require
a separate mass matrix and complex-valued computations.

<li>
Currently, there is no %LOCA::Petsc support provided, but this will be
added as soon as there is demand.
</ul>

\section loca_user_info_step_3 Step 3: Call LOCA from your code

See \ref loca_continuation_tutorial for a tutorial describing how to use 
%LOCA to perform a simple continuation. The
example problems in the examples-lapack and examples-epetra
directories can serve as templates as well.

\section loca_user_info_step_4 Step 4: Link your code to LOCA

Linking with %LOCA is very similar to linking with %NOX.  You must additionally
link against -lloca and one of the concrete implementations, either 
-llocalapack (if using the LOCA::LAPACK), -llocaepetra (if using 
LOCA::Epetra),or a custom implementation.

\section loca_user_info_problems If you have problems...

If you have problems with %LOCA, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_continuation_tutorial LOCA Continuation Tutorial

\section loca_continuation_tutorial_overview Overview

Here we provide a brief tutorial for using %LOCA to perform a simple
continuation using the Chan problem (see ChanProblemInterface).  The 
code fragements discussed below can be found in ChanContinuation.C in
the Chan subdirectory of the %LOCA %LAPACK examples directory.  

The ChanProblemInterface implements a 1-D finite-difference discretization
of the Chan problem:
\f[
      \frac{d^2 T}{d x^2} + \alpha s \left(1 + \frac{x + 0.5 x^2}{1 + 0.01 x^2}\right) = 0
  \f]
subject to the boundary conditions \f$T(0) = T(1) = \beta\f$.  The parameters
are \f$\alpha\f$, \f$\beta\f$, \f$s\f$, and \f$n\f$, the size of the 
discretization.  The scaling factor \f$s\f$ is used to test continuation
at different scales of \f$\alpha\f$.  Below we will track the solution
\f$T(x,\alpha,\beta)\f$ as \f$\alpha\f$ is varied from \f$0\f$ to \f$5\f$.

\section loca_continuation_tutorial_code ChanContinuation.C Line by Line

\code
#include "LOCA.H"
#include "LOCA_LAPACK.H"
#include "ChanProblemInterface.H"
\endcode

To use %LOCA in your code, you must always include the header LOCA.H.  Since
this is a %LAPACK problem, we also include LOCA_LAPACK.H.  Finally, we include
the header for the Chan problem, ChanProblemInterface.H.

\code
int main()
{
  int n = 100;
  double alpha = 0.0;
  double beta = 0.0;
  double scale = 1.0;
  int maxNewtonIters = 20;

  alpha = alpha / scale;

  try {
\endcode

Next we set up the basic problem parameters.  \em maxNewtonIters is the
maximum number of nonlinear solver iterations we wish to take at each
continuation step.  %NOX and %LOCA do throw exceptions (but only when
serious errors occur), so all %NOX and %LOCA calculations should be placed
in a try block.

\code
    // Create output file to save solutions
    ofstream outFile("ChanContinuation.dat");
    outFile.setf(ios::scientific, ios::floatfield);
    outFile.precision(14);

    // Save size of discretizations
    outFile << n << endl;
\endcode

Here we set up a file for saving the solutions computed at each continuation
step.  The printSolution method of the ChanProblemInterface is set up in
such a way that if a file is provided, the current continuation parameter
and solution vector are appended to the file, in addition to printing a portion
of the solution vector to the screen.  The format of this file is a series
of rows, with each row containing \f$n+1\f$ numbers, the first is the 
continuation parameter with the remaining \f$n\f$ consisting of each 
component of the solution vector.

\code
    // Set up the problem interface
    ChanProblemInterface chan(n, alpha, beta, scale, outFile);
    LOCA::ParameterVector p;
    p.addParameter("alpha",alpha);
    p.addParameter("beta",beta);
    p.addParameter("scale",scale);
\endcode

Next we instantiate the ChanProblemInterface and create a parameter vector
to store the values of the problem parameters.  Note that it is not 
necessary to put every problem parameter into the parameter vector, only those
that serve as possible continuation or bifurcation parameters need to be
supplied.

\code
    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    LOCA::LAPACK::Group grp(chan);
    
    grp.setParams(p);
\endcode

Next we instantiate the %LAPACK group with the Chan problem and then set
the parameter vector in the group.  From this point on, the LOCA::Stepper,
via the LOCA::LAPACK::Group, will take care of setting parameters in
the problem interface.

\code
   // Create parameter list
    NOX::Parameter::List paramList;

    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.setParameter("Continuation Method", "Natural");
    stepperList.setParameter("Continuation Method", "Arc Length");
    stepperList.setParameter("Continuation Parameter", "alpha");
    stepperList.setParameter("Initial Value", alpha);
    stepperList.setParameter("Max Value", 5.0/scale);
    stepperList.setParameter("Min Value", 0.0/scale);
    stepperList.setParameter("Max Steps", 50);
    stepperList.setParameter("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.setParameter("Enable Arc Length Scaling", true);
    stepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
    stepperList.setParameter("Initial Scale Factor", 1.0);
    stepperList.setParameter("Min Scale Factor", 1.0e-8);
    stepperList.setParameter("Enable Tangent Factor Step Size Scaling",true);
    stepperList.setParameter("Min Tangent Factor", -1.0);
    stepperList.setParameter("Tangent Factor Exponent",1.0);
    stepperList.setParameter("Compute Eigenvalues",false);

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "None");

    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Tangent");
    //predictorList.setParameter("Method", "Secant");

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.setParameter("Method", "Constant");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1/scale);
    stepSizeList.setParameter("Min Step Size", 1.0e-3/scale);
    stepSizeList.setParameter("Max Step Size", 10.0/scale);
    stepSizeList.setParameter("Aggressiveness", 0.5);

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("Output Information", 
			       LOCA::Utils::Warning +
			       LOCA::Utils::StepperIteration +
   			       LOCA::Utils::StepperDetails +
			       LOCA::Utils::Solver +
			       LOCA::Utils::Parameters +
			       LOCA::Utils::SolverDetails);
\endcode

Next we set up the %LOCA parameters.  We are setting up the problem to 
perform arc-length continuation in the parameter "alpha" from 0 to 5 with 
a maximum of 50 continuation steps and \em maxNewtonIters nonlinear iterations
per step.  The "Enable Arc Length Scaling", "Goal/Max Arc Length Parameter Contribution" and "Scale Factor" parameters deal with arc-length
parameter scaling to equilibrate relative contributions between solution
and parameter components in the arc-length equations.  The "Tangent Factor" parameters deal with reducing
the step size when the solution changes by a large amount from step to 
step (see \ref loca_parameters for more details).  Since we are doing an
equilibrium continuation, we set the bifurcation method to "None".  We use a tangent predictor
and adaptive step size control with an initial step size of 0.1, maximum of 
10.0 and minimum of 0.001.  The "Aggressiveness" parameter determines how
quickly the step size is increased every successful continuation step.
Finally, we specify how much information %LOCA should print out during
the continuation run.

\code
    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList.sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("Output Information", 
			  NOX::Utils::Details +
			  NOX::Utils::OuterIteration + 
			  NOX::Utils::InnerIteration + 
			  NOX::Utils::Warning);
\endcode

Next we set up the %NOX parameters.  We use a simple full-step Newton method 
for the nonlinear solve at each continuation step.

\code
    // Set up the status tests
    NOX::StatusTest::NormF normF(1.0e-8);
    NOX::StatusTest::MaxIters maxIters(maxNewtonIters);
    NOX::StatusTest::Combo comboOR(NOX::StatusTest::Combo::OR, 
				   normF, 
				   maxIters);
\endcode

Next we create appropriate status tests for the problem.  For convergence
at each continuation step, we require the extended (solution and parameter
components) residual norm be smaller than 1.0e-8 and the number of nonlinear
iterations be smaller than \em maxNewtonIters.

\code
    // Create the stepper  
    LOCA::Stepper stepper(grp, comboOR, paramList);

    // Perform continuation run
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished)
      cout << "Stepper failed to converge!" << endl;
\endcode

Finally we instantiate the stepper, run the continuation, and check the
returned status.

\code
    // Get the final solution from the stepper
    const LOCA::LAPACK::Group& finalGroup = 
      dynamic_cast<const LOCA::LAPACK::Group&>(stepper.getSolutionGroup());
    const NOX::LAPACK::Vector& finalSolution = 
      dynamic_cast<const NOX::LAPACK::Vector&>(finalGroup.getX());

    // Output the parameter list
    if (LOCA::Utils::doPrint(LOCA::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList().print(cout);
      cout << endl;
    }

    outFile.close();
  }

  catch (string& s) {
    cout << s << endl;
  }
  catch (char *s) {
    cout << s << endl;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
  }

  return 0;
}
\endcode

Lastly we copy the final solution out of the stepper, print out the final
parameter list, close the output file, and catch any thrown exceptions.

After running the example and plotting the maximum of the temperature versus the continuation parameter \f$\alpha\f$ at each step, we  obtain the following continuation curve with two turning points:

\image html chan.png

\image latex chan.eps

This plot was generated via MATLAB using the output file ChanContinuation.dat specified above.  For those interested, the MATLAB commands used to generate this plot are shown below.
\code

% open output file
fid = fopen('ChanContinuation.dat');

% read dimension of discretization
n = fscanf(fid, '%d', 1);
  
alpha = []; % array of continuation parameter values at each step
x = [];     % array of solution components at each step
  
% read values from output file
while ~feof(fid)
  
  % read alpha
  alpha = [alpha fscanf(fid, '%g', 1)];
  
  % read x
  x = [x fscanf(fid, '%g', n)];
  
end

% close output file
fclose(fid);

% compute maximum of each temperature profile
maxT = max(x);

plot(alpha,maxT,'bo-');
xlabel('\alpha');
ylabel('T_{max}   ','Rotation',0);
title('Arc-length Continuation:  \beta = 0');

\endcode
*/

