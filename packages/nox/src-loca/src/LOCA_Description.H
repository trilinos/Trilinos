//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
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

To learn more about the theory behind %LOCA algorithms, please read the %LOCA v1.0 manual (<a href="../../src-loca/doc/loca_v1.0_manual/src_latex/loca_book.ps">ps</a>, <a href="../../src-loca/doc/loca_v1.0_manual/src_latex/loca_book.pdf">pdf</a>, <a href="../../src-loca/doc/loca_v1.0_manual/src_latex/loca_book/index.html">html</a>).  Please note the theory and implementation described in this manual is for the original version of %LOCA available <a href="http://www.cs.sandia.gov/loca">here</a>.  The implementation here is completely different, but the theory is more-or-less the same.

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
<li> "%NOX" sublist -- See <a href="../../doc/html/parameters.html">%NOX Parameter Reference</a> for more details
<li> "LOCA: sublist
     <ul>
     <li> "Stepper" sublist -- used by LOCA::Stepper
          <ul>
          <li> "Continuation Method" -- ["Natural", "Arc Length"] Continuation type
	  <li> "Continuation Parameter" -- [string] Name of continuation parameter
	  <li> "Initial Value" -- [double] Initial value of continuation parameter
	  <li> "Max Value" -- [double] Maximum value of continuation parameter
	  <li> "Min Value" -- [double] Minimum value of continuation parameter
	  <li> "Max Steps" -- [int] Maximum number of continuation steps
	  <li> "Max Nonlinear Iterations" -- [int] Maximum number of nonlinear iterations per continuation step
	  <li> "Goal g" -- [double] Goal for parameter contribution to unit arclength
	  <li> "Max g" -- [double] Max for parameter contribution to unit arclength, triggering rescaling
	  <li> "Initial Scale Factor" -- [double] Initial factor for scaling solution term of arclength eq.
	  <li> "Min Scale Factor" -- [double] Minimum factor for scaling solution term of arclength eq.
	  <li> "Min Tangent Factor" -- [double] Minumum cosine between two consecutive tangent vectors, below which the continuation step is failed
	  <li> "Tangent Factor Exponent" -- [double] Exponent on the cosine between two consecutive tangent vectors, which then modifies the step size
	  <li> "Compute Eigenvalues" -- [bool] Flag for requesting eigenvalue calculation after each continuation step
	  </ul>
     <li> "Predictor" sublist -- used by LOCA::Predictor::Manager
          <ul>
	  <li> "Method" -- ["Constant", "Secant", "Tangent", "Random"] Predictor method
	  <li> "Epsilon" -- [double] Size of perturbation for LOCA::Predictor::Random predictor
	  <li> "First Step Predictor" sublist -- used by LOCA::Predictor::Secant -- Replicates "Predictor" sublist structure
	  </ul>
     <li> "Step Size" sublist -- used by LOCA::StepSize::Manager
          <ul>
	  <li> "Method" ["Constant", "Adaptive"]
	  <li> "Initial Step Size" -- [double] Initial parameter step size
	  <li> "Min Step Size" -- [double] Minimum parameter step size
	  <li> "Max Step Size" -- [double] Maximum parameter step size
	  <li> "Aggressiveness" -- [double] Aggressiveness factor in adaptive step size adjustment
	  </ul>
     <li> "Utilities" sublist -- used by LOCA::Utils
          <ul>
	  <li> "Output Information" -- [int] Sum of message types (see LOCA::Utils)
	  <li> "Output Processor" -- [int] ID of printing processor
	  <li> "MyPID" -- [int] Processor ID
	  <li> "Output Precision" -- [int] Precision of output
	  </ul>
     <li> "Anasazi" sublist -- used by LOCA::Epetra::Group to compute eigenvalues using the Anasazi package
          <ul>
	  <li> "Block Size" -- [int] block size
	  <li> "Arnoldi Size" -- [int] maximum length of the Arnoldi factorization
	  <li> "NEV" -- [int] number of requested eigenvalues
	  <li> "Tol" -- [double] tolerance for the converged eigenvalues
	  <li> "Convergence Check" -- [int] checks convergence every so many steps
	  <li> "Restarts" -- [int] number of restarts allowed
	  <li> "Frequency" -- [int] how often to recalculate eigenvalues
	  <li> "Debug Level" -- [int] Anasazi Debug level
	  </ul>
     </ul>
</ul>

\section group_params LOCA::Abstract::Group parameters

Parameters can also be passed to the LOCA::Abstract::Group constructor.  Currently, the only parameters used by the group are passed to the LOCA::Bifurcation::TPBord::SingularSolveGroup parent class, which provides techniques for approximating the solution to a (nearly) singular set of linear equations.

<ul> 
<li> "Singular Solve" sublist -- Used by LOCA::Bifurcation::TPBord::SingularSolveGroup
     <ul> 
     <li> "Method" ["Default", "Nic", "Nic-Day", "Iterative Refinement"]
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
option.  %LOCA %LAPACK and Epetra support and examples are automatically 
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
The LOCA::Epetra implementation extenders the NOX::Epetra implementation to
%LOCA.  Epetra support and examples are automatically enabled if they are
enabled in %NOX.  All levels of functionality provided by LOCA are supported
except Hopf calculations and generalized eigenvalue calculations, which require
a separate mass matrix and complex-valued computations.

<li>
Currently, there is no %LOCA::Petsc support provided, but this will be
added as soon as there is demand.
</ul>

\section loca_user_info_step_3 Step 3: Call LOCA from your code

A tutorial describing how to use %LOCA is coming shortly. The
example problems in the examples-lapack and examples-epetra
directories can serve as templates.

\section loca_user_info_step_4 Step 4: Link your code to LOCA

Linking with %LOCA is very similar to linking with %NOX.  You must additionally
link against -lloca and one of the concrete implementations, either -lloca-lapack 
(if using the LOCA::LAPACK), -lloca-epetra (if using LOCA::Epetra), or a
custom implementation.

\section loca_user_info_problems If you have problems...

If you have problems with %LOCA, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

*/
