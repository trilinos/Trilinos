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
\mainpage %LOCA:  Library of Continuation Algorithms

\section loca_intro Introduction
  
%LOCA is short for "Library of Continuation Algorithms". Its objective is to compute families of solutions to
\f[
F(x,p)=0
\f]
and their bifurcations, where \f$F:\Re^n\times\Re^m\rightarrow\Re^n\f$.  %LOCA is designed to be used with the <a href="../../../doc/html/index.html">%NOX</a> nonlinear solver package.

To learn more about the theory behind %LOCA algorithms, please read the %LOCA v1.0 manual (<a href="../../doc/loca_v1.0_manual/src_latex/loca_book.ps">ps</a>, <a href="../../doc/loca_v1.0_manual/src_latex/loca_book.pdf">pdf</a>, <a href="../../doc/loca_v1.0_manual/src_latex/loca_book/index.html">html</a>).

\section loca_classes Class Overview

See \ref loca_overview.

\section loca_params Parameter Specifications

See \ref loca_parameters.
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_overview LOCA Class Overview

  %LOCA is generic continuation package based upon the <a href="../../../doc/html/index.html">NOX</a> nonlinear solver package.  The foundation for %NOX consists of the NOX::Abstract::Group and NOX::Abstract::Vector which provide abstract interfaces to the linear algebra, data structures, and nonlinear equations to be solved.  %LOCA extends these concepts via additional abstract groups that provide the interface needed for continuation and bifurcation tracking such as setting parameters and computing derivatives with respect to parameters.  Also, %LOCA provides several generic groups that implement the extended sets of nonlinear equations representing various forms of continuation and bifurcations.  Finally, %LOCA provides a stepper class that drives a %NOX nonlinear solver to compute points along a continuation curve, and several support classes that compute predictors, step sizes, etc.  Each of these are discussed in more detail below.

Unlike %NOX which can provide a range of nonlinear solvers using a single abstract interface to the nonlinear equations and linear algebra, %LOCA provides several different levels of functionality, each requiring additional information from the underlying problem.  Therefore, the interface to %LOCA is split among several abstract classes each encapsulating a different level of functionality.  To interface to %LOCA, the user need only provide implementations of those abstract classes for the functionality the user is interested in.  %LOCA provides two such interfaces:

 - LOCA::LAPACK::Group
 - LOCA::Epetra::Group

both of which implement the interface required for all levels of functionality provided by %LOCA.  

\section continuation Continuation Classes

The lowest level of functionality provided by %LOCA is that for continuation.  The LOCA::Continuation::AbstractGroup provides the interface required for single parameter continuation and is derived from the NOX::Abstract::Group.  The additional methods provided primarily consist of setting/retrieving parameter values and computing the derivative \f$\partial F/\partial p\f$ for a given parameter \f$p\f$.  For application codes that cannot provide this derivative directly, the LOCA::Continuation::FiniteDifferenceGroup provides a generic implementation using finite differencing.

Various continuation algorithms have been implemented.  These include

- Natural Continuation (LOCA::Continuation::NaturalGroup)
- Arc-length Continuation (LOCA::Continuation::ArcLengthGroup)

each of which implements the extended set of nonlinear equations to be solved at each continuation step.  

\section bifurcation Bifurcation Classes

%LOCA can locate and track several codimension 1 bifurcations:

- Turning point (LOCA::Bifurcation::TPBord::ExtendedGroup)
- Pitchfork (LOCA::Bifurcation::PitchforkBord::ExtendedGroup)
- Hopf (LOCA::Bifurcation::HopfBord::ExtendedGroup).

For each bifurcation, the ExtendedGroup for that bifurcation implements a set of nonlinear equations representing that bifurcation.  Each of these groups provide a complete implementation of the LOCA::Continuation::AbstractGroup and therefore can be used in continuation to track the bifurcation in a second parameter.  The interfaces for these bifurcation groups are encapsulated in the following abstract groups:

- LOCA::Bifurcation::TPBord::AbstractGroup (turning points and pitchforks)
- LOCA::Bifurcation::HopfBord::AbstractGroup (Hopfs)

See the relevant class documentation for more details on the additional interfaces provided by these abstract classes.

\section stepping Parameter Stepping Classes

The final component of %LOCA consits of several parameter stepping classes.  These consist of the LOCA::Stepper which drives a %NOX nonlinear solver at each parameter step, and several support classes that compute predictors (classes belonging to the LOCA::Predictor namespace) and encapsulate various step size control strategies (classes in the LOCA::StepSize namespace).

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page loca_parameters LOCA Parameter Reference Page

As in %NOX, all parameters are passed to %LOCA through parameter lists defined by the NOX::Parameter::List class.  Nearly all parameters are passed to the LOCA::Stepper class during construction, although a few can be passed to the LOCA::Abstract::Group.  Below is a listing of all of the parameters that can be passed to the stepper with links to relevant classes that use those parameters, followed by parameters that can be passed to the LOCA::Abstract::Group constructor.

\section stepper_params Stepper Parameters

The parameter list passed to the Stepper has two sublists, "LOCA" and "NOX".  The "LOCA" sublist contains parameters and sublists for the Stepper and various %LOCA support classes, while the "NOX" sublist contains parameters and sublists specific to %NOX solvers.  The %NOX sublist is passed directly to an instance of the NOX::Solver::Manager class in the Stepper and therefore has the same structure as any %NOX solver parameter list hierarchy.  Please see the <a href="../../../doc/html/parameters.html">%NOX Parameter Reference</a> page for more details on the structure of this list.  

<ul> 
<li> "NOX" sublist -- See <a href="../../../doc/html/parameters.html">%NOX Parameter Reference</a> for more details
<li> "LOCA: sublist
     <ul>
     <li> "Stepper" sublist -- used by LOCA::Stepper
          <ul>
          <li> "Continuation Method"    ["Natural", "Arc Length"]
	  <li> "Continuation Parameter" -- String name of continuation parameter
	  <li> "Initial Value" -- Initial value of continuation parameter
	  <li> "Max Value" -- Maximum value of continuation parameter
	  <li> "Min Value" -- Minimum value of continuation parameter
	  <li> "Max Steps" -- Maximum number of continuation steps
	  <li> "Max Nonlinear Iterations" -- Maximum number of nonlinear iterations per continuation step
	  <li> "Goal g"
	  <li> "Max g"
	  <li> "Initial Scale Factor"
	  <li> "Min Scale Factor"
	  <li> "Min Tangent Factor"
	  <li> "Tangent Factor Exponent"
	  <li> "Compute Eigenvalues" -- Whether to compute eigenvalues after each continuation step
	  </ul>
     <li> "Predictor" sublist -- used by LOCA::Predictor::Manager
          <ul>
	  <li> "Method" ["Constant", "Secant", "Tangent", "Random"]
	  <li> "Epsilon" -- Size of perturbation for LOCA::Predictor::Random predictor
	  <li> "First Step Predictor" sublist -- used by LOCA::Predictor::Secant -- Replicates "Predictor" sublist structure
	  </ul>
     <li> "Step Size" sublist -- used by LOCA::StepSize::Manager
          <ul>
	  <li> "Method" ["Contant", "Adaptive"]
	  <li> "Initial Step Size" -- Initial parameter step size
	  <li> "Min Step Size" -- Minimum parameter step size
	  <li> "Max Step Size" -- Maximum parameter step size
	  <li> "Aggressiveness" -- Agressiveness factor in adaptive step size adjustement
	  </ul>
     <li> "Utilities" sublist -- used by LOCA::Utils
          <ul>
	  <li> "Output Information"
	  </ul>
     <li> "Anasazi" sublist -- used by LOCA::Epetra::Group to compute eigenvalues using the Anasazi package
          <ul>
	  <li> "Block Size" -- block size
	  <li> "Arnoldi Size" -- maximum length of the Arnoldi factorization
	  <li> "NEV" -- number of requested eigenvalues
	  <li> "Tol" -- tolerance for the converged eigenvalues
	  <li> "Convergence Check" -- checks convergence every so many steps
	  <li> "Restarts" -- number of restarts allowed
	  <li> "Frequency" -- how often to recalculate eigenvalues
	  <li> "Debug Level" -- Anasazi Debug level
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
