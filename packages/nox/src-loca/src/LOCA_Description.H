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
  
  %LOCA is short for "Library of Continuation Algorithms". Its 
  objective is to compute families of solutions to
  \f[
  F(x,p)=0
  \f]
  and their bifurcations, where \f$F:\Re^n\times\Re^m\rightarrow\Re^n\f$.  %LOCA is designed to be used with the <a href="../../../doc/html/index.html">NOX</a> nonlinear solver package.

  \section loca_classes Class Overview

  See \ref loca_overview.
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
