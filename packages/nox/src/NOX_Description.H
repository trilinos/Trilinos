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
  \mainpage %NOX: An Object-Oriented Nonlinear Solver Package

  \section intro Introduction
  
  %NOX is short for "Nonlinear Object-Oriented Solutions". Its 
  objective is to enable efficient solution of the equation:
  \f[
  F(x)=0.
  \f]

  \section started Getting Started

  To configure and compile %NOX, follow the 
  <A HREF="../configure.html">Configuration and compilation instructions</A>. 

  \section bugs Reporting Bugs and Making Enhancement Requests

  To reports bugs or make enhancement reuqests, visit 
  <A HREF="http://csmr.ca.sandia.gov/bugzilla/">NOX's Bugzilla (Bug Tracking) Database</A>,
  and use the following instructions.
      <UL>
      <LI>Click on "Enter a new bug report"
      <LI>Choose "%NOX"
      <LI>Either login or create a new account
      <LI>Submit your bug report
      </UL>
  </UL>

  \section email Mailing Lists

  <p>We recommend that all %NOX users subscribe to the 
  <a href="http://csmr.ca.sandia.gov/mailman/listinfo/nox-users">NOX-Users Mailing List</a>.
  This mailing list is used to keep users informed of new releases and patches. 
  Developers should subscribe to both the 
  <a href="http://csmr.ca.sandia.gov/mailman/listinfo/nox-developers">NOX-Developers Mailing List</a>
  and the 
  <a href="http://csmr.ca.sandia.gov/mailman/listinfo/nox-cvs">NOX-CVS Mailing List</a>.
  The first is a mailing list dedicated to discussion of the
  development of %NOX and is archived; the second is the mailing list
  to which all CVS commit messages are sent.

  \section classes Class Overview

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
  the class decription.

  Included with %NOX are three such interfaces:

  - NOX::Example::Vector and NOX::Example::Group
  - NOX::Epetra::Vector and NOX::Epetra::Group
  - NOX::Petsc::Vector and NOX::Petsc::Group

  The \b Example interface is a very simple interface that is not
  intended for actual computations, but to serve as an
  easy-to-understand example of how one might interface to %NOX.
  
  The \b Epetra interface is an interface to the
  Trilinos/Epetra library.

  The \b PETSc interface is an interface with the PETSc library. PETSc
  was developed at Argonne National Labs; for more information see the
  <a href="http://www-fp.mcs.anl.gov/petsc">PETSc Home Page</a>

  All solvers are in the NOX::Solver namespace. The solvers are
  accessed via the NOX::Solver::Manager. The recommended solver is
  NOX::Solver::LineSearchBased, which is a basic nonlinear solver
  based on a line search. Each solver has a number of options that can
  be specified, as documented in each class.

  The search directions are in the NOX::Direction namespace and
  accessed via the NOX::Direction::Manager. The default search
  direction for a line-search based method in NOX::Direction::Newton.

  Several line searches are available, as defined in the
  NOX::LineSearch, and accessed via the NOX::LineSearch::Manager
  class. Examples include

  - NOX::LineSearch::FullStep
  - NOX::LineSearch::Backtrack
  - NOX::LineSearch::MoreThuente

  Convergence or failure of a given solver method is determined by 
  the status tests defined in the NOX::StatusTest namespace. Various
  status tests may be combined via the NOX::StatusTest::Combo
  object. Users are free to create additional status tests that derive
  from the NOX::StatusTest::Generic class.

  All parameters are specified via parameter lists as defined by the
  NOX::Parameter::List class.

  The NOX::Utils class describes some parameters that can be set to
  control output.

*/
 

