// $Id$ 
// $Source$ 

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

#ifndef _NOX_SOLVER_MANAGER_H
#define _NOX_SOLVER_MANAGER_H

#include "NOX_Solver_Generic.H"	// base class
#include "NOX_Common.H"		// class data element (string)
namespace NOX {
namespace Solver{

/*! 

  \brief %Manager class to control the instatiation of the objects
  derived from the NOX::Solver::Generic object.

  <B>Parameters</B>
  
  The following entries may be specified in the parameter list.

   <ul>
   <li> "Nonlinear %Solver" - Name of the solver method. Valid choices are
   <ul> 
   <li> "Line Search Based" (NOX::Solver::LineSearchBased) [Default]
   <li> "Trust Region Based" (NOX::Solver::TrustRegionBased)
   </ul>
   </ul>
   
   \deprecated The "Nonlinear %Solver" choices "Newton" and "Line
   Search" are deprecated and revert to "Line Search Based". Likewise,
   the choice "Trust Region" is deprecated and reverts to "Trust
   Region Based".

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class Manager : public Generic {

public:

  //! Empty constructor - reset called later to really construct it
  Manager();

  /*!
    \brief Constructor

    See reset() for a full description.
  */
  Manager(NOX::Abstract::Group& grp, NOX::StatusTest::Generic &t, NOX::Parameter::List& p);

  //! Destructor.
  ~Manager();

  bool reset(NOX::Abstract::Group& grp, NOX::StatusTest::Generic &t, NOX::Parameter::List& p);

  NOX::StatusTest::StatusType getStatus();

  NOX::StatusTest::StatusType iterate();

  NOX::StatusTest::StatusType solve();

  const NOX::Abstract::Group& getSolutionGroup() const;
  
  const NOX::Abstract::Group& getPreviousSolutionGroup() const;

  int getNumIterations() const;

  const NOX::Parameter::List& getParameterList() const;

private:

  //! Print a warning message that oldName is deprecated and newName should be used instead
  void deprecated(const string& oldName, const string& newName) const;

  /*! 
    \brief Check that ptr is non-null. If it is null, print a
    message that prints fname and throw an error. 
  */
  void checkNullPtr(const string& fname) const;

private:

  //! Name of method being used
  string method;

  //! Pointer to the solver object
  NOX::Solver::Generic* ptr;

};
}
}
#endif
