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

#ifndef NOX_SOLVER_LINESEARCHBASED_H
#define NOX_SOLVER_LINESEARCHBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Manager.H"  // class data element
#include "NOX_Parameter_List.H"	    // class data element

namespace NOX {
namespace Solver {

/*!
  \brief Basic line search based solver
  
  The following parameters are valid for this solver:

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::LineSearch::Manager constructor. Defaults to an empty list.

  - "%Direction" - Sublist of the direction parameters, passed to
    the NOX::Direction::Manager constructor. Defaults to an empty list.

  - "Linear %Solver" - Sublist of the linear solver paramaters, passed
    to NOX::Abstract::Group::computeNewton(). Furthermore, the "Tolerance"
    within this list may be modified by the
    resetForcingTerm(). Defaults to an empty list.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class LineSearchBased : public Generic {

public:

  //! %LineSearchBased-method constructor. 
  LineSearchBased(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p);

  //! Destructor
  virtual ~LineSearchBased();

  virtual bool reset(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p);
  virtual StatusTest::StatusType getStatus();
  virtual StatusTest::StatusType iterate();
  virtual StatusTest::StatusType solve();
  virtual const Abstract::Group& getSolutionGroup() const;
  virtual const Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Parameter::List& getOutputParameters() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Current solution.
  Abstract::Group* solnptr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Abstract::Group* oldsolnptr;	
  //! Previous solution reference.
  Abstract::Group& oldsoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Abstract::Vector* dirptr;
  //! Current search direction.reference.
  Abstract::Vector& dir;

  //! Stopping test.
  StatusTest::Generic* testptr;		

  //! Input parameters.
  Parameter::List iparams;	

  //! Output paramters.
  mutable Parameter::List oparams; 

  //! Linesearch. 
  NOX::LineSearch::Manager linesearch; 

  //! %Search %Direction. 
  Direction::Manager direction; 

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int niter;                    

  //! %Status of nonlinear solver.
  StatusTest::StatusType status;

};
} // namespace Solver
} // namespace NOX

#endif

