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
#include "NOX_Utils.H"		    // class data element

namespace NOX {
namespace Solver {

/*!
  \brief Nonlinear solver based on a line search (i.e., damping)

  Solves \f$F(x)=0\f$ using an iterative line-search-based method.

  Each iteration, the solver does the following.

  <ul> 
  <li>Compute a search direction \f$d\f$ via a NOX::Direction method

  <li>Compute a step length \f$\alpha\f$ and update \f$x\f$ as \f$x_{\rm new}
      = x_{\rm old} + \alpha d\f$ via a NOX::LineSearch method.

  </ul>
  
  The iterations progress until the status tests determine either
  failure or convergence.

  <B>Input Parameters</B>

  The following parameter list entries are valid for this solver:

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::LineSearch::Manager constructor. Defaults to an empty list.

  - "Direction" - Sublist of the direction parameters, passed to
    the NOX::Direction::Manager constructor. Defaults to an empty list.

  <B>Output Parameters</B>

  A sublist for output parameters called "Output" will be created and
  contain the following parameters:

  - "Nonlinear Iterations" - Number of nonlinear iterations
  
  - "2-Norm or Residual" - Two-norm of final residual

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class LineSearchBased : public Generic {

public:

  //! Constructor
  /*!
    See reset() for description
   */
  LineSearchBased(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, NOX::Parameter::List& p);

  //! Destructor
  virtual ~LineSearchBased();

  virtual bool reset(NOX::Abstract::Group& xgrp, NOX::StatusTest::Generic& t, NOX::Parameter::List& p);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType iterate();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const NOX::Parameter::List& getParameterList() const;

  //! Return the line search step size from the current iteration
  virtual double getStepSize() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Current solution.
  NOX::Abstract::Group* solnPtr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Group* oldSolnPtr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldSoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  NOX::Abstract::Vector* dirPtr;
  //! Current search direction.reference.
  NOX::Abstract::Vector& dir;

  //! Stopping test.
  NOX::StatusTest::Generic* testPtr;		

  //! Input parameters.
  NOX::Parameter::List* paramsPtr;	

  //! Utils
  NOX::Utils utils;
  
  //! Linesearch. 
  NOX::LineSearch::Manager lineSearch; 

  //! %Search %Direction. 
  NOX::Direction::Manager direction; 

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

};
} // namespace Solver
} // namespace NOX

#endif

