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

#ifndef NOX_SOLVER_NEWTON_H
#define NOX_SOLVER_NEWTON_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_Linesearch_Manager.H" // line search
#include "NOX_Status_Test.H"	    // status test for convergence

namespace NOX {
namespace Solver {

//! Newton-like solver with a line search.
class Newton : public Generic {

public:

  //! Newton-method constructor. 
  Newton(Abstract::Group& xgrp, Status::Test& t, Parameter::List& p);

  //! Destructor
  virtual ~Newton();

  virtual void resetInputParameters(Parameter::List& p);
  virtual Status::StatusType getStatus();
  virtual Status::StatusType iterate();
  virtual Status::StatusType solve();
  virtual Abstract::Group& getSolutionGroup() const;
  virtual Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Parameter::List& getOutputParameters() const;

protected:
  
  virtual void resetForcingTerm();
  virtual void printUpdate();

protected:
  
  //! Current solution.
  Abstract::Group& soln;		

  //! Previous solution.
  Abstract::Group* oldsolnptr;	
  Abstract::Group& oldsoln;	

  //! Current search direction.
  Abstract::Vector* dirptr;
  Abstract::Vector& dir;

  //! Stopping test.
  Status::Test& test;		

  //! Input parameters.
  Parameter::List& iparams;	

  //! Output paramters.
  mutable Parameter::List oparams; 

  //! Linesearch. 
  Linesearch::Manager linesearch; 

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int niter;                    

  //! Status of nonlinear solver.
  Status::StatusType status;

};
} // namespace Solver
} // namespace NOX

#endif

