// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_SOLVER_NEWTON_H
#define _NOX_SOLVER_NEWTON_H

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
}
}

#endif

