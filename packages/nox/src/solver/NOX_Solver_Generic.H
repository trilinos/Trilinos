// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_SOLVER_GENERIC_H
#define _NOX_SOLVER_GENERIC_H

#include "NOX_Parameter_List.H"	
#include "NOX_Abstract_Group.H"	
#include "NOX_Status_Test.H" 

namespace NOX {

//! NOX Nonlinear Solvers Namespace.
namespace Solver {

//! Abstract nonlinear solver method interface.
class Generic {

public:

  //! Constructor.
  Generic() {};

  //! Destructor.
  virtual ~Generic() {};

  //! Reset parameters for the nonlinear method.
  virtual void resetInputParameters(Parameter::List& p) = 0;

  //! Check current convergence/failure status. 
  virtual Status::StatusType getStatus() = 0;

  //! Do one nonlinear method iteration and return status.
  virtual Status::StatusType iterate() = 0;

  //! Solve the nonlinear problem and return final status.
  virtual Status::StatusType solve() = 0;

  //! Return a reference to the current solution group.
  virtual Abstract::Group& getSolutionGroup() const = 0;

  //! Return a reference to the previous solution group.
  virtual Abstract::Group& getPreviousSolutionGroup() const = 0;

  //! Get number of iterations.
  virtual int getNumIterations() const = 0;

  //! Return list containing various output parameters.
  virtual const Parameter::List& getOutputParameters() const = 0;

};
}
}


#endif
