// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_SOLVER_MANAGER_H
#define _NOX_SOLVER_MANAGER_H

#include "NOX_Solver_Generic.H"		// base class

#include "NOX_Status_Test.H"

namespace NOX {
namespace Solver{

//! Concrete class to control the instatiation of the derived 
//! NLS_Method class.  
/*! This is the class that takes control from the user code and hands
 *  it over to the nonlinear solver component. */

class Manager : public Generic {

public:
  //! Generic interface.
  Manager(Abstract::Group& initialGuess, Status::Test &t, Parameter::List& p);

  //! Destructor.
  ~Manager();

  void resetInputParameters(Parameter::List& p);

  Status::StatusType getStatus();

  Status::StatusType iterate();

  Status::StatusType solve();

  const Abstract::Group& getSolutionGroup() const;
  
  const Abstract::Group& getPreviousSolutionGroup() const;

  int getNumIterations() const;

  const Parameter::List& getOutputParameters() const;

private:

  //! Pointer to the solver object
  Generic* ptr;

};
}
}
#endif
