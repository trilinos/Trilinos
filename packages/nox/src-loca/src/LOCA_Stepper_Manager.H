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

#ifndef LOCA_STEPPER_MANAGER_H
#define LOCA_STEPPER_MANAGER_H

#include "LOCA_Stepper_Generic.H"    // base class

#include "NOX_Common.H"              // for <string> library

// Forward Declarations
namespace NOX {
  namespace Parameter {
    class List;
  }
}
namespace LOCA {
  namespace Abstract {
    class Vector;
    class Group;
  }
}


namespace LOCA {
  
namespace Stepper {
  
  using namespace NOX::StatusTest;

/*!
  \brief Manager (handle) for the LOCA::Stepper derived classes.
  
  The following parameters are valid for the LOCA::Stepper::Generic base class:

  - "Stepper Method" - Choice of continuation algorithm to choose:
    <ul>
    <li> "Zero Order"
    <li> "First Order"
    <li> "Arc-Length" </ul> </li>

  - "Bifurcation Parameters" - sublist containing the names and values of the parameters used in the user's application.

  - "Continuation Parameter" - String with the name of the parameter from the "Bifurcation Parameter" sublist that will be used as the continuation parameter.

  - "Initial Value" - Initial value of continuation parameter.

  - "Final Value" - Final value of continuation parameter.

  - "Initial Step Size" - Initial step size of the continuation parameter.

  - "Min Step Size" - Minimum step size of continuation parameter.

  - "Max Step Size" - Maximum step size of continuation parameter.

  - "Step Size Aggressiveness" - Controls the aggressiveness in computing the next step size.  Must be a positive number and can be greater than 1.0.  Setting this to zero yields a constant step size.

  - "Max Continuation Steps" - Maximum number of continuation steps.

  - "Max Nonlinear Iterations" - Maximum number of nonlinear iterations allowed in the nonlinear solve.  This is used in computing a new step size based on the difficulty of the nonlinear solve at the last step.

  \authors Roger Pawlowski (SNL 9233), Andrew Salinger (SNL 9233)
*/

class Manager : Generic {

public:

  //! %Stepper method constructor.
  Manager(Solver::Generic& s);

  //! Destructor
  virtual ~Manager();

  virtual bool reset(Solver::Generic& s);

  virtual bool resetFromFailedStep();

  virtual StatusType getStatus();

  virtual StatusType solve();

  virtual StatusType step();

  virtual const Abstract::Group& getSolutionGroup() const;

  virtual const Abstract::Group& getPreviousSolutionGroup() const;

  //! Returns the number of continuation steps taken by LOCA
  virtual int getNumContinuationSteps() const;

  virtual const NOX::Parameter::List& getParameterList() const;

protected:

  //! Name of the method that the stepper pointer has instantiated
  string method;

  //! Pointer to the stepper object
  Generic* ptr;

};
} // namespace Stepper
} // namespace LOCA

#endif

