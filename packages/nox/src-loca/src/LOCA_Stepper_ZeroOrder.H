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

#ifndef LOCA_SOLVER_ZERO_ORDER_H
#define LOCA_SOLVER_ZERO_ORDER_H

#include "LOCA_Stepper_Generic.H"    // base class

namespace LOCA {

  namespace Solver {
    class Generic;
  }

namespace Stepper {

  using namespace NOX::StatusTest;

/*!
  \brief %Zero Order continuation method for LOCA.
  
  The following parameters are valid for this solver:

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::Linesearch::Manager constructor. Defaults to an empty list.

  - "%Direction" - Sublist of the direction parameters, passed to
    the NOX::Direction::Manager constructor. Defaults to an empty list.

  - "Linear %Solver" - Sublist of the linear solver paramaters, passed
    to Abstract::Group::computeNewton(). Furthermore, the "Tolerance"
    within this list may be modified by the
    resetForcingTerm(). Defaults to an empty list.

  \author Roger Pawlowski (SNL 9233)
*/

class ZeroOrder : public Generic {

public:

  //! %ZeroOrder method constructor. 
  ZeroOrder(Solver::Generic& s);

  //! Destructor
  virtual ~ZeroOrder();

  virtual bool reset(Solver::Generic& s);

  virtual bool resetFromFailedStep();

  virtual StatusType getStatus();

  virtual StatusType solve();

  virtual StatusType step();

  virtual const Abstract::Group& getSolutionGroup() const;

  virtual const Abstract::Group& getPreviousSolutionGroup() const;

  virtual int getNumContinuationSteps() const;

  virtual const NOX::Parameter::List& getParameterList() const;

protected:

};
} // namespace Solver
} // namespace NOX

#endif
