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

#ifndef NOX_SOLVER_GENERIC_H
#define NOX_SOLVER_GENERIC_H

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
} // namespace Solver
} // namespace NOX


#endif
