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

#ifndef LOCA_SOLVER_GENERIC_H
#define LOCA_SOLVER_GENERIC_H

#include "NOX_Solver_Generic.H"      // Base class

#include "NOX_StatusTest_Generic.H"  // for enum StatusType
#include "NOX_Parameter_List.H"      // for NOX::PArameter::List
#include "NOX_Common.H"              // for <string>

// Forward declarations
namespace LOCA {
  namespace Abstract {
    class Group;
  }
  namespace Stepper {
    class Generic;
  }
  class ParameterVector;
}

namespace LOCA {

//! %NOX Nonlinear Solvers Namespace.
namespace Solver {

  //! Type that defines how to reset the solution in the LOCA::Solver.
  enum ResetType {InitialGuess, PreviousSolution};

//! Abstract nonlinear solver method interface.
class Generic {

public:

  //! Constructor
  Generic(const NOX::Parameter::List& params);

  //! Destructor
  virtual ~Generic();

  /*! 
    \brief Reset the LOCA::Solver for a new nonlinear solve.  This
    is an alternative to constructing a new solver which may result in
    some internal data being reconstructed. 
  */
  virtual bool reset(ResetType t) = 0;

  //! Register the ParameterVector with the solver.
  /*! These parameter values will be used in the next set of solve()/step()
    calls made by the solver.
  */
  virtual bool setParams(ParameterVector& p) = 0;

  //! Check current convergence/failure status. 
  virtual NOX::StatusTest::StatusType getStatus() = 0;

  //! Do one nonlinear method iteration and return status.
  virtual NOX::StatusTest::StatusType iterate() = 0;

  //! Solve the nonlinear problem and return final status.
  virtual NOX::StatusTest::StatusType solve() = 0;

  //! Return a reference to the current solution group.
  virtual const Abstract::Group& getSolutionGroup() const = 0;

  //! Return a reference to the previous solution group.
  virtual const Abstract::Group& getPreviousSolutionGroup() const = 0;

  //! Get number of iterations.
  virtual int getNumIterations() const = 0;

  //! Return list containing various output parameters.
  virtual const NOX::Parameter::List& getParameterList() const = 0;

  //! Return non-const parameter list containing various output parameters.  
  /*! NOTE: This is not "const" since other groups such as NOX::Stepper::Generic need to get at this list and set output parameters.
   */
  virtual NOX::Parameter::List& getParameterList() = 0;

  //! Returns a string describing the solver method.
  virtual const string getLabel() const = 0;

  //! This allows the steper routines access to the parameter list so that it can set output parameters
  friend class LOCA::Stepper::Generic;

protected:
  
  //! A copy of the main parameter list passed in by user.
  /*! Mutable because we need to copy the NOX::Solver::Manager parameter
    list into this list even when returning a const version.
  */
  mutable NOX::Parameter::List params;

private:

  //! Disallow the default constructor
  Generic();

};
} // namespace Solver
} // namespace LOCA


#endif
