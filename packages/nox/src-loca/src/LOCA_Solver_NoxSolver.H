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

#ifndef LOCA_SOLVER_NEWTON_H
#define LOCA_SOLVER_NEWTON_H

#include "LOCA_Solver_Generic.H"	    // base class
#include "NOX_Solver_Manager.H"	            // class data element

namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace StatusTest {
    class Generic;
  }
}

namespace LOCA {
namespace Solver {

  using namespace NOX::StatusTest;

/*!
  \brief Handle for a NOX::Solver.  
  
  This class contains via composition a NOX::Solver.  No restrictions are placed on the type of solver although we expect Newton's method will be the typical choice.

  \author Roger Pawlowski (SNL 9233)
*/

class NoxSolver : public Generic {

public:

  //! %LineSearchBased-method constructor. 
  NoxSolver(Abstract::Group& xgrp, 
	 NOX::StatusTest::Generic& t, 
	 const NOX::Parameter::List& p);

  //! Destructor
  virtual ~NoxSolver();

  // Derived methods
  virtual bool reset(ResetType t);
  virtual bool setParams(ParameterVector& p);
  virtual StatusType getStatus();
  virtual StatusType iterate();
  virtual StatusType solve();
  virtual const Abstract::Group& getSolutionGroup() const;
  virtual const Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual NOX::Parameter::List& getParameterList();
  virtual const NOX::Parameter::List& getParameterList() const;
  virtual const string getLabel() const;

protected:

  
  //@{ \name Reset Information 
  
  //! Parameter list supplied by the user.  Stored in case a reset(InitialGuess) is called.
  NOX::Parameter::List initialParams;

  //! Convergence test supplied by the user. Stored in case a reset(InitialGuess) is called.
  NOX::StatusTest::Generic& statusTest;

  //! Deep copy of the group supplied by the user.  Stored in case a reset(InitialGuess) is called.
  NOX::Abstract::Group* initialGroupPtr;
  
  //@}

  //! The Group passed in by the user and used in the solver.
  NOX::Abstract::Group& solverGroup;

  //! The NOX Solver
  NOX::Solver::Manager* noxSolverPtr;

  //! Name of the solver
  string label;

};
} // namespace Solver
} // namespace LOCA

#endif

