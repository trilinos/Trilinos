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

#ifndef LOCA_SOLVER_WRAPPER_H
#define LOCA_SOLVER_WRAPPER_H

#include "NOX_Solver_Generic.H"	// base class

namespace LOCA {

  namespace Solver{

    class Wrapper : public NOX::Solver::Generic {

    public:

      Wrapper(NOX::Solver::Generic& solver);

      Wrapper(const NOX::Solver::Generic& solver);

      ~Wrapper();

      bool reset(NOX::Abstract::Group& grp, 
		 NOX::StatusTest::Generic& tests, 
		 NOX::Parameter::List& params);
  
      bool reset(NOX::Abstract::Group& grp, 
		 NOX::StatusTest::Generic& tests);
      
      NOX::StatusTest::StatusType getStatus();
      
      NOX::StatusTest::StatusType iterate();
      
      NOX::StatusTest::StatusType solve();
      
      const NOX::Abstract::Group& getSolutionGroup() const;
      
      const NOX::Abstract::Group& getPreviousSolutionGroup() const;
      
      int getNumIterations() const;
      
      const NOX::Parameter::List& getParameterList() const;

    protected:

      void resetWrapper();

    protected:
      
      //! Pointer to the solver object
      NOX::Solver::Generic* solverPtr;

      //! Pointer to solution group
      const NOX::Abstract::Group* solnGrpPtr;

      //! Pointer to previous solution group
      const NOX::Abstract::Group* oldSolnGrpPtr;

    };
  }
}
#endif
