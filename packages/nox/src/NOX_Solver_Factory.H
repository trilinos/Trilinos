// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_SOLVER_FACTORY_H
#define NOX_SOLVER_FACTORY_H

#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace NOX {

  namespace Abstract {
    class Group;
  }
  namespace StatusTest {
    class Generic;
  }

  namespace Solver{

    class Generic;
    
    /*! \brief %Factory class to control the creation of solvers derived from the NOX::Solver::Generic object.
    
    <B>Parameters</B>
    
    The following entries may be specified in the parameter list.
    
    <ul>
    <li> "Nonlinear Solver" <std::string> Name of the solver method. Valid choices are
    <ul> 
    <li> "Line Search Based" or "Newton" (NOX::Solver::LineSearchBased) [Default]
    <li> "Trust Region Based" (NOX::Solver::TrustRegionBased)
    <li> "Inexact Trust Region Based" (NOX::Solver::InexactTrustRegionBased)
    <li> "Tensor Based" (NOX::Solver::TensorBased)
    </ul>
    <li> "Line Search Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::LineSearchBased solver.
    <li> "Trust Region Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TrustRegionBased solver.
    <li> "Inexact Trust Region Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::InexactTrustRegionBased solver.
    <li> "Tensor Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TensorBased solver.
    <li> "Tensor-Krylov Based" <Teuchos::ParameterList> Sublist for the NOX::Solver::TensorBasedTest solver. (Prerelease only)
    </ul>
    
    \note Solvers can be constructed using a nonmember function NOX::Solver::buildSolver instead of using this object directly.

    \author Roger Pawlowski (SNL 1416)
    */
    
    class Factory {
      
    public:
      
      //! Constructor.
      Factory();
      
      //! Destructor.
      ~Factory();
      
      Teuchos::RCP<NOX::Solver::Generic> 
      buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp, 
		  const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
		  const Teuchos::RCP<Teuchos::ParameterList>& params);
      
    };
    
    /*! \brief Nonmember helper function for the NOX::Solver::Factory.
      
    \relates NOX::Solver::Factory
    
    */

    Teuchos::RCP<NOX::Solver::Generic> 
    buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp, 
		const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
		const Teuchos::RCP<Teuchos::ParameterList>& params);
    
  }  // namespace Solver
}  // namespace NOX
#endif
