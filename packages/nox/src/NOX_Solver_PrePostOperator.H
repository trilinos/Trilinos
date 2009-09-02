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

#ifndef NOX_SOLVER_PREPOSTOPERATOR_H
#define NOX_SOLVER_PREPOSTOPERATOR_H

#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"
#include "NOX_Abstract_PrePostOperator.H"

// Forward declaration
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  class Utils;
  namespace Solver {
    class Generic;
  }
  namespace Parameter {
    class PrePostOperator;
  }
}

namespace NOX {

namespace Solver {

/*!
  @brief Functor to process the pre/post operator object in the parameter list.

  This is a wrapper class for a user derived
  NOX::Abstract::PrePostOperator (ppo) object.  All solvers use this
  class so we don't have to repeat all parsing code in each
  NOX::Solver.  This class searches the "Solver Options" parameter
  list passed into the ctor and if a ppo is found will wrap the
  object.

  For instructions on how to implement a PrePostOperator, see
  NOX::Abstract::PrePostOperator.

*/
class PrePostOperator {

private:

  //! Disallow default ctor.
  PrePostOperator();

  //! Disallow copy ctor.
  PrePostOperator(const PrePostOperator& p);

  //! Disallow assignment operator.
  PrePostOperator& operator=(const PrePostOperator& p);

public:

  //! Ctor.
  PrePostOperator(const Teuchos::RCP<NOX::Utils>& utils, 
		  Teuchos::ParameterList& solverOptionsSubList);

  //! Destructor.
  virtual ~PrePostOperator();

  //! Resets the pre/post operator.
  virtual void
  reset(const Teuchos::RCP<NOX::Utils>& utils, 
	Teuchos::ParameterList& solverOptionsSublist);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
  virtual void runPreIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
  virtual void runPostIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  virtual void runPreSolve(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  virtual void runPostSolve(const NOX::Solver::Generic& solver);

protected:

  //! Flag that determines if a pre/post operator has been supplied by user.
  bool havePrePostOperator;

  //! Points to user defined pre/post operator.
  Teuchos::RCP<NOX::Abstract::PrePostOperator> prePostOperatorPtr;

};
} // namespace Solver
} // namespace NOX

inline void NOX::Solver::PrePostOperator::
runPreIterate(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreIterate(solver);
}

inline void NOX::Solver::PrePostOperator::
runPostIterate(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPostIterate(solver);
}

inline void NOX::Solver::PrePostOperator::
runPreSolve(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPreSolve(solver);
}

inline void NOX::Solver::PrePostOperator::
runPostSolve(const NOX::Solver::Generic& solver)
{
  if (havePrePostOperator)
    prePostOperatorPtr->runPostSolve(solver);
}

#endif
