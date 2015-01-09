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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
