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

#ifndef NOX_ABSTRACT_PREPOSTOPERATOR_H
#define NOX_ABSTRACT_PREPOSTOPERATOR_H

#include "NOX_Common.H"  // for NOX_Config.h

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}

namespace NOX {
namespace Abstract {

/*! 

  \brief %NOX's pure virtual class to allow users to insert pre and
  post operations into nox's solvers (before and after the
  NOX::Solver::Generic::iterate() and NOX::Solver::Generic::solve()
  methods).

  The user should implement their own concrete implementation of this
  class and register it as a
  Teuchos::RCP<NOX::Abstract::PrePostoperator> in the "Solver
  Options" sublist.

  
  To create and use a user defined pre/post operators:

  <ol>

  <li> Create a pre/post operator that derives from
  NOX::Abstract::PrePostOperator. For example, the pre/post operator \c
  Foo might be defined as shown below.

  \code
  class Foo : public NOX::Abstract::PrePostOperator {
  // Insert class definition here
  }
  \endcode

  <li> Create the appropriate entries in the parameter list, as follows.

  \code
  Teuchos::RCP<Foo> foo = Teuchos::rcp(new Foo);
  params.sublist("Sovler Options").set<NOX::Abstract::PrePostOperator>("User Defined Pre/Post Operator", foo);
  \endcode

  </ol>

  \author Roger Pawlowski (SNL 9233)
*/

class PrePostOperator {

public:
  
  //! %Abstract %Vector constructor (does nothing)
  PrePostOperator() {};

  //! Copy constructor (doesnothing)
  PrePostOperator(const NOX::Abstract::PrePostOperator& /*source*/) {};

  //! %Abstract %Vector destructor (does nothing)
  virtual ~PrePostOperator() {};

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
  virtual void runPreIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
  virtual void runPostIterate(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  virtual void runPreSolve(const NOX::Solver::Generic& solver);

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  virtual void runPostSolve(const NOX::Solver::Generic& solver);

}; // class PrePostOperator
} // namespace Abstract
} // namespace NOX

inline void NOX::Abstract::PrePostOperator::
runPreIterate(const NOX::Solver::Generic& /*solver*/) {
  return;
}

inline void NOX::Abstract::PrePostOperator::
runPostIterate(const NOX::Solver::Generic& /*solver*/) {
  return;
}

inline void NOX::Abstract::PrePostOperator::
runPreSolve(const NOX::Solver::Generic& /*solver*/) {
  return;
}

inline void NOX::Abstract::PrePostOperator::
runPostSolve(const NOX::Solver::Generic& /*solver*/) {
  return;
}

#endif
