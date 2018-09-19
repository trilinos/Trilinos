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
//@HEADER

#ifndef NOX_OBSERVER_HPP
#define NOX_OBSERVER_HPP

#include "NOX_Common.H"

// Forward Declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
  namespace Abstract {
    class Vector;
  }
}

namespace NOX {

/** \brief %NOX's pure virtual class to allow users to insert user
  defined operations into nox's solvers (before and after the
  NOX::Solver::Generic::step() and NOX::Solver::Generic::solve()
  methods). This is an Observer from GoF design pattern book.

  The user should implement their own concrete implementation of this
  class and register it as a
  Teuchos::RCP<NOX::Abstract::PrePostoperator> in the "Solver
  Options" sublist.

  To create and register a user defined pre/post operator:

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
  Teuchos::RCP<NOX::Abstract::PrePostOperator> foo = Teuchos::rcp(new Foo);
  params.sublist("Solver Options").set("User Defined Pre/Post Operator", foo);
  \endcode

  </ol>
*/

class Observer {

public:

  //! Constructor
  Observer() {}

  //! Destructor
  virtual ~Observer() {}

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::step().
  virtual void runPreIterate(const NOX::Solver::Generic& solver) {}

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::step().
  virtual void runPostIterate(const NOX::Solver::Generic& solver) {}

  //! User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
  virtual void runPreSolve(const NOX::Solver::Generic& solver) {}

  //! User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
  virtual void runPostSolve(const NOX::Solver::Generic& solver) {}

  /** \brief User defined method that will be executed prior to the
      update of the solution vector during a call to
      NOX::Solver::Generic::step(). This is intended to allow users to
      adjust the direction before the solution update, typically based
      on knowledge of the problem formulation. The direction is const
      as we can't guarantee that changes to the direction won't
      violate assumptions of the solution algorithm. Users can change
      the update/direciton after a const cast, but NOX may not
      function as expected. Use at your own risk!

      \param [in] update - the direction vector that will be used to update the solution.
      \param [in] solver - the nox solver
   */
  virtual void runPreSolutionUpdate(const NOX::Abstract::Vector& update, const NOX::Solver::Generic& solver) {}

  /** \brief User defined method that will be executed after the
      update of the solution vector during a call to
      NOX::Solver::Generic::step(). This is intended to allow users to
      adjust the direction after the solution update, typically based
      on knowledge of the problem formulation (e.g. clipping negative
      mass fractions). The direction is const as we can't guarantee
      that changes to the direction won't violate assumptions of the
      solution algorithm. Users can change the update/direciton after
      a const cast, but NOX may not function as expected. Use at your
      own risk!

      \param [in] solver - the nox solver
   */
  virtual void runPostSolutionUpdate(const NOX::Solver::Generic& solver) {}

  //! User defined method that will be executed before a call to NOX::LineSearch::Generic::compute(). Only to be used in NOX::Solver::LineSearchBased!
  virtual void runPreLineSearch(const NOX::Solver::Generic& solver) {}

  //! User defined method that will be executed after a call to NOX::LineSearch::Generic::compute(). Only to be used in NOX::Solver::LineSearchBased!
  virtual void runPostLineSearch(const NOX::Solver::Generic& solver) {}
}; // class PrePostOperator

} // namespace NOX

#endif
