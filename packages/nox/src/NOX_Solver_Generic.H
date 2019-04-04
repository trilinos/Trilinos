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

#ifndef NOX_SOLVER_GENERIC_H
#define NOX_SOLVER_GENERIC_H

#include "NOX_StatusTest_Generic.H"  // for NOX::StatusTest::StatusType
#include "Teuchos_RCP.hpp"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  struct SolverStats;
  namespace Abstract {
    class Vector;
    class Group;
  }
}

namespace NOX {

/*!
  \brief %NOX nonlinear solvers namespace.

  The Generic class spells out the interface. The Manager class can
  instantiate any other concrete implementation. The remaining classes
  are concrete implementations.
*/
namespace Solver {

/*!
  \brief %Abstract nonlinear solver method interface.

  Defines the type of access methods into the iterative nonlinear solvers.

  <ul>
  <li>Instantiate or reset() the solver.
  <li>Find the solution via solve() or perform a single iterations via iterate().
  <li>Get information about the current solver state via
      getSolutionGroup(), getPreviousSolutionGroup(),
      getNumIterations(), and getList() --- particularily
      useful for NOX::StatusTest methods.
  <li>Get the current status of the solver via getStatus().
  </ul>
 */
class Generic {

public:

  //! Constructor (does nothing)
  Generic() {};

  //! Destructor (does nothing)
  virtual ~Generic() {};

  /** \brief Resets the solver for another solve. This resets the
      counters and status only. Uses the final solution from the last
      solve as the initial guess for the next solve.

      NOTE: All NOX solvers will call reset() automatically at teh
      beginning of the solve() method. We add the reset() method to
      the solver interface for the application to call in case the
      application needs to reset counters and status manually before
      the next call to solve() is made.
  */
  virtual void reset() = 0;

  //! Resets the solver and sets a new initial guess.
  virtual void reset(const NOX::Abstract::Vector& initial_guess) = 0;

  //! Resets the solver, sets a new status test, and sets a new initial guess.
  virtual void reset(const NOX::Abstract::Vector& initial_guess,
             const Teuchos::RCP<NOX::StatusTest::Generic>& test) = 0;

  //! Do one nonlinear step in the iteration sequence and return status.
  virtual NOX::StatusTest::StatusType step() = 0;

  /*!
    \brief Solve the nonlinear problem and return final status.

    By "solve", we call iterate() until the NOX::StatusTest value is
    either NOX::StatusTest::Converged or NOX::StatusTest::Failed.
  */
  virtual NOX::StatusTest::StatusType solve() = 0;

  //! Return a reference to the current solution group.
  virtual const NOX::Abstract::Group& getSolutionGroup() const = 0;

  //! Return a reference to the previous solution group.
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const = 0;

  //! Returns the current status of the solver.
  virtual NOX::StatusTest::StatusType getStatus() const = 0;

  //! Get number of iterations.
  virtual int getNumIterations() const = 0;

  //! Return a reference to the solver parameters.
  virtual const Teuchos::ParameterList& getList() const = 0;

  //@{ \name PyTrilinos Interface functions (require RCPs)

  //! Return a RCP to the solution group
  virtual Teuchos::RCP< const NOX::Abstract::Group > getSolutionGroupPtr() const = 0;

  //! Return a RCP to the previous solution group
  virtual Teuchos::RCP< const NOX::Abstract::Group > getPreviousSolutionGroupPtr() const = 0;

  //! Return a RCP to the solver parameters.
  virtual Teuchos::RCP< const Teuchos::ParameterList > getListPtr() const = 0;

  //@}

  //! Return a RCP to the solver statistics.
  virtual Teuchos::RCP<const NOX::SolverStats> getSolverStatistics() const = 0;

};
} // namespace Solver
} // namespace NOX


#endif
