/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file   Ifpack2_Details_LinearSolver_def.hpp
/// \author Mark Hoemmen
/// \brief  Definition of Ifpack2::Details::LinearSolver.

#ifndef IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP
#define IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP

#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Tpetra_MultiVector.hpp"

// Ifpack2: key is for Ifpack2's factory to have subordinate
// factories.  That way, each package still has one factory, but we
// don't have to worry about intrapackage circular dependencies (e.g.,
// relating to AdditiveSchwarz).  There are two approaches:
//
//   1. Reuse existing Ifpack2::Details::OneLevelFactory
//   2. Have each Ifpack2 solver register itself with Ifpack2's factory

namespace Ifpack2 {
namespace Details {

template<class SC, class LO, class GO, class NT>
LinearSolver<SC, LO, GO, NT>::
LinearSolver (const Teuchos::RCP<prec_type>& solver, const std::string& solverName) :
  solver_ (solver),
  solverName_ (solverName)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  const char prefix[] = "Ifpack2::Details::LinearSolver: ";
  TEUCHOS_TEST_FOR_EXCEPTION(solver.is_null (), std::invalid_argument,
                             prefix << "Input solver is NULL.");

  typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef ::Ifpack2::Details::CanChangeMatrix<row_matrix_type> mixin_type;
  RCP<mixin_type> innerSolver = rcp_dynamic_cast<mixin_type> (solver);
  TEUCHOS_TEST_FOR_EXCEPTION
    (innerSolver.is_null (), std::invalid_argument, prefix << "The input "
     "solver does not implement the setMatrix() feature.  Only Ifpack2 solvers "
     "that inherit from Ifpack2::Details::CanChangeMatrix implement this feature.");
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
setMatrix (const Teuchos::RCP<const OP>& A)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Tpetra::RowMatrix<SC, LO, GO, NT> row_matrix_type;
  typedef ::Ifpack2::Details::CanChangeMatrix<row_matrix_type> mixin_type;
  const char prefix[] = "Ifpack2::Details::LinearSolver::setMatrix: ";

  // It's OK for the input matrix to be null.  Ifpack2 solvers may
  // interpret this as a hint to clear out their state.  It's also a
  // way to keep the preconditioner around, but disassociate it from a
  // particular matrix.  (The code that uses the preconditioner might
  // not want to or be able to recreate it.)
  RCP<const row_matrix_type> A_row;
  if (! A.is_null ()) {
    A_row = rcp_dynamic_cast<const row_matrix_type> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_row.is_null (), std::invalid_argument, prefix << "The input matrix A, "
       "if not null, must be a Tpetra::RowMatrix.");
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "Solver is NULL.  "
     "This should never happen!  Please report this bug to the Ifpack2 "
     "developers.");

  RCP<mixin_type> innerSolver = rcp_dynamic_cast<mixin_type> (solver_);
  TEUCHOS_TEST_FOR_EXCEPTION
    (innerSolver.is_null (), std::logic_error, prefix << "The solver does not "
     "implement the setMatrix() feature.  Only input preconditioners that "
     "inherit from Ifpack2::Details::CanChangeMatrix implement this.  We should"
     " never get here!  Please report this bug to the Ifpack2 developers.");
  innerSolver->setMatrix (A_row);

  A_ = A; // keep a pointer to A, so that getMatrix() works
}

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<const typename LinearSolver<SC, LO, GO, NT>::OP>
LinearSolver<SC, LO, GO, NT>::
getMatrix () const {
  return A_; // may be null
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
solve (MV& X, const MV& B)
{
  const char prefix[] = "Ifpack2::Details::LinearSolver::solve: ";
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "The solver is NULL!  "
     "This should never happen.  Please report this bug to the Ifpack2 "
     "developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
     "set yet.  You must call setMatrix() with a nonnull matrix before you "
     "may call this method.");
  solver_->apply (B, X);
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  solver_->setParameters (*params);
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
symbolic ()
{
  const char prefix[] = "Ifpack2::Details::LinearSolver::symbolic: ";
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "The solver is NULL!  "
     "This should never happen.  Please report this bug to the Ifpack2 "
     "developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
     "set yet.  You must call setMatrix() with a nonnull matrix before you "
     "may call this method.");
  solver_->initialize ();
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
numeric ()
{
  const char prefix[] = "Ifpack2::Details::LinearSolver::numeric: ";
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "The solver is NULL!  "
     "This should never happen.  Please report this bug to the Ifpack2 "
     "developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The matrix has not been "
     "set yet.  You must call setMatrix() with a nonnull matrix before you "
     "may call this method.");
  solver_->compute ();
}

template<class SC, class LO, class GO, class NT>
std::string
LinearSolver<SC, LO, GO, NT>::
description () const
{
  const char prefix[] = "Ifpack2::Details::LinearSolver::description: ";
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "The solver is NULL!  "
     "This should never happen.  Please report this bug to the Ifpack2 "
     "developers.");
  return solver_->description ();
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  const char prefix[] = "Ifpack2::Details::LinearSolver::describe: ";
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::logic_error, prefix << "The solver is NULL!  "
     "This should never happen.  Please report this bug to the Ifpack2 "
     "developers.");
  solver_->describe (out, verbLevel);
}

} // namespace Details
} // namespace Ifpack2

// Explicit template instantiation macro for LinearSolver.  This is
// generally not for users!  It is used by automatically generated
// code, and perhaps by expert Trilinos developers.
#define IFPACK2_DETAILS_LINEARSOLVER_INSTANT(SC, LO, GO, NT) \
  template class Ifpack2::Details::LinearSolver<SC, LO, GO, NT>;

#endif // IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP
