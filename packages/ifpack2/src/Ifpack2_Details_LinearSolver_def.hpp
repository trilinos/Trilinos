#ifndef IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP
#define IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP

#include <Ifpack2_Factory.hpp>

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
LinearSolver (const std::string& solverName) :
  solverName_ (solverName)
{}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
setMatrix (const Teuchos::RCP<const OP>& A)
{
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  typedef row_matrix_type MAT;
  const char prefix[] = "Ifpack2::Details::LinearSolver::setMatrix: ";

  // This needs to invoke Ifpack2's (internal) factory, in order to
  // create the solver on demand if necessary.  That's because
  // Ifpack2::Preconditioner has no setMatrix() method.
  if (A.is_null ()) {
    // Setting the matrix to null is a way for users to signal that
    // they want to free up the solver's resources.
    solver_ = null;
  }
  else if (solver_.is_null () || A_ != A) {
    // Ifpack2 solvers expect a RowMatrix.
    RCP<const MAT> A_mat = Teuchos::rcp_dynamic_cast<const MAT> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_mat.is_null (), std::invalid_argument, prefix <<
       "The input Tpetra::Operator A must be a Tpetra::RowMatrix.");

    // Ifpack2::Preconditioner instances must be created with a
    // nonnull matrix.  Thus, we don't actually create the solver
    // until setMatrix is called for the first time with a nonnull
    // matrix.  Furthermore, Ifpack2::Preconditioner has no
    // setMatrix method.  Thus, if we change the matrix, we have to
    // free and recreate the solver.  We risk extra memory usage in
    // favor of the strong exception guarantee.
    RCP<prec_type> solver;
    try {
      solver = Ifpack2::Factory::template create<MAT> (solverName_, A_mat);
    }
    catch (std::exception& e) {
      TEUCHOS_TEST_FOR_EXCEPTION
        (true, std::invalid_argument, prefix << "Failed to create Ifpack2 "
         "preconditioner named \"" << solverName_ << "\", for the following "
         "template parameters: "
         << "SC = " << TypeNameTraits<SC>::name ()
         << ", LO = " << TypeNameTraits<LO>::name ()
         << ", GO = " << TypeNameTraits<GO>::name ()
         << ", NT = " << TypeNameTraits<NT>::name ()
         << ".  Ifpack2::Factory::create threw an exception: " << e.what ());
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver.is_null (), std::invalid_argument, prefix << "Failed to create "
       "Ifpack2 preconditioner named \"" << solverName_ << "\", for the "
       "following template parameters: "
       << "SC = " << TypeNameTraits<SC>::name ()
       << ", LO = " << TypeNameTraits<LO>::name ()
       << ", GO = " << TypeNameTraits<GO>::name ()
       << ", NT = " << TypeNameTraits<NT>::name ()
       << ".  Ifpack2::Factory::create returned null.");

    // Use same parameters as before, if user set parameters.
    if (! params_.is_null ()) {
      // Ifpack2::Preconditioner expects a const ParameterList&.
      solver->setParameters (*params_);
    }
    solver_ = solver;
  }

  if (A_ != A) {
    A_ = A; // keep a pointer to A, so that getMatrix() works
  }
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
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver"
     "::solve: The solver does not exist yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver::"
     "solve: The matrix has not been set yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  solver_->apply (B, X);
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  if (! solver_.is_null ()) {
    solver_->setParameters (*params);
  }
  // Remember them, so that if the solver gets recreated, we'll have
  // the original parameters.
  params_ = params;
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
symbolic ()
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver"
     "::symbolic: The solver does not exist yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver::"
     "symbolic: The matrix has not been set yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  solver_->initialize ();
}

template<class SC, class LO, class GO, class NT>
void
LinearSolver<SC, LO, GO, NT>::
numeric ()
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver"
     "::numeric: The solver does not exist yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, "Ifpack2::Details::LinearSolver::"
     "numeric: The matrix has not been set yet.  You must call setMatrix() "
     "with a nonnull matrix before you may call this method.");
  solver_->compute ();
}

} // namespace Details
} // namespace Ifpack2

// Explicit template instantiation macro for LinearSolver.  This is
// generally not for users!  It is used by automatically generated
// code, and perhaps by expert Trilinos developers.
#define IFPACK2_DETAILS_LINEARSOLVER_INSTANT(SC, LO, GO, NT) \
  template class Ifpack2::Details::LinearSolver<SC, LO, GO, NT>;

#endif // IFPACK2_DETAILS_LINEARSOLVER_DEF_HPP
