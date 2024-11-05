// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_DETAILS_LINEARSOLVER_HPP
#define BELOS_DETAILS_LINEARSOLVER_HPP

/// \file Belos_Details_LinearSolver.hpp
/// \brief Implementation of Trilinos::Details::LinearSolver.

#include "BelosSolverFactory.hpp"
#include "Trilinos_Details_LinearSolver.hpp"

namespace Belos {
namespace Details {

/// \class LinearSolver
/// \brief Belos' implementation of Trilinos::Details::LinearSolver.
///
/// Note to developers: The main question for developers writing a
/// package's implementation of LinearSolver, is whether it can just
/// wrap the package's solver interface directly.  (This means that
/// the LinearSolver subclass' constructor just takes an RCP to the
/// package's solver.)  If not, then the LinearSolver needs a
/// "factory" inside that can create solvers from that package.  This
/// complicates the code a bit and may add to build time.
///
/// Belos lets users create a SolverManager without a matrix, and
/// change the matrix after creation.  However, Belos solvers can't
/// handle changes to the matrix's domain or range Maps.  (Supposedly
/// they should be able to do this, but implementations of
/// SolverManager::reset tend just to call setProblem() without
/// reallocating the basis.)  It's safest in this case to destroy the
/// solver and start over.  SolverManager instances don't know how to
/// clone themselves in an uninitialized state.  This means that
/// LinearSolver cannot wrap the SolverManager directly; it must be
/// able to create Belos solvers inside.
template<class MV, class OP, class ScalarType, class NormType>
class LinearSolver :
    public Trilinos::Details::LinearSolver<MV, OP, NormType>
{
private:
  //! Belos::LinearProblem specialization used by this class.
  typedef Belos::LinearProblem<ScalarType, MV, OP> problem_type;
  //! Belos' own solver type.
  typedef Belos::SolverManager<ScalarType, MV, OP> solver_type;

public:

  /// \brief Constructor
  ///
  /// \param solverName The name of the Belos::SolverManager instance
  ///   to create.
  LinearSolver (const std::string& solverName) :
    solverName_ (solverName)
  {
    // In alignment with Belos philosophy, we delay initialization
    // (which in this case means creation of the solver) until needed.
  }


  //! Destructor (virtual for memory safety).
  virtual ~LinearSolver () {}

  /// \brief Set the solver's matrix.
  ///
  /// \param A [in] Pointer to the matrix A in the linear system(s)
  ///   AX=B to solve.
  void setMatrix (const Teuchos::RCP<const OP>& A) {
    if (problem_.is_null ()) {
      problem_ = Teuchos::rcp (new problem_type (A, Teuchos::null, Teuchos::null));
    } else if (A != problem_->getOperator ()) {
      problem_->setOperator (A);
    }
  }

  //! Get the solver's matrix.
  Teuchos::RCP<const OP> getMatrix () const {
    if (problem_.is_null ()) {
      return Teuchos::null;
    } else {
      return problem_->getOperator ();
    }
  }

  //! Solve the linear system AX=B for X.
  void solve (MV& X, const MV& B) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (problem_.is_null () || problem_->getOperator ().is_null (),
       std::runtime_error, "Belos::Details::LinearSolver::solve: "
       "The matrix A in the linear system to solve has not yet been set.  "
       "Please call setMatrix() with a nonnull input before calling solve().");
    Teuchos::RCP<MV> X_ptr = Teuchos::rcpFromRef (X);
    Teuchos::RCP<const MV> B_ptr = Teuchos::rcpFromRef (B);

    problem_->setLHS (X_ptr);
    problem_->setRHS (B_ptr);
    problem_->setProblem ();

    // We can delay creating the Belos solver until the moment when we
    // actually need it.  This aligns with Belos' preference for lazy
    // initialization.
    if (solver_.is_null ()) {
      Belos::SolverFactory<ScalarType, MV, OP> factory;
      solver_ = factory.create (solverName_, params_);
      solver_->setProblem (problem_);
    }

    //Belos::ReturnType ret = solver_->solve ();
    (void) solver_->solve ();
  }

  //! Set the solver's parameters.
  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) {
    if (! solver_.is_null () && ! params.is_null ()) {
      solver_->setParameters (params);
    }
    params_ = params;
  }

  //! Precompute for matrix structure changes.
  void symbolic () {
    TEUCHOS_TEST_FOR_EXCEPTION
      (problem_.is_null () || problem_->getOperator ().is_null (),
       std::runtime_error, "Belos::Details::LinearSolver::symbolic: "
       "The matrix A in the linear system to solve has not yet been set.  "
       "Please call setMatrix() with a nonnull input before calling this method.");

    // Belos solvers can't handle changes to the matrix's domain or
    // range Maps.  It's best in this case to destroy the solver and
    // start over.
    solver_ = Teuchos::null;
  }

  //! Precompute for matrix values' changes.
  void numeric () {
    TEUCHOS_TEST_FOR_EXCEPTION
      (problem_.is_null () || problem_->getOperator ().is_null (),
       std::runtime_error, "Belos::Details::LinearSolver::numeric: "
       "The matrix A in the linear system to solve has not yet been set.  "
       "Please call setMatrix() with a nonnull input before calling this method.");
    // NOTE (mfh 23 Aug 2015) For the seed or recycling solvers, it
    // would make sense to do something special here.  However, the
    // line below is always correct.  It recomputes the initial
    // residual vector, which is what Belos expects before a solve if
    // the matrix or right-hand side may have changed.
    problem_->setProblem ();
  }

private:
  //! The name of the Belos solver to create.
  std::string solverName_;
  //! The LinearProblem instance to give to the Belos solver.
  Teuchos::RCP<problem_type> problem_;
  //! The Belos solver (SolverManager instance).
  Teuchos::RCP<solver_type> solver_;
  //! The Belos solver's list of parameters.
  Teuchos::RCP<Teuchos::ParameterList> params_;
};

} // namespace Details
} // namespace Belos

#endif /* BELOS_DETAILS_LINEARSOLVER_HPP */
