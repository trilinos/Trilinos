//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#ifndef BELOS_TPETRA_SOLVER_MANAGER_BASE_HPP
#define BELOS_TPETRA_SOLVER_MANAGER_BASE_HPP

#include "BelosLinearProblem.hpp"
#include "BelosSolverManager.hpp"
#include "Belos_Tpetra_Krylov.hpp"

namespace BelosTpetra {
namespace Impl {

namespace ErrMsgs {
  static constexpr char linearProblemNotSet[] = "The linear problem has not "
    "yet been set.  Please call setProblem with a nonnull argument before "
    "calling this method.";
  static constexpr char solverNull[] = "The solver implementation (solver_) "
    "is null.  This should never happen.  "
    "Please report this bug to the Belos developers.";
}

template<class SC = Tpetra::MultiVector<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class SolverManagerBase :
    public Belos::SolverManager<SC, MV, OP>
{
public:
  using linear_problem_type = Belos::LinearProblem<SC, MV, OP>;
  
  SolverManagerBase () = delete;

  SolverManagerBase (const Teuchos::RCP<Krylov<SC, MV, OP>>& solver,
		     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) :
    solver_ (solver),
    lastSolverOutput_ {}
  {
    const char prefix[] = "SolverManagerBase::SolverManagerBase: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (solver.get () == nullptr, std::invalid_argument, prefix << ErrMsgs::solverNull);
    if (params.get () != nullptr) {
      this->setParameters (params);
    }
  }

  virtual ~SolverManagerBase () = default;

  const linear_problem_type& getProblem () const override {
    const char prefix[] = "SolverManagerBase::getProblem: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->problem_.get () == nullptr, std::runtime_error,
       prefix << ErrMsgs::linearProblemNotSet);
    return *(this->problem_);
  }

  //! Get valid parameters
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters () const override {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    const char prefix[] = "SolverManagerBase::getValidParameters: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (this->solver_.get () == nullptr, std::logic_error, prefix << ErrMsgs::solverNull);
    RCP<ParameterList> params (new ParameterList ("SolverManagerBase"));
    const bool defaultValues = true;
    this->solver_->getParameters (*params, defaultValues);
    return params;
  }

  //! Get current parameters
  Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters () const override {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    const char prefix[] = "SolverManagerBase::getCurrentParameters: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (this->solver_.get () == nullptr, std::logic_error, prefix << ErrMsgs::solverNull);
    RCP<ParameterList> params (new ParameterList ("SolverManagerBase"));
    const bool defaultValues = false;
    this->solver_->getParameters (*params, defaultValues);
    return params;
  }

  //! Get number of iterations
  int getNumIters () const override {
    return this->lastSolverOutput_.numIters;
  }

  //! Get number of restarts
  int getNumRests () const {
    return this->lastSolverOutput_.numRests;
  }

  bool isLOADetected () const override {
    return false; // this solver doesn't attempt to detect loss of accuracy
  }

  void
  setProblem (const Teuchos::RCP<linear_problem_type>& problem) override
  {
    const char prefix[] = "SolverManagerBase::setProblem: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->solver_.get () == nullptr, std::logic_error, prefix << ErrMsgs::solverNull);

    if (problem.is_null ()) {
      this->solver_->setMatrix (Teuchos::null);
    }
    else {
      if (this->solver_->getMatrix ().get () != problem->getOperator ().get ()) {
        // setMatrix resets state, so only call if necessary.
        this->solver_->setMatrix (problem->getOperator ());
      }
      if (problem->getRightPrec ().get () != nullptr && 
	  this->solver_->getPreconditioner () != problem->getRightPrec ()) {
        this->solver_->setRightPrec (problem->getRightPrec ());
      }
      if ((problem->getLeftPrec ()).get() != nullptr &&
          this->solver_->getPreconditioner () != problem->getLeftPrec ()) {
        this->solver_->setLeftPrec (problem->getLeftPrec ());
      }
    }
    this->problem_ = problem;
  }

  void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& params) override {
    const char prefix[] = "SolverManagerBase::setParameters: ";
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->solver_.get () == nullptr, std::logic_error, prefix << ErrMsgs::solverNull);
    if (! params.is_null ()) {
      this->solver_->setParameters (*params);
    }
  }

  void reset (const Belos::ResetType /* type */ ) override {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
  }

  Belos::ReturnType solve () override {
    using Teuchos::RCP;
    const char prefix[] = "SolverManagerBase::solve: ";

    TEUCHOS_TEST_FOR_EXCEPTION
      (this->solver_.get () == nullptr, std::logic_error,
       prefix << ErrMsgs::solverNull);
    TEUCHOS_TEST_FOR_EXCEPTION
      (this->problem_.is_null (), std::runtime_error,
       prefix << ErrMsgs::linearProblemNotSet);
    
    RCP<const MV> B = problem_->getRHS ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (B.is_null (), std::runtime_error, "The linear problem's right-hand "
       "side(s) B has/have not yet been set.  Please call setProblem with "
       "a nonnull argument before calling this method.");
    RCP<MV> X = problem_->getLHS ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.is_null (), std::runtime_error, "The linear problem's left-hand "
       "side(s) X has/have not yet been set.  Please call setProblem with "
       "a nonnull argument before calling this method.");

    lastSolverOutput_ = solver_->solve (*X, *B);
    return lastSolverOutput_.converged ? Belos::Converged : Belos::Unconverged;
  }

private:
  Teuchos::RCP<Krylov<SC, MV, OP>> solver_;
  //! Output of the last solve, not including the solution (multi)vector.
  SolverOutput<SC> lastSolverOutput_;
  Teuchos::RCP<linear_problem_type> problem_;
};

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_SOLVER_MANAGER_BASE_HPP

