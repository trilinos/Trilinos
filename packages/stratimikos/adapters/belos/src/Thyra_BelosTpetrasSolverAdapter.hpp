// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef __Thyra_BelosTpetraSolverAdapter_hpp
#define __Thyra_BelosTpetraSolverAdapter_hpp

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
// BelosThyraAdapter.hpp only includes this file if HAVE_BELOS_TSQR is
// defined.  Thus, it's OK to include TSQR header files here.

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "Tpetra_CrsMatrix.hpp"
#include "BelosTpetraAdapter.hpp"
#include "solvers/Belos_Tpetra_Krylov_parameters.hpp"
#include "solvers/Belos_Tpetra_Krylov.hpp"
#include "solvers/Belos_Tpetra_Gmres.hpp"
#include "solvers/Belos_Tpetra_GmresPipeline.hpp"
#include "solvers/Belos_Tpetra_GmresSingleReduce.hpp"
#include "solvers/Belos_Tpetra_GmresSstep.hpp"

namespace Thyra {

  /// \class base Krylov
  ///
  template<class SC, class MV, class OP>
  class BelosTpetraKrylov : public Belos::SolverManager<SC, MV, OP> {
  public:
    using tpetra_base_solver_type = BelosTpetra::Impl::Krylov<SC>;
    using converter = Thyra::TpetraOperatorVectorExtraction<SC>;

    //! constructor
    BelosTpetraKrylov() = default;

    //! clone for Inverted Injection (DII)
    virtual Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const = 0;

    //! set/get problem
    void setProblem( const Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> > &problem ) override {
      problem_ = problem;
    }

    const Belos::LinearProblem<SC, MV, OP>& getProblem() const override {
      return *problem_;
    }

    //! set/get parameters
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList> &params ) override {
      params_ = params;
      tpetra_solver->setParameters(*params);
    }

    Teuchos::RCP<const Teuchos::ParameterList> getCurrentParameters() const override {
      return params_;
    }

    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override {

      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::parameterList();

      bool defaultValues = true;
      tpetra_solver->getParameters (*params, defaultValues);

      return params;
    }

    //! Get the iteration count for the most recent call to \c solve().
    int  getNumIters() const override { return solver_output.numIters; }
    bool isLOADetected() const override { return false; }
    void reset( const Belos::ResetType type ) override {
      if ((type & Belos::Problem) && !Teuchos::is_null(problem_))
        problem_->setProblem();
    }

    //! solve
    Belos::ReturnType solve() override {
      auto  A = problem_->getOperator();
      auto tpetraA = converter::getConstTpetraOperator(A);
      tpetra_solver->setMatrix(tpetraA);

      auto lM = problem_->getLeftPrec();
      if (lM != Teuchos::null) {
        auto tpetraLM = converter::getConstTpetraOperator(lM);
        tpetra_solver->setLeftPrec(tpetraLM);
      }
      auto rM = problem_->getRightPrec();
      if (rM != Teuchos::null) {
        auto tpetraRM = converter::getConstTpetraOperator(rM);
        tpetra_solver->setRightPrec(tpetraRM);
      }

      auto  B = this->problem_->getRHS();
      auto  X = this->problem_->getLHS();
      auto tpetraB = converter::getConstTpetraMultiVector(B);
      auto tpetraX = converter::getTpetraMultiVector(X);

      solver_output = tpetra_solver->solve(*tpetraX, *tpetraB);

      // copy output
      Belos::ReturnType belos_output = (solver_output.converged ? Belos::Converged : Belos::Unconverged);
      return belos_output;
    }

  protected:

    Teuchos::RCP<tpetra_base_solver_type> tpetra_solver;
    BelosTpetra::Impl::SolverOutput<SC> solver_output;

    //! Current parameter list.
    Teuchos::RCP<Teuchos::ParameterList> params_;

    //! The linear problem to solve.
    Teuchos::RCP<Belos::LinearProblem<SC, MV, OP> > problem_;
  };

  //
  // Gmres
  template<class SC, class MV, class OP>
  class BelosTpetraGmres : public BelosTpetraKrylov<SC, MV, OP> {
  public:
    using krylov_base_solver_type = BelosTpetraKrylov<SC, MV, OP>;
    using tpetra_solver_type = BelosTpetra::Impl::Gmres<SC>;

    //! constructor
    BelosTpetraGmres() :
    krylov_base_solver_type ()
    {
      this->tpetra_solver = Teuchos::rcp(new tpetra_solver_type ());
    }

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const override {
      return Teuchos::rcp(new BelosTpetraGmres<SC, MV, OP>);
    }
  };

  //
  // Pipelined Gmres 
  template<class SC, class MV, class OP>
  class BelosTpetraGmresPipeline : public BelosTpetraKrylov<SC, MV, OP> {
  public:
    using krylov_base_solver_type = BelosTpetraKrylov<SC, MV, OP>;
    using tpetra_solver_type = BelosTpetra::Impl::GmresPipeline<SC>;

    //! constructor
    BelosTpetraGmresPipeline() :
    krylov_base_solver_type ()
    {
      this->tpetra_solver = Teuchos::rcp(new tpetra_solver_type ());
    }

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const override {
      return Teuchos::rcp(new BelosTpetraGmresPipeline<SC, MV, OP>);
    }
  };

  //
  // SingleReduce Gmres 
  template<class SC, class MV, class OP>
  class BelosTpetraGmresSingleReduce : public BelosTpetraKrylov<SC, MV, OP> {
  public:
    using krylov_base_solver_type = BelosTpetraKrylov<SC, MV, OP>;
    using tpetra_solver_type = BelosTpetra::Impl::GmresSingleReduce<SC>;

    //! constructor
    BelosTpetraGmresSingleReduce() :
    krylov_base_solver_type ()
    {
      this->tpetra_solver = Teuchos::rcp(new tpetra_solver_type ());
    }

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const override {
      return Teuchos::rcp(new BelosTpetraGmresSingleReduce<SC, MV, OP>);
    }
  };

  //
  // s-step Gmres 
  template<class SC, class MV, class OP>
  class BelosTpetraGmresSstep : public BelosTpetraKrylov<SC, MV, OP> {
  public:
    using krylov_base_solver_type = BelosTpetraKrylov<SC, MV, OP>;
    using tpetra_solver_type = BelosTpetra::Impl::GmresSstep<SC>;

    //! constructor
    BelosTpetraGmresSstep() :
    krylov_base_solver_type ()
    {
      this->tpetra_solver = Teuchos::rcp(new tpetra_solver_type ());
    }

    //! clone for Inverted Injection (DII)
    Teuchos::RCP<Belos::SolverManager<SC, MV, OP> > clone () const override {
      return Teuchos::rcp(new BelosTpetraGmresSstep<SC, MV, OP>);
    }
  };

} // namespace Thyra

#endif // __Thyra_TpetraSolverAdapter_hpp

