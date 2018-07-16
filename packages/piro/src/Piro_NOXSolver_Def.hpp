// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_NOXSOLVER_DEF_HPP
#define PIRO_NOXSOLVER_DEF_HPP

#include "Piro_NOXSolver.hpp"
#include "Piro_MatrixFreeDecorator.hpp"

#include "Thyra_ModelEvaluatorHelpers.hpp"

#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

#include <stdexcept>
#include <cstddef>
#include <ostream>

template <typename Scalar>
Piro::NOXSolver<Scalar>::
NOXSolver(const Teuchos::RCP<Teuchos::ParameterList> &appParams_,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model_,
    const Teuchos::RCP<ObserverBase<Scalar> > &observer_) :
    SteadyStateSolver<Scalar>(model_),
    appParams(appParams_),
    observer(observer_),
    solver(new Thyra::NOXNonlinearSolver),
    out(Teuchos::VerboseObjectBase::getDefaultOStream()),
    model(model_),
    writeOnlyConvergedSol(appParams_->get("Write Only Converged Solution", true)),
    current_iteration(-1)
    {
  using Teuchos::RCP;

  const RCP<Teuchos::ParameterList> noxParams =
      Teuchos::sublist(appParams, "NOX", /*mustAlreadyExist =*/ false);

  solver->setParameterList(noxParams);

  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");

  exitUponFailedNOXSolve = appParams->get("Exit on Failed NOX Solve", false); 

  if (jacobianSource == "Matrix-Free") {
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model_,
          appParams->get<double>("Matrix-Free Perturbation")));
    }
    else 
      model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model_));
  }
  solver->setModel(model);
    }

template <typename Scalar>
void Piro::NOXSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
    {
  using Teuchos::RCP;
  bool observeFinalSolution = true;
  const int num_p = this->num_p();
  const int num_g = this->num_g();

  //For Analysis problems we typically do not want to write the solution at each NOX solver.
  //Instead we write at every "write interval" iterations of optimization solver.
  //If write_interval == -1, we print after every successful NOX solver
  //If write_inteval == 0 we ever print.
  //This relies on the fact that sensitivities are always called by ROL at each iteration to asses whether the solver is converged
  //TODO: when write_interval>1, at the moment there is no guarantee that the final iteration of the optimization (i.e. the converged solution) gets printed

  //When write_interval>0 we print only after computing the sensitivities,
  //to make sure to print them updated if required.
  bool solving_sensitivities = false;
  for (int i=0; i<num_p; i++) {
    for (int j=0; j<=num_g; j++) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() && !outArgs.get_DgDp(j,i).isEmpty()) {
        solving_sensitivities = true;
        break;
      }
    }
  }
  if(appParams->isSublist("Analysis")){
    auto analysisParams = appParams->sublist("Analysis");
    if(analysisParams.isSublist("Optimization Status")){
      auto optimizationParams = analysisParams.sublist("Optimization Status");
      if(optimizationParams.isParameter("Optimizer Iteration Number")) {
        int iteration = optimizationParams.get<int>("Optimizer Iteration Number");
        int write_interval = analysisParams.get("Write Interval",1);
        Teuchos::RCP<Teuchos::ParameterList> writeParams = Teuchos::rcp(new Teuchos::ParameterList());
        if (write_interval > 0) {
          observeFinalSolution = (iteration >= 0 && solving_sensitivities && iteration != current_iteration) ? (iteration%write_interval == 0) : false;
          if(observeFinalSolution)
            current_iteration = iteration;
        }
        else if (write_interval == 0)
          observeFinalSolution = false;
        else if (write_interval == -1)
          observeFinalSolution = true;
      }
    }
  }



  // Forward all parameters to underlying model
  Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = this->getModel().createInArgs();
  for (int l = 0; l < this->num_p(); ++l) {
    modelInArgs.set_p(l, inArgs.get_p(l));
  }

  // Find the solution of the implicit underlying model
  Thyra::SolveStatus<Scalar> solve_status;
  const Thyra::SolveCriteria<Scalar> solve_criteria;

  {
    solver->setBasePoint(modelInArgs);

    const RCP<const Thyra::VectorBase<Scalar> > modelNominalState =
        this->getModel().getNominalValues().get_x();

    const RCP<Thyra::VectorBase<Scalar> > initial_guess = modelNominalState->clone_v();

    solve_status = solver->solve(initial_guess.get(), &solve_criteria, /*delta =*/ NULL);

    //  MPerego: I think it is better not to throw an error when the solver does not converge.
    //  One can look at the solver status to check whether the solution is converged.
    if (exitUponFailedNOXSolve == true) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          solve_status.solveStatus != ::Thyra::SOLVE_STATUS_CONVERGED,
          std::runtime_error,
          "Nonlinear solver failed to converge");
    }
  }

  // Retrieve final solution to evaluate underlying model
  const RCP<const Thyra::VectorBase<Scalar> > finalSolution = solver->get_current_x();
  modelInArgs.set_x(finalSolution);

  if (Teuchos::nonnull(this->observer) && (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED || !writeOnlyConvergedSol) && observeFinalSolution) {
    this->observer->observeSolution(*finalSolution);
  }

  std::string sensitivity_method = appParams->get("Sensitivity Method",
      "Forward");

  //TODO: Handle forward and adjoint cases in a similar way.
  //      At the moment we keep the current implementation of Forward sensitivities
  //      and add code for adjoint sensitivities below.
  if(sensitivity_method == "Forward")
    return this->evalConvergedModel(modelInArgs, outArgs);
  else if(sensitivity_method != "Adjoint") {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
        Teuchos::Exceptions::InvalidParameter,
        std::endl <<
        "Piro::NOXSolver::evalModel():  " <<
        "Unknown sensitivity method" << sensitivity_method <<
        ".  Valid choices are \"Forward\" and \"Adjoint\"."
        << std::endl);
  }

  //Adjoint sensitivities:

  //
  // Do Sensitivity Calc, if requested. See 3 main steps
  //

  // Set inargs and outargs
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = this->getModel().createOutArgs();

  // Compute adjoint layouts of df/dp, dg/dx depending on
  bool do_sens = false;
  for (int i=0; i<num_p; i++) {
    // p
    modelInArgs.set_p(i, inArgs.get_p(i));

    // df/dp
    do_sens = false;
    for (int j=0; j<=num_g; j++) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() &&
          !outArgs.get_DgDp(j,i).isEmpty()) {
        do_sens = true;
      }
    }
    if (do_sens) {
      auto p_space = this->getModel().get_p_space(i);
      auto f_space = this->getModel().get_f_space();
      int num_params = p_space->dim();
      int num_resids = f_space->dim();
      auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(p_space);
      auto f_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(f_space);
      bool p_dist = !p_space_plus->isLocallyReplicated();//p_space->DistributedGlobal();
      bool f_dist = !f_space_plus->isLocallyReplicated();//f_space->DistributedGlobal();
      Thyra::ModelEvaluatorBase::DerivativeSupport ds =  modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp,i);
      // Determine which layout to use for df/dp.  Ideally one would look
      // at num_params, num_resids, what is supported by the underlying
      // model evaluator, and the sensitivity method, and make the best
      // choice to minimze the number of solves.  However this choice depends
      // also on what layout of dg/dx is supported (e.g., if only the operator
      // form is supported for forward sensitivities, then df/dp must be
      // DERIV_MV_BY_COL).  For simplicity, we order the conditional tests
      // to get the right layout in most situations.
      DerivativeLayout dfdp_layout;
      { // if (sensitivity_method == "Adjoint")
        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP))
          dfdp_layout = OP;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !f_dist)
          dfdp_layout = ROW;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !p_dist)
          dfdp_layout = COL;
        else
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              std::endl << "Piro::NOXSolver::evalModel():  " <<
              "For df/dp(" << i <<") with adjoint sensitivities, " <<
              "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
              "DERIV_MV_BY_COL with p not distributed, or "
              "DERIV_TRANS_MV_BY_ROW with f not distributed." <<
              std::endl);
      }
      if (dfdp_layout == COL) {
        auto dfdp = Thyra::createMembers(f_space, num_params);
        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        modelOutArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == ROW) {
        auto dfdp = Thyra::createMembers(p_space, num_resids);

        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
        dmv_dfdp(dfdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        modelOutArgs.set_DfDp(i,dmv_dfdp);
      }
      else if (dfdp_layout == OP) {
        auto dfdp_op = this->getModel().create_DfDp_op(i);
        TEUCHOS_TEST_FOR_EXCEPTION(
            dfdp_op == Teuchos::null, std::logic_error,
            std::endl << "Piro::NOXSolver::evalModel():  " <<
            "Needed df/dp operator (" << i << ") is null!" << std::endl);
        modelOutArgs.set_DfDp(i,dfdp_op);
      }
    }
  }

  for (int j=0; j<num_g; j++) {
    // g
    const RCP<Thyra::VectorBase<Scalar> > g_out = outArgs.get_g(j);
    if (g_out != Teuchos::null) {
      Thyra::put_scalar(0.0, g_out.ptr());
      modelOutArgs.set_g(j, g_out);
    }


    // dg/dx
    do_sens = false;
    for (int i=0; i<num_p; i++) {
      if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none() &&
          !outArgs.get_DgDp(j,i).isEmpty()) {
        do_sens = true;
      }
    }
    if (solving_sensitivities) {
      auto g_space = this->getModel().get_g_space(j);
      auto x_space = this->getModel().get_x_space();
      int num_responses = g_space->dim();
      int num_solution = x_space->dim();
      auto g_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(g_space);
      auto x_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(x_space);
      bool g_dist = !g_space_plus->isLocallyReplicated();//g_space->DistributedGlobal();
      bool x_dist = !x_space_plus->isLocallyReplicated();//x_space->DistributedGlobal();
      Thyra::ModelEvaluatorBase::DerivativeSupport ds =  modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx,j);
      DerivativeLayout dgdx_layout;
      { //if (sensitivity_method == "Adjoint")
        if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !g_dist)
          dgdx_layout = ROW;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !x_dist)
          dgdx_layout = COL;
        else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP))
          dgdx_layout = OP;
        else
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::logic_error,
              std::endl << "Piro::NOXSolver::evalModel():  " <<
              "For dg/dx(" << j <<") with adjoint sensitivities, " <<
              "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
              "DERIV_MV_BY_COL with x not distributed, or "
              "DERIV_TRANS_MV_BY_ROW with g not distributed." <<
              std::endl);
      }

      if (dgdx_layout == OP) {
        auto dgdx_op = this->getModel().create_DgDx_op(j);
        TEUCHOS_TEST_FOR_EXCEPTION(
            dgdx_op == Teuchos::null, std::logic_error,
            std::endl << "Piro::NOXSolver::evalModel():  " <<
            "Needed dg/dx operator (" << j << ") is null!" << std::endl);
        modelOutArgs.set_DgDx(j,dgdx_op);
      }
      else if (dgdx_layout == ROW) {
        auto dgdx = Thyra::createMembers(x_space, num_responses);
        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
        dmv_dgdx(dgdx, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
        modelOutArgs.set_DgDx(j,dmv_dgdx);
      }
      else if (dgdx_layout == COL) {
        auto dgdx = Thyra::createMembers(g_space, num_solution);
        Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
        dmv_dgdx(dgdx, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        modelOutArgs.set_DgDx(j,dmv_dgdx);
      }

      // dg/dp
      for (int i=0; i<num_p; i++) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,j,i).none()) {
          Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp = outArgs.get_DgDp(j,i);
          if (dgdp.getLinearOp() != Teuchos::null) {
            auto g_space = this->getModel().get_g_space(j);
            auto p_space = this->getModel().get_p_space(i);
            int num_responses = g_space->dim();
            int num_params = p_space->dim();
            auto g_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(g_space);
            auto p_space_plus = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(p_space);
            bool g_dist = !g_space_plus->isLocallyReplicated();//g_space->DistributedGlobal();
            bool p_dist = !p_space_plus->isLocallyReplicated();//p_space->DistributedGlobal();
            Thyra::ModelEvaluatorBase::DerivativeSupport ds = modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,j,i);
            if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
              auto dgdp_op =
                  this->getModel().create_DgDp_op(j,i);
              TEUCHOS_TEST_FOR_EXCEPTION(
                  dgdp_op == Teuchos::null, std::logic_error,
                  std::endl << "Piro::NOXSolver::evalModel():  " <<
                  "Needed dg/dp operator (" << j << "," << i << ") is null!" <<
                  std::endl);
              modelOutArgs.set_DgDp(j,i,dgdp_op);
            }
            else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) && !p_dist) {
              auto dgdp = createMembers(g_space, num_params);
              Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
              dmv_dgdp(dgdp, Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
              modelOutArgs.set_DgDp(j,i,dmv_dgdp);
            }
            else if (ds.supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) && !g_dist) {
              auto dgdp = createMembers(p_space, num_responses);
              Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar>
              dmv_dgdp(dgdp, Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
              modelOutArgs.set_DgDp(j,i,dmv_dgdp);
            }
            else
              TEUCHOS_TEST_FOR_EXCEPTION(
                  true, std::logic_error,
                  std::endl << "Piro::NOXSolver::evalModel():  " <<
                  "For dg/dp(" << j << "," << i <<
                  ") with operator sensitivities, "<<
                  "underlying ModelEvaluator must support DERIV_LINEAR_OP, " <<
                  "DERIV_MV_BY_COL with p not distributed, or "
                  "DERIV_TRANS_MV_BY_ROW with g not distributed." <<
                  std::endl);
          }
          else
            modelOutArgs.set_DgDp(j,i,outArgs.get_DgDp(j,i));
        }
      }
    }
  }

  // (1) Calculate g, df/dp, dg/dp, dg/dx
  this->getModel().evalModel(modelInArgs, modelOutArgs);


  // Ensure Jacobian is up-to-date
  if (solving_sensitivities) // Adjoint sensitivities
  {

    // Create implicitly transpose Jacobian and preconditioner
    Teuchos::RCP< const Thyra::LinearOpWithSolveFactoryBase<double> > lows_factory = model->get_W_factory();
    TEUCHOS_ASSERT(Teuchos::nonnull(lows_factory));
    Teuchos::RCP< Thyra::LinearOpBase<double> > lop = model->create_W_op();
    Teuchos::RCP< const ::Thyra::DefaultLinearOpSource<double> > losb = Teuchos::rcp(new ::Thyra::DefaultLinearOpSource<double>(lop));
    Teuchos::RCP< ::Thyra::PreconditionerBase<double> > prec;

    auto in_args = model->createInArgs();
    auto out_args = model->createOutArgs();

    Teuchos::RCP< ::Thyra::PreconditionerFactoryBase<double> > prec_factory =  lows_factory->getPreconditionerFactory();
    if (Teuchos::nonnull(prec_factory)) {
      prec = prec_factory->createPrec();
    } else if (out_args.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) {
      prec = model->create_W_prec();
    }

    const RCP<Thyra::LinearOpWithSolveBase<Scalar> > jacobian =
        lows_factory->createOp();

    {
      in_args.set_x(finalSolution);
      out_args.set_W_op(lop);
      model->evalModel(in_args, out_args);
      in_args.set_x(Teuchos::null);
      out_args.set_W_op(Teuchos::null);
    }

    if (Teuchos::nonnull(prec_factory))
      prec_factory->initializePrec(losb, prec.get());
    else if ( Teuchos::nonnull(prec) && (out_args.supports( Thyra::ModelEvaluatorBase::OUT_ARG_W_prec)) ) {
      in_args.set_x(finalSolution);
      out_args.set_W_prec(prec);
      model->evalModel(in_args, out_args);
    }

    if(Teuchos::nonnull(prec))
      Thyra::initializePreconditionedOp<double>(*lows_factory,
          Thyra::transpose<double>(lop),
          Thyra::unspecifiedPrec<double>(::Thyra::transpose<double>(prec->getUnspecifiedPrecOp())),
          jacobian.ptr());
    else
      Thyra::initializeOp<double>(*lows_factory,
          Thyra::transpose<double>(lop),
          jacobian.ptr());


    for (int j=0; j<num_g; j++) {

      // See if there are any forward sensitivities we need to do
      // that aren't handled by the operator
      do_sens = false;
      for (int i=0; i<num_p; i++) {
        if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none()) {
          if (outArgs.get_DgDp(j,i).getMultiVector() != Teuchos::null)
            do_sens = true;
        }
      }
      if (!do_sens)
        continue;

      if (!modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j).none()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
            modelOutArgs.get_DgDx(j).getLinearOp()!=Teuchos::null,
            std::logic_error,
            std::endl << "Piro::NOXSolver::evalModel():  " <<
            "Can\'t use dg/dx operator " << j << " with non-operator " <<
            "adjoint sensitivities." << std::endl);
        auto dgdx  = modelOutArgs.get_DgDx(j).getMultiVector();
        if (dgdx != Teuchos::null) {
          int num_cols = dgdx->domain()->dim();

          // (2) Calculate xbar multivector from -(J^{-T}*dg/dx)

          auto xbar = createMembers(dgdx->range(), num_cols);
          xbar->assign(0);
          Thyra::solve(
              *jacobian,
              Thyra::NOTRANS,
              *dgdx,
              xbar.ptr(),
              Teuchos::ptr(&solve_criteria));

          Thyra::scale(-1.0, xbar.ptr());

          // (3) Calculate dg/dp^T = df/dp^T*xbar + dg/dp^T
          for (int i=0; i<num_p; i++) {
            if (!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, i).none()) {
              auto dgdp_out = outArgs.get_DgDp(j,i).getMultiVector();
              if (dgdp_out != Teuchos::null) {
                Thyra::ModelEvaluatorBase::Derivative<Scalar> dfdp_dv = modelOutArgs.get_DfDp(i);
                auto dfdp_op = dfdp_dv.getLinearOp();
                auto dfdp = dfdp_dv.getMultiVector();
                Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_dv = modelOutArgs.get_DgDx(j);
                Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient =
                    outArgs.get_DgDp(j,i).getMultiVectorOrientation();
                if (dfdp_op != Teuchos::null) {
                  bool transpose = false;
                  if (dgdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)
                    transpose = true;
                  auto tmp = createMembers(dfdp_op->domain(),
                      xbar->domain()->dim());

                  dfdp_op->apply(Thyra::TRANS,*xbar, tmp.ptr(),1.0, 0.0);
                  auto dgdp_range = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceDefaultBase<Scalar>>(dgdp_out->range());

                  if (transpose) {
                    TEUCHOS_TEST_FOR_EXCEPTION(
                        !dgdp_range->isLocallyReplicated(),
                        //dgdp_out->Map().DistributedGlobal(),
                        std::logic_error,
                        std::endl <<
                        "Piro::NOXSolver::evalModel():  " <<
                        "Can\'t handle special case:  " <<
                        " df/dp operator, " <<
                        " transposed, distributed dg/dp. " << std::endl);
                    Thyra::DetachedMultiVectorView<Scalar> dgdp_out_view(dgdp_out);
                    Thyra::DetachedMultiVectorView<Scalar> tmp_view(tmp);
                    for (int j=0; j<dgdp_out_view.numSubCols(); j++)
                      for (int i=0; i<dgdp_out_view.subDim(); i++)
                        dgdp_out_view(i,j) += tmp_view(j,i);
                  }
                  else {
                    Thyra::update(1.0,  *tmp, dgdp_out.ptr());
                  }
                }
                else {
                  Teuchos::RCP<Thyra::MultiVectorBase<Scalar>> arg1, arg2;
                  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient =
                      modelOutArgs.get_DgDp(j,i).getMultiVectorOrientation();
                  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dfdp_orient =
                      modelOutArgs.get_DfDp(i).getMultiVectorOrientation();
                  Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdx_orient =
                      modelOutArgs.get_DgDx(j).getMultiVectorOrientation();


                  Thyra::DetachedMultiVectorView<Scalar> dgdp_out_view(dgdp_out);
                  int num_g = (dgdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) ?
                      dgdp_out_view.subDim() : dgdp_out_view.numSubCols();
                  int num_p = (dgdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) ?
                      dgdp_out_view.numSubCols() : dgdp_out_view.subDim();
                  if(dgdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM){
                    //dgdp (np columns of g_space vectors)

                    if (dgdx_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                      Thyra::DetachedMultiVectorView<Scalar> xbar_view(xbar);
                      Thyra::DetachedMultiVectorView<Scalar> dfdp_view(dfdp);
                      if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                        for (int ip=0; ip<num_p; ip++)
                          for (int ig=0; i<num_g; ig++) {
                            for (int ix=0; ix<xbar_view.numSubCols(); ix++)
                              dgdp_out_view(ip,ig) += xbar_view(ix,ig)*dfdp_view(ip,ix);
                          }
                      }
                      else {
                        for (int ip=0; ip<num_p; ip++)
                          for (int ig=0; ig<num_g; ig++)
                            for (int ix=0; ix<xbar_view.numSubCols(); ix++)
                              dgdp_out_view(ip,ig) += xbar_view(ix,ig)*dfdp_view(ix,ip);
                      }
                    }
                    else if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                      for (int ip=0; ip<num_p; ip++)
                        for (int ig=0; ig<num_g; ig++)
                          dgdp_out_view(ip,ig) += Thyra::scalarProd(*xbar->col(ig),*dfdp->col(ip));
                    }
                    else {
                      Thyra::DetachedMultiVectorView<Scalar> xbar_view(xbar);
                      Thyra::DetachedMultiVectorView<Scalar> dfdp_view(dfdp);
                      for (int ip=0; ip<dgdp_out_view.numSubCols(); ip++)
                        for (int ig=0; ig<dgdp_out_view.subDim(); ig++)
                          for (int ix=0; ix<xbar_view.subDim(); ix++)
                            dgdp_out_view(ip,ig) += xbar_view(ig,ix)*dfdp_view(ix,ip);
                    }
                  }
                  else {
                    //dgdp (ng columns of p_space vectors)
                    if (dgdx_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                      Thyra::DetachedMultiVectorView<Scalar> xbar_view(xbar);
                      Thyra::DetachedMultiVectorView<Scalar> dfdp_view(dfdp);
                      if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                        for (int ip=0; ip<num_p; ip++)
                          for (int ig=0; i<num_g; ig++) {
                            for (int ix=0; ix<xbar_view.numSubCols(); ix++)
                              dgdp_out_view(ig,ip) += xbar_view(ix,ig)*dfdp_view(ip,ix);
                          }
                      }
                      else {
                        for (int ip=0; ip<num_p; ip++)
                          for (int ig=0; ig<num_g; ig++)
                            for (int ix=0; ix<xbar_view.numSubCols(); ix++)
                              dgdp_out_view(ig,ip) += xbar_view(ix,ig)*dfdp_view(ix,ip);
                      }
                    }
                    else if (dfdp_orient == Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) {
                      for (int ip=0; ip<num_p; ip++)
                        for (int ig=0; ig<num_g; ig++)
                          dgdp_out_view(ig,ip) += Thyra::scalarProd(*xbar->col(ig),*dfdp->col(ip));
                    }
                    else {
                      Thyra::DetachedMultiVectorView<Scalar> xbar_view(xbar);
                      Thyra::DetachedMultiVectorView<Scalar> dfdp_view(dfdp);
                      for (int ip=0; ip<dgdp_out_view.numSubCols(); ip++)
                        for (int ig=0; ig<dgdp_out_view.subDim(); ig++)
                          for (int ix=0; ix<xbar_view.subDim(); ix++)
                            dgdp_out_view(ig,ip) += xbar_view(ig,ix)*dfdp_view(ix,ip);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

    }

#endif /*PIRO_NOXSOLVER_DEF_HPP*/
