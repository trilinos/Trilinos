// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_LOCAADAPTIVESOLVER_DEF_HPP
#define PIRO_LOCAADAPTIVESOLVER_DEF_HPP

#include "Piro_LOCAAdaptiveSolver.hpp"

#include "Piro_ObserverToLOCASaveDataStrategyAdapter.hpp"

#include "Piro_MatrixFreeDecorator.hpp"

#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

#include "NOX_StatusTest_Factory.H"

#include "Teuchos_as.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"

#include <stdexcept>
#include <ostream>

namespace Piro {

namespace Detail {

class ModelEvaluatorParamName {
public:
  explicit ModelEvaluatorParamName(const Teuchos::RCP<const Teuchos::Array<std::string> > &p_names);
  std::string operator()(Teuchos_Ordinal k) const;

private:
  Teuchos::RCP<const Teuchos::Array<std::string> > p_names_;
  enum { Default, OneShared, FullList } type_;
};

} // namespace Detail

} // namespace Piro


template <typename Scalar>
Piro::LOCAAdaptiveSolver<Scalar>::LOCAAdaptiveSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy) :
  SteadyStateSolver<Scalar>(model, adjointModel, model->Np() > 0), // Only one parameter supported
  piroParams_(piroParams),
  saveDataStrategy_(saveDataStrategy),
  globalData_(LOCA::createGlobalData(piroParams)),
  paramVector_(),
  solMgr_(solMgr),
  locaStatusTests_(),
  noxStatusTests_(),
  stepper_(),
  model_(model)
{
  const int l = 0; // TODO: Allow user to select parameter index
  const Detail::ModelEvaluatorParamName paramName(this->getModel().get_p_names(l));
  const Thyra::Ordinal p_entry_count = this->getModel().get_p_space(l)->dim();
  for (Teuchos_Ordinal k = 0; k < p_entry_count; ++k) {
    (void) paramVector_.addParameter(paramName(k));
  }
  
  std::string jacobianSource = piroParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    if (piroParams->isParameter("Matrix-Free Perturbation")) {
      model_ = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model,
                           piroParams->get<double>("Matrix-Free Perturbation")));
    }
    else model_ = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(model));
  }

  solMgr_->initialize(Teuchos::rcp(new Thyra::LOCAAdaptiveState(model_, saveDataStrategy_, globalData_, 
            Teuchos::rcpFromRef(paramVector_), l)));

  // TODO: Create non-trivial stopping criterion for the stepper
  locaStatusTests_ = Teuchos::null;

  // Create stopping criterion for the nonlinear solver
  const Teuchos::RCP<Teuchos::ParameterList> noxStatusParams =
    Teuchos::sublist(Teuchos::sublist(piroParams_, "NOX"), "Status Tests");
  noxStatusTests_ = NOX::StatusTest::buildStatusTests(*noxStatusParams, *(globalData_->locaUtils));

  stepper_ = Teuchos::rcp(new LOCA::AdaptiveStepper(piroParams_, solMgr_, globalData_, noxStatusTests_));

  if (piroParams_->isSublist("NOX") &&
      piroParams_->sublist("NOX").isSublist("Printing"))
    utils_.reset(piroParams_->sublist("NOX").sublist("Printing"));

  this->setSensitivityMethod("Forward");
}

template<typename Scalar>
Piro::LOCAAdaptiveSolver<Scalar>::~LOCAAdaptiveSolver()
{
  LOCA::destroyGlobalData(globalData_);
}

template <typename Scalar>
void
Piro::LOCAAdaptiveSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  const int l = 0; // TODO: Allow user to select parameter index
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > p_inargs = inArgs.get_p(l);

  // Forward parameter values to the LOCAAdaptive stepper
  {
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > p_inargs_or_nominal =
      Teuchos::nonnull(p_inargs) ? p_inargs : this->getNominalValues().get_p(l);
    const Thyra::ConstDetachedVectorView<Scalar> p_init_values(p_inargs_or_nominal);
    const Teuchos_Ordinal p_entry_count = p_init_values.subDim();
    TEUCHOS_ASSERT(p_entry_count == Teuchos::as<Teuchos_Ordinal>(paramVector_.length()));

    for (Teuchos_Ordinal k = 0; k < p_entry_count; ++k) {
      paramVector_[k] = p_init_values[k];
    }

//    solMgr_->getSolutionGroup()->setParams(paramVector_);
    Teuchos::rcp_dynamic_cast< ::Thyra::LOCAAdaptiveState >(solMgr_->getState())
                 ->getSolutionGroup()->setParams(paramVector_);
  }

  LOCA::Abstract::Iterator::IteratorStatus status;

  status = stepper_->run();

  if (status == LOCA::Abstract::Iterator::Finished) {
    utils_.out() << "Continuation Stepper Finished.\n";
  } else if (status == LOCA::Abstract::Iterator::NotFinished) {
    utils_.out() << "Continuation Stepper did not reach final value.\n";
  } else {
    utils_.out() << "Nonlinear solver failed to converge.\n";
    outArgs.setFailed();
  }

  // The time spent
  globalData_->locaUtils->out() << std::endl <<
    "#### Statistics ########" << std::endl;

  // Check number of steps
  int numSteps = stepper_->getStepNumber();
  globalData_->locaUtils->out() << std::endl <<
    " Number of continuation Steps = " << numSteps << std::endl;

  // Check number of failed steps
  int numFailedSteps = stepper_->getNumFailedSteps();
  globalData_->locaUtils->out() << std::endl <<
    " Number of failed continuation Steps = " << numFailedSteps << std::endl;

  globalData_->locaUtils->out() << std::endl;


  // Note: the last g is used to store the final solution. It can be null - if it is just
  // skip the store. If adaptation has occurred, g is not the correct size.

  const Teuchos::RCP<Thyra::VectorBase<Scalar> > x_outargs = outArgs.get_g(this->num_g());
  Teuchos::RCP<Thyra::VectorBase<Scalar> > x_final;

  int x_args_dim = 0;
  int f_sol_dim = 0;

  // Pardon the nasty cast to resize the last g in outArgs - need to fit the solution
  Thyra::ModelEvaluatorBase::OutArgs<Scalar>* mutable_outArgsPtr =
    const_cast<Thyra::ModelEvaluatorBase::OutArgs<Scalar>* >(&outArgs);

  if(Teuchos::nonnull(x_outargs)){ // g has been allocated, calculate the sizes of g and the solution

    x_args_dim = x_outargs->space()->dim();
//    f_sol_dim = solMgr_->getSolutionGroup()->getX().length();
    f_sol_dim = Teuchos::rcp_dynamic_cast< ::Thyra::LOCAAdaptiveState >(solMgr_->getState())
          ->getSolutionGroup()->getX().length();


  }

  if(Teuchos::is_null(x_outargs) || (x_args_dim != f_sol_dim)){ // g is not the right size

      x_final = Thyra::createMember(this->get_g_space(this->num_g()));

      mutable_outArgsPtr->set_g(this->num_g(), x_final);

  }
  else { // g is OK, use it
    x_final = x_outargs;
  }

  {
    // Deep copy final solution from LOCA group
    NOX::Thyra::Vector finalSolution(x_final);
//    finalSolution = solMgr_->getSolutionGroup()->getX();
    finalSolution = Teuchos::rcp_dynamic_cast< ::Thyra::LOCAAdaptiveState >(solMgr_->getState())
                      ->getSolutionGroup()->getX();

  }

  // If the arrays need resizing
  if(x_args_dim != f_sol_dim){

    const int parameterCount = this->Np();

    for (int pc = 0; pc < parameterCount; ++pc) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, this->num_g(), pc);
      const Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation dgdp_orient =
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM;
      if (dgdp_support.supports(dgdp_orient)) {
        const Thyra::ModelEvaluatorBase::DerivativeMultiVector<Scalar> dgdp =
          Thyra::create_DgDp_mv(*this, this->num_g(), pc, dgdp_orient);
        mutable_outArgsPtr->set_DgDp(this->num_g(), pc, dgdp);
      }
    }
  }

  // Compute responses for the final solution
  {
    Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs =
      this->getModel().createInArgs();
    {
      modelInArgs.set_x(x_final);
      modelInArgs.set_p(l, p_inargs);
    }

    Teuchos::ParameterList analysisParams;
    this->evalConvergedModelResponsesAndSensitivities(modelInArgs, outArgs, analysisParams);

    // Save the final solution TODO: this needs to be redone

    Teuchos::RCP<Thyra::ModelEvaluatorBase::InArgs<Scalar> > fp
         = Teuchos::rcp_const_cast<Thyra::ModelEvaluatorBase::InArgs<Scalar> >(finalPoint_);
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> ia;
    ia.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x, true);
    *fp = ia;
    fp->set_x(x_final);

  }
}


template <typename Scalar>
Teuchos::RCP<Piro::LOCAAdaptiveSolver<Scalar> >
Piro::observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Thyra::AdaptiveSolutionManager> &solMgr,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy =
    Teuchos::nonnull(observer) ?
    Teuchos::rcp(new Piro::ObserverToLOCASaveDataStrategyAdapter(observer)) :
    Teuchos::null;

  return Teuchos::rcp(new Piro::LOCAAdaptiveSolver<Scalar>(appParams, model, adjointModel, solMgr, saveDataStrategy));
}

template <typename Scalar>
Teuchos::RCP<LOCA::AdaptiveStepper>
Piro::LOCAAdaptiveSolver<Scalar>::getStepper()
{
  return stepper_;
}

#endif /* PIRO_LOCAADAPTIVESOLVER_DEF_HPP */
