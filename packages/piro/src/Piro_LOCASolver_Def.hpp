// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_LOCASOLVER_DEF_HPP
#define PIRO_LOCASOLVER_DEF_HPP

#include "Piro_LOCASolver.hpp"

#include "Piro_ObserverToLOCASaveDataStrategyAdapter.hpp"

#include "Thyra_DetachedVectorView.hpp"

#include "NOX_StatusTest_Factory.H"

#include "Piro_MatrixFreeDecorator.hpp" 

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
Piro::LOCASolver<Scalar>::LOCASolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy) :
  SteadyStateSolver<Scalar>(model, adjointModel, model->Np() > 0), // Only one parameter supported
  piroParams_(piroParams),
  saveDataStrategy_(saveDataStrategy),
  globalData_(LOCA::createGlobalData(piroParams)),
  paramVector_(),
  group_(),
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

  const NOX::Thyra::Vector initialGuess(*model_->getNominalValues().get_x());

  // Don't set a scaling vector here; it'll make calling LOCA impossible with a
  // Jacobian that is not a Tpetra::RowMatrix (but only a Tpetra::Operator).
  Teuchos::RCP<Thyra::VectorBase<double> > scaling_vector_ = Teuchos::null;

  group_ = Teuchos::rcp(new LOCA::Thyra::Group(globalData_, initialGuess, model_, paramVector_, l, false, scaling_vector_));
  group_->setSaveDataStrategy(saveDataStrategy_);

  // TODO: Create non-trivial stopping criterion for the stepper
  locaStatusTests_ = Teuchos::null;

  // Create stopping criterion for the nonlinear solver
  const Teuchos::RCP<Teuchos::ParameterList> noxStatusParams =
    Teuchos::sublist(Teuchos::sublist(piroParams_, "NOX"), "Status Tests");
  noxStatusTests_ = NOX::StatusTest::buildStatusTests(*noxStatusParams, *(globalData_->locaUtils));
  

  stepper_ = Teuchos::rcp(new LOCA::Stepper(globalData_, group_, locaStatusTests_, noxStatusTests_, piroParams_));
  first_ = true;

  if (piroParams_->isSublist("NOX") &&
      piroParams_->sublist("NOX").isSublist("Printing"))
    utils_.reset(piroParams_->sublist("NOX").sublist("Printing"));

  this->setSensitivityMethod("Forward");
}

template<typename Scalar>
Piro::LOCASolver<Scalar>::~LOCASolver()
{
  LOCA::destroyGlobalData(globalData_);
}

template<typename Scalar>
Teuchos::RCP<NOX::Solver::Generic>
Piro::LOCASolver<Scalar>::getSolver()
{
  return stepper_->getSolver();
}

template<typename Scalar>
Teuchos::ParameterList &
Piro::LOCASolver<Scalar>::getStepperParams()
{
  return stepper_->getParams();
}

template<typename Scalar>
Teuchos::ParameterList &
Piro::LOCASolver<Scalar>::getStepSizeParams()
{
  return stepper_->getStepSizeParams();
}

template <typename Scalar>
Teuchos::RCP<LOCA::Stepper>
Piro::LOCASolver<Scalar>::getStepper()
{
  return stepper_;
}

template <typename Scalar>
void
Piro::LOCASolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  const int l = 0; // TODO: Allow user to select parameter index
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > p_inargs = inArgs.get_p(l);

  // Forward parameter values to the LOCA stepper
  {
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> > p_inargs_or_nominal =
      Teuchos::nonnull(p_inargs) ? p_inargs : this->getNominalValues().get_p(l);
    const Thyra::ConstDetachedVectorView<Scalar> p_init_values(p_inargs_or_nominal);
    const Teuchos_Ordinal p_entry_count = p_init_values.subDim();
    TEUCHOS_ASSERT(p_entry_count == Teuchos::as<Teuchos_Ordinal>(paramVector_.length()));

    for (Teuchos_Ordinal k = 0; k < p_entry_count; ++k) {
      paramVector_[k] = p_init_values[k];
    }

    group_->setParams(paramVector_);
  }

  if (first_) {
    // No need to call reset. The call can result in long stdout output, so it's
    // nice to avoid it since we can.
    first_ = false;
  } else
    stepper_->reset(globalData_, group_, locaStatusTests_, noxStatusTests_, piroParams_);
  const LOCA::Abstract::Iterator::IteratorStatus status = stepper_->run();

  if (status == LOCA::Abstract::Iterator::Finished) {
    utils_.out() << "Continuation Stepper Finished.\n";;
  } else if (status == LOCA::Abstract::Iterator::NotFinished) {
    utils_.out() << "Continuation Stepper did not reach final value.\n";
  } else {
    utils_.out() << "Nonlinear solver failed to converge.\n";
    outArgs.setFailed();
  }

  const Teuchos::RCP<Thyra::VectorBase<Scalar> > x_outargs = outArgs.get_g(this->num_g());
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > x_final =
    Teuchos::nonnull(x_outargs) ? x_outargs : Thyra::createMember(this->get_g_space(this->num_g()));

  {
    // Deep copy final solution from LOCA group
    NOX::Thyra::Vector finalSolution(x_final);
    finalSolution = group_->getX();
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
  }
}


template <typename Scalar>
Teuchos::RCP<Piro::LOCASolver<Scalar> >
Piro::observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &adjointModel,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy =
    Teuchos::nonnull(observer) ?
    Teuchos::rcp(new Piro::ObserverToLOCASaveDataStrategyAdapter(observer)) :
    Teuchos::null;

  return Teuchos::rcp(new Piro::LOCASolver<Scalar>(appParams, model, adjointModel, saveDataStrategy));
}

#endif /* PIRO_LOCASOLVER_DEF_HPP */
