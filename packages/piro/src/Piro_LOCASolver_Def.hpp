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

#ifndef PIRO_LOCASOLVER_DEF_HPP
#define PIRO_LOCASOLVER_DEF_HPP

#include "Piro_LOCASolver.hpp"

#include "Piro_ObserverToLOCASaveDataStrategyAdapter.hpp"

#include "Thyra_DetachedVectorView.hpp"

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
Piro::LOCASolver<Scalar>::LOCASolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy) :
  SteadyStateSolver<Scalar>(model, model->Np() > 0), // Only one parameter supported
  piroParams_(piroParams),
  saveDataStrategy_(saveDataStrategy),
  globalData_(LOCA::createGlobalData(piroParams)),
  paramVector_(),
  group_(),
  locaStatusTests_(),
  noxStatusTests_(),
  stepper_()
{
  const int l = 0; // TODO: Allow user to select parameter index
  const Detail::ModelEvaluatorParamName paramName(this->getModel().get_p_names(l));
  const Thyra::Ordinal p_entry_count = this->getModel().get_p_space(l)->dim();
  for (Teuchos_Ordinal k = 0; k < p_entry_count; ++k) {
    (void) paramVector_.addParameter(paramName(k));
  }

  const NOX::Thyra::Vector initialGuess(*model->getNominalValues().get_x());
  group_ = Teuchos::rcp(new LOCA::Thyra::Group(globalData_, initialGuess, model, paramVector_, l));
  group_->setSaveDataStrategy(saveDataStrategy_);

  // TODO: Create non-trivial stopping criterion for the stepper
  locaStatusTests_ = Teuchos::null;

  // Create stopping criterion for the nonlinear solver
  const Teuchos::RCP<Teuchos::ParameterList> noxStatusParams =
    Teuchos::sublist(Teuchos::sublist(piroParams_, "NOX"), "Status Tests");
  noxStatusTests_ = NOX::StatusTest::buildStatusTests(*noxStatusParams, *(globalData_->locaUtils));

  stepper_ = Teuchos::rcp(new LOCA::Stepper(globalData_, group_, locaStatusTests_, noxStatusTests_, piroParams_));
  first_ = true;
}

template<typename Scalar>
Piro::LOCASolver<Scalar>::~LOCASolver()
{
  LOCA::destroyGlobalData(globalData_);
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
    std::cerr << "Continuation Stepper Finished.\n";
  } else if (status == LOCA::Abstract::Iterator::NotFinished) {
    std::cerr << "Continuation Stepper did not reach final value.\n";
  } else {
    std::cerr << "Nonlinear solver failed to converge.\n";
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

    this->evalConvergedModel(modelInArgs, outArgs);
  }
}


template <typename Scalar>
Teuchos::RCP<Piro::LOCASolver<Scalar> >
Piro::observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer)
{
  const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy =
    Teuchos::nonnull(observer) ?
    Teuchos::rcp(new Piro::ObserverToLOCASaveDataStrategyAdapter(observer)) :
    Teuchos::null;

  return Teuchos::rcp(new Piro::LOCASolver<Scalar>(appParams, model, saveDataStrategy));
}

#endif /* PIRO_LOCASOLVER_DEF_HPP */
