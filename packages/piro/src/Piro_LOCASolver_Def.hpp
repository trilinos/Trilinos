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

  stepper_->reset(globalData_, group_, locaStatusTests_, noxStatusTests_, piroParams_);
  stepper_->setCallerWillCallPrintSolution(true);

  LOCA::Abstract::Iterator::IteratorStatus istat;
  LOCA::Abstract::Iterator::StepStatus sstat;
  Thyra::ModelEvaluatorBase::InArgs<Scalar>
    modelInArgs = this->getModel().createInArgs();
  stepper_->initializeIterateOnce(istat);
  for (;;) {
    std::cout.precision(18);
    const double cpb = stepper_->getContinuationParameter();
    const bool keep_iterating = stepper_->iterateOnce(istat, sstat);
    const double cpa = stepper_->getContinuationParameter();
    // Call evalConvergedModel()?
    bool
      // Call evalConvergedModel() if the step was successful ...
      call_ecm = sstat == LOCA::Abstract::Iterator::Successful
      // and if work was done. LOCA::Stepper::finish() may or may not do an
      // extra step. (finish() was called if !keep_iterating.) Distinguish
      // between these two cases by comparing the continuation parameter before
      // and after. If they are not the same, then finish() did an extra step,
      // and we call evalConvergedModel() on the result; if they are the same,
      // it was already evaluated.
      && (keep_iterating || cpa != cpb);

    if (!keep_iterating) {
      // Just after the call to finish().
      if (istat == LOCA::Abstract::Iterator::Finished) {
        std::cout << "Continuation Stepper Finished.\n";
      } else if (istat == LOCA::Abstract::Iterator::NotFinished) {
        std::cout << "Continuation Stepper did not reach final value.\n";
      } else {
        std::cout << "Nonlinear solver failed to converge.\n";
        outArgs.setFailed();
        // The behavior of the previous implementation of evalModelImpl was to
        // call evalConvergedModel regardless of the iterator status. Do that in
        // the new implementation, too. This has the effect of calling
        // evalConvergedModel() twice on the same point. But a subtlety is that
        // outArgs.setFailed() sets get_g(num_g) to null. Only the call to
        // get_g_space in the following code restores it. If it is not restored,
        // then exit tests based on response functions will show nan or
        // equivalent. In any case, the redundant call occurs only in the case
        // of failure.
        call_ecm = true;
      }
    }

    if (call_ecm) {
      // Compute responses at the end of this continuation step.
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >
        x_outargs = outArgs.get_g(this->num_g());
      const Teuchos::RCP<Thyra::VectorBase<Scalar> >
        x_current = (Teuchos::nonnull(x_outargs) ?
                     x_outargs :
                     Thyra::createMember(this->get_g_space(this->num_g())));

      { // Deep copy solution for this continuation step from LOCA
        // group. Vector::operator= does a deep copy into the memory that
        // x_current points to. v_current itself does not need to be kept; it's
        // used just for its operator=.
        NOX::Thyra::Vector v_current(x_current);
        v_current = group_->getX(); }

      modelInArgs.set_x(x_current);
      modelInArgs.set_p(l, p_inargs);

      this->evalConvergedModel(modelInArgs, outArgs);
      stepper_->callPrintSolution();
    }

    if (!keep_iterating) break;
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
