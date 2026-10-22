// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_ROL_OBSERVEDOBJECTIVE_HPP
#define PIRO_THYRAPRODUCTME_ROL_OBSERVEDOBJECTIVE_HPP

#include <functional>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Thyra_VectorBase.hpp"

#include "ROL_Objective.hpp"
#include "ROL_UpdateType.hpp"
#include "ROL_Vector.hpp"

#include "Piro_ROL_ObserverBase.hpp"

namespace Piro {

/** \brief Decorator of the reduced objective of a transient optimization problem.
 *
 * It provides the observer calls and the objective recovery value that the
 * steady-state optimization provides through Piro::ThyraProductME_Objective_SimOpt
 * and Piro::ThyraProductME_Constraint_SimOpt:
 *  - observer->parametersChanged() is called when the parameters change;
 *  - observer->observeSolution() is called with the state at the final time,
 *    followed by observer->observeResponse(), for the initial guess and for
 *    the accepted iterates (every "Write Interval" iterations);
 *  - if "Objective Recovery Value" is set in the "Optimization Status" sublist
 *    and the forward solve did not converge ("State Solve Converged" = false in
 *    the "Optimization Status" sublist), value() returns the recovery value.
 * The "Optimizer Iteration Number" in the "Optimization Status" sublist is
 * updated at every update() call.
 *
 * The final state is provided by the finalState function, which returns the
 * state at the final time computed by the last forward solve if this was
 * performed with the given parameters, and null otherwise. In the latter case
 * the decorator recomputes the objective value at the given parameters.
 */
template <typename Real>
class ThyraProductME_ROL_ObservedObjective : public ROL::Objective<Real> {
public:
  using FinalStateFunction =
    std::function<Teuchos::RCP<const Thyra::VectorBase<Real>>(const ROL::Vector<Real>&)>;

  ThyraProductME_ROL_ObservedObjective(
    const Teuchos::RCP<ROL::Objective<Real>>& obj,
    const FinalStateFunction& finalState,
    Teuchos::ParameterList& piroParams,
    Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_HIGH,
    const Teuchos::RCP<ROL_ObserverBase<Real>>& observer = Teuchos::null) :
    obj_(obj),
    finalState_(finalState),
    optParams_(piroParams.sublist("Optimization Status")),
    out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
    verbosityLevel_(verbLevel),
    observer_(observer),
    updateType_(ROL::UpdateType::Temp),
    initialStateObserved_(false)
  {
    write_interval_ = optParams_.get("Write Interval", 1);
    optParams_.set<int>("Optimizer Iteration Number", -1);
    useObjectiveRecoveryValue_ = optParams_.isParameter("Objective Recovery Value");
    objectiveRecoveryValue_ = useObjectiveRecoveryValue_ ? optParams_.get<Real>("Objective Recovery Value") : Real(0.0);
  }

  virtual ~ThyraProductME_ROL_ObservedObjective() {}

  using ROL::Objective<Real>::update;

  void update(const ROL::Vector<Real>& x, ROL::UpdateType type, int iter = -1) override {
    if(verbosityLevel_ >= Teuchos::VERB_HIGH)
      *out_ << "Piro::ThyraProductME_ROL_ObservedObjective::update, UpdateType::" << ROL::UpdateTypeToString(type) << " Iter " << iter << std::endl;

    obj_->update(x, type, iter);
    updateType_ = type;
    optParams_.set<int>("Optimizer Iteration Number", iter);

    if (observer_.is_null()) return;

    if (type != ROL::UpdateType::Accept)
      observer_->parametersChanged();
    else if ((write_interval_ != 0) && (iter%write_interval_ == 0))
      observe(x, iter);
  }

  Real value(const ROL::Vector<Real>& x, Real& tol) override {
    Real val = obj_->value(x, tol);

    const bool converged = !optParams_.isParameter("State Solve Converged") ||
                           optParams_.get<bool>("State Solve Converged");

    if((updateType_ == ROL::UpdateType::Initial) && !initialStateObserved_ && observer_ != Teuchos::null) {
      observe(x, 0);
      initialStateObserved_ = true;
    }

    //set value to (large) recovery value if the forward solve did not converge
    if(useObjectiveRecoveryValue_ && !converged) {
      if(verbosityLevel_ >= Teuchos::VERB_LOW)
        *out_ << "Piro::ThyraProductME_ROL_ObservedObjective::value, Setting objective value to recovery value " << objectiveRecoveryValue_ << std::endl;
      val = objectiveRecoveryValue_;
    }
    return val;
  }

  void gradient(ROL::Vector<Real>& g, const ROL::Vector<Real>& x, Real& tol) override {
    obj_->gradient(g, x, tol);
  }

  void hessVec(ROL::Vector<Real>& hv, const ROL::Vector<Real>& v, const ROL::Vector<Real>& x, Real& tol) override {
    obj_->hessVec(hv, v, x, tol);
  }

  void invHessVec(ROL::Vector<Real>& hv, const ROL::Vector<Real>& v, const ROL::Vector<Real>& x, Real& tol) override {
    obj_->invHessVec(hv, v, x, tol);
  }

  void precond(ROL::Vector<Real>& Pv, const ROL::Vector<Real>& v, const ROL::Vector<Real>& x, Real& tol) override {
    obj_->precond(Pv, v, x, tol);
  }

private:

  void observe(const ROL::Vector<Real>& x, int iter) {
    Teuchos::RCP<const Thyra::VectorBase<Real>> state = finalState_(x);
    if (state.is_null()) {
      // The last forward solve was not performed at x, recompute it.
      Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
      obj_->value(x, tol);
      state = finalState_(x);
    }
    if (state.is_null()) {
      if(verbosityLevel_ >= Teuchos::VERB_LOW)
        *out_ << "WARNING: Piro::ThyraProductME_ROL_ObservedObjective: the state at the final time is not available, "
              << "skipping the observation of iteration " << iter << std::endl;
      return;
    }
    observer_->observeSolution(iter, *state, Teuchos::null, Teuchos::null, Teuchos::null);
    observer_->observeResponse(iter);
  }

  const Teuchos::RCP<ROL::Objective<Real>> obj_;
  const FinalStateFunction finalState_;
  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;
  ROL::UpdateType updateType_;
  bool initialStateObserved_;
  int write_interval_;
  bool useObjectiveRecoveryValue_;
  Real objectiveRecoveryValue_;
};

} // namespace Piro

#endif
