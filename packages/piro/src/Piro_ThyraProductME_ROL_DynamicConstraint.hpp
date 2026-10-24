// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_THYRAPRODUCTME_ROL_DYNAMICCONSTRAINT_HPP
#define PIRO_THYRAPRODUCTME_ROL_DYNAMICCONSTRAINT_HPP

#include <algorithm>
#include <cmath>

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_Integrator.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"

#include "ROL_DynamicConstraint.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"
#include "Piro_TempusIntegrator.hpp"
#include "Piro_ROL_ObserverBase.hpp"


namespace Piro{

/** \brief ROL dynamic constraint given by one time step of a Tempus stepper.
 *
 * The residual and its derivatives are computed through the
 * Tempus::StepperOptimizationInterface of the stepper (currently only
 * implemented by the Backward Euler stepper).
 *
 * The convention is that the adjoint model, if provided, provides the
 * application of the adjoint Jacobian and its inverse; all other operations
 * are provided by the forward model. Without an adjoint model, the adjoint
 * Jacobian is applied as the transpose of the Jacobian, and its inverse is
 * applied by solving with the transposed Jacobian, using the linear solver
 * (and preconditioner) factory of the forward model.
 *
 * The time steps are solved with the Newton solver of ROL::DynamicConstraint,
 * configured by the "Dynamic Constraint"->"Solve" sublist passed to
 * setSolveParameters(). After each time step the "State Solve Converged" flag
 * of the "Optimization Status" sublist is set to whether all the time steps
 * solved so far in the current trajectory converged; time stamp 1 is the first
 * time step of a trajectory. The state at the time stamp set with
 * setFinalTimeStampIndex() is stored and can be retrieved with getFinalState().
 */
template<class Real>
class ThyraProductME_ROL_DynamicConstraint : public ROL::DynamicConstraint<Real> {
private:
  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;

  // Forward and (optional) adjoint model evaluator.
  ROL::Ptr<const Thyra::ModelEvaluator<Real>> model_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Real>> modelAdjoint_;
  // Forward and (optional) adjoint steppers.
  ROL::Ptr<Tempus::StepperOptimizationInterface<Real>> stepper_;
  ROL::Ptr<Tempus::StepperOptimizationInterface<Real>> stepperAdjoint_;
  // Forward linear operators.
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> Ju_;
  ROL::Ptr<Thyra::LinearOpBase<Real>> Ju_op_;
  ROL::Ptr<Thyra::LinearOpBase<Real>> Jz_op_;
  // (Optional) Adjoint linear operators.
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> adjointJu_;
  ROL::Ptr<Thyra::LinearOpBase<Real>> adjointJu_op_;
  // Solver and preconditioner of the transposed Jacobian (used without adjoint model).
  mutable ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> transposedJu_;
  mutable ROL::Ptr<Thyra::PreconditionerBase<Real>> transposedJuPrec_;

  // Flag indicating the use of the adjoint evaluator.
  bool usingAdjoint_;

  // Copies of the time step solver tolerances, used to check convergence.
  Real solveAtol_;
  Real solveRtol_;
  bool solveZeroInitialGuess_;

  // Convergence of the current trajectory and state at the final time stamp.
  bool trajectoryConverged_;
  int finalTimeStampIndex_;
  ROL::Ptr<ROL::Vector<Real>> finalState_;
  ROL::Ptr<ROL::Vector<Real>> finalControl_;

  using SMT=Teuchos::ScalarTraits<typename Thyra::ModelEvaluator<Real>::ScalarMag>;

  void initialize(const Teuchos::RCP<Tempus::Integrator<Real>> & forward_integrator,
                  const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
                  const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> & modelAdjoint);

  ROL::Ptr<Tempus::StepperOptimizationInterface<Real>>
  getOptimizationStepper(const Teuchos::RCP<Tempus::Integrator<Real>> & integrator) const;

  void solveTransposedJacobian(const Thyra::VectorBase<Real>& v, Thyra::VectorBase<Real>& jv,
                               const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>>& x,
                               const Teuchos::Array<Real>& t,
                               const Thyra::VectorBase<Real>& z) const;

  // States and times in the format of Tempus::StepperOptimizationInterface:
  // x[0], t[0] are the NEW state and time, x[1], t[1] the OLD ones.
  static void dynamicConstraintStencil(const ROL::Vector<Real> &uold, const ROL::Vector<Real> &unew,
                                       const ROL::TimeStamp<Real> &ts,
                                       Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>>& x,
                                       Teuchos::Array<Real>& t) {
    x.resize(2); t.resize(2);
    x[0] = dynamic_cast<const ROL::ThyraVector<Real>&>(unew).getVector(); t[0] = ts.t[1];
    x[1] = dynamic_cast<const ROL::ThyraVector<Real>&>(uold).getVector(); t[1] = ts.t[0];
  }

  static void checkDirection(const ROL::Vector<Real> &v, const std::string& function) {
    TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(dynamic_cast<const ROL::ThyraVector<Real>&>(v).getVector()->norm_2()), std::logic_error,
      "ThyraProductME_ROL_DynamicConstraint<Real>::" << function << ", the norm of rtv_v is not a number." << std::endl);
  }

public:

  //! Constructor without adjoint model: the adjoint Jacobian is the transpose of the Jacobian.
  ThyraProductME_ROL_DynamicConstraint(const Teuchos::RCP<Tempus::Integrator<Real>> & integrator,
    Teuchos::ParameterList& piroParams,
    Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  //! Constructor with forward and adjoint models. If the adjoint integrator or the
  //! adjoint model is null, this is equivalent to the constructor without adjoint model.
  ThyraProductME_ROL_DynamicConstraint(const Teuchos::RCP<Tempus::Integrator<Real>> & forward_integrator,
    const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoint,
    Teuchos::ParameterList& piroParams,
    Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  virtual ~ThyraProductME_ROL_DynamicConstraint() {}

  //! Index of the time stamp whose state is stored (see getFinalState()).
  void setFinalTimeStampIndex(int index) {
    finalTimeStampIndex_ = index;
  }

  /** \brief State at the final time stamp computed by the last solve of that
   *  time step, or null if it was not computed with control z. */
  Teuchos::RCP<const Thyra::VectorBase<Real>> getFinalState(const ROL::Vector<Real>& z) const;

  void setSolveParameters(Teuchos::ParameterList& parlist) override;

  void solve(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &uold,
             ROL::Vector<Real> &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts) override;

  void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &uold,
             const ROL::Vector<Real> &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts) const override;

  void applyJacobian_uo(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const override;

  void applyJacobian_un(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const override;

  void applyJacobian_z(ROL::Vector<Real> &jv,
                       const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &uold,
                       const ROL::Vector<Real> &unew,
                       const ROL::Vector<Real> &z,
                       const ROL::TimeStamp<Real> &ts) const override;

  void applyAdjointJacobian_uo(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const override;

  void applyAdjointJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const override;

  void applyAdjointJacobian_z(ROL::Vector<Real> &jv,
                              const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &uold,
                              const ROL::Vector<Real> &unew,
                              const ROL::Vector<Real> &z,
                              const ROL::TimeStamp<Real> &ts) const override;

  void applyInverseJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const override;

  void applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,
                                      const ROL::Vector<Real> &v,
                                      const ROL::Vector<Real> &uold,
                                      const ROL::Vector<Real> &unew,
                                      const ROL::Vector<Real> &z,
                                      const ROL::TimeStamp<Real> &ts) const override;

};


//----------------------------------------------------------------------------
// Constructors

template<class Real>
ThyraProductME_ROL_DynamicConstraint<Real>::ThyraProductME_ROL_DynamicConstraint(
  const Teuchos::RCP<Tempus::Integrator<Real>>& integrator,
  Teuchos::ParameterList& piroParams,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> /* observer */) :
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel) {
  initialize(integrator, Teuchos::null, Teuchos::null);
}

template<class Real>
ThyraProductME_ROL_DynamicConstraint<Real>::ThyraProductME_ROL_DynamicConstraint(
  const Teuchos::RCP<Tempus::Integrator<Real>>& forward_integrator,
  const Teuchos::RCP<Tempus::Integrator<Real>>& adjoint_integrator,
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>>& modelAdjoint,
  Teuchos::ParameterList& piroParams,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> /* observer */) :
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel) {
  initialize(forward_integrator, adjoint_integrator, modelAdjoint);
}

template<class Real>
ROL::Ptr<Tempus::StepperOptimizationInterface<Real>>
ThyraProductME_ROL_DynamicConstraint<Real>::getOptimizationStepper(
  const Teuchos::RCP<Tempus::Integrator<Real>> & integrator) const {
  auto stepper = integrator->getStepper();
  TEUCHOS_TEST_FOR_EXCEPTION(stepper.is_null(), std::logic_error,
    std::endl << "Piro::ThyraProductME_ROL_DynamicConstraint: the Tempus integrator has no stepper." << std::endl);
  auto optStepper = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(stepper);
  TEUCHOS_TEST_FOR_EXCEPTION(optStepper.is_null(), Teuchos::Exceptions::InvalidParameter,
    std::endl << "Piro::ThyraProductME_ROL_DynamicConstraint: the Tempus stepper \"" << stepper->getStepperType() <<
    "\" does not implement Tempus::StepperOptimizationInterface, which is currently only implemented by the " <<
    "\"Backward Euler\" stepper. Use a Backward Euler stepper, or set the ROL parameter \"Tempus Driver\" to true." << std::endl);
  return optStepper;
}

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::initialize(
  const Teuchos::RCP<Tempus::Integrator<Real>> & forward_integrator,
  const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
  const Teuchos::RCP<const Thyra::ModelEvaluator<Real>> & modelAdjoint) {
  model_   = forward_integrator->getStepper()->getModel();
  stepper_ = getOptimizationStepper(forward_integrator);
  Ju_      = model_->create_W();
  Ju_op_   = model_->create_W_op();
  Jz_op_   = model_->create_DfDp_op(0);

  usingAdjoint_ = Teuchos::nonnull(adjoint_integrator) && Teuchos::nonnull(modelAdjoint);
  if (usingAdjoint_) {
    modelAdjoint_   = modelAdjoint;
    stepperAdjoint_ = getOptimizationStepper(adjoint_integrator);
    adjointJu_      = modelAdjoint_->create_W();
    adjointJu_op_   = modelAdjoint_->create_W_op();
  }

  // Same defaults as ROL::DynamicConstraint
  solveAtol_ = 1e-4*std::sqrt(ROL::ROL_EPSILON<Real>());
  solveRtol_ = 1.0;
  solveZeroInitialGuess_ = false;

  trajectoryConverged_ = true;
  finalTimeStampIndex_ = -1;
}

//----------------------------------------------------------------------------
// Time step solve

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::setSolveParameters(Teuchos::ParameterList& parlist) {
  ROL::DynamicConstraint<Real>::setSolveParameters(parlist);
  Teuchos::ParameterList& list = parlist.sublist("Dynamic Constraint").sublist("Solve");
  solveAtol_ = list.get("Absolute Residual Tolerance", solveAtol_);
  solveRtol_ = list.get("Relative Residual Tolerance", solveRtol_);
  solveZeroInitialGuess_ = list.get("Zero Initial Guess", solveZeroInitialGuess_);
}

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::solve(ROL::Vector<Real> &c,
                                                       const ROL::Vector<Real> &uold,
                                                       ROL::Vector<Real> &unew,
                                                       const ROL::Vector<Real> &z,
                                                       const ROL::TimeStamp<Real> &ts) {
  // Initial residual, with the same initial guess used by ROL::DynamicConstraint::solve
  if (solveZeroInitialGuess_) unew.zero();
  this->update(uold, unew, z, ts);
  value(c, uold, unew, z, ts);
  const Real cnorm0 = c.norm();

  ROL::DynamicConstraint<Real>::solve(c, uold, unew, z, ts);

  // At the end of the solve c is the residual at unew. Same convergence criterion as ROL.
  const Real cnorm = c.norm();
  const Real ctol = std::min(solveAtol_, solveRtol_*cnorm0);
  const bool converged = (cnorm <= ctol);

  // Time stamp 1 is the first time step of a trajectory in ROL::ReducedDynamicObjective
  const bool firstFailure = !converged && ((ts.k <= 1) || trajectoryConverged_);
  trajectoryConverged_ = (ts.k <= 1) ? converged : (trajectoryConverged_ && converged);
  optParams_.set("State Solve Converged", trajectoryConverged_);

  // Report only the first time step of the trajectory that did not converge
  if (firstFailure && verbosityLevel_ >= Teuchos::VERB_LOW)
    *out_ << "WARNING: Piro::ThyraProductME_ROL_DynamicConstraint::solve: the solve of the time step [" << ts.t.front()
          << ", " << ts.t.back() << "] did not converge, residual norm " << cnorm << " > tolerance " << ctol << std::endl;

  if (static_cast<int>(ts.k) == finalTimeStampIndex_) {
    if (finalState_.is_null()) finalState_ = unew.clone();
    if (finalControl_.is_null()) finalControl_ = z.clone();
    finalState_->set(unew);
    finalControl_->set(z);
  }
}

template<class Real>
Teuchos::RCP<const Thyra::VectorBase<Real>>
ThyraProductME_ROL_DynamicConstraint<Real>::getFinalState(const ROL::Vector<Real>& z) const {
  if (finalState_.is_null() || finalControl_.is_null())
    return Teuchos::null;
  auto diff = finalControl_->clone();
  diff->set(*finalControl_);
  diff->axpy(Real(-1.0), z);
  if (diff->norm() != Real(0.0))
    return Teuchos::null;
  return dynamic_cast<const ROL::ThyraVector<Real>&>(*finalState_).getVector();
}

//----------------------------------------------------------------------------
// Value

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::value(ROL::Vector<Real> &c,
                                          const ROL::Vector<Real> &uold,
                                          const ROL::Vector<Real> &unew,
                                          const ROL::Vector<Real> &z,
                                          const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&          rtv_c = dynamic_cast<ROL::ThyraVector<Real>& >(c);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  stepper_->computeStepResidual(*(rtv_c.getVector()), x, t, *(rtv_z.getVector()), 0);
} // value


//----------------------------------------------------------------------------
// Partial Jacobians

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_uo(ROL::Vector<Real> &jv,
                                                     const ROL::Vector<Real> &v,
                                                     const ROL::Vector<Real> &uold,
                                                     const ROL::Vector<Real> &unew,
                                                     const ROL::Vector<Real> &z,
                                                     const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyJacobian_uo");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  int deriv_index = 1; // 1 = old state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_uo

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_un(ROL::Vector<Real> &jv,
                                                     const ROL::Vector<Real> &v,
                                                     const ROL::Vector<Real> &uold,
                                                     const ROL::Vector<Real> &unew,
                                                     const ROL::Vector<Real> &z,
                                                     const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyJacobian_un");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  int deriv_index = 0; // 0 = new state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_un

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_z(ROL::Vector<Real> &jv,
                                                    const ROL::Vector<Real> &v,
                                                    const ROL::Vector<Real> &uold,
                                                    const ROL::Vector<Real> &unew,
                                                    const ROL::Vector<Real> &z,
                                                    const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyJacobian_z");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  stepper_->computeStepParamDeriv(*Jz_op_, x, t, *(rtv_z.getVector()), 0);

  rtv_jv.getVector()->assign(0.0);
  Jz_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_z


//----------------------------------------------------------------------------
// Adjoint partial Jacobians

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_uo(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyAdjointJacobian_uo");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  int deriv_index = 1; // 1 = old state
  if (!usingAdjoint_) {
    stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
  else {
    stepperAdjoint_->computeStepJacobian(*adjointJu_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    rtv_jv.getVector()->assign(0.0);
    adjointJu_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
} // applyAdjointJacobian_uo

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_un(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyAdjointJacobian_un");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  int deriv_index = 0; // 0 = new state
  if (!usingAdjoint_) {
    stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
  else {
    stepperAdjoint_->computeStepJacobian(*adjointJu_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    rtv_jv.getVector()->assign(0.0);
    adjointJu_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
} // applyAdjointJacobian_un

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_z(ROL::Vector<Real> &jv,
                                                           const ROL::Vector<Real> &v,
                                                           const ROL::Vector<Real> &uold,
                                                           const ROL::Vector<Real> &unew,
                                                           const ROL::Vector<Real> &z,
                                                           const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyAdjointJacobian_z");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  stepper_->computeStepParamDeriv(*Jz_op_, x, t, *(rtv_z.getVector()), 0);
  Jz_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyAdjointJacobian_z


//----------------------------------------------------------------------------
// Inverses

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyInverseJacobian_un(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyInverseJacobian_un");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria = Teuchos::null;
  stepper_->computeStepSolver(*Ju_, x, t, *(rtv_z.getVector()), 0);
  rtv_jv.getVector()->assign(0.0);
  Ju_->solve(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
} // applyInverseJacobian_un

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,
                                                                   const ROL::Vector<Real> &v,
                                                                   const ROL::Vector<Real> &uold,
                                                                   const ROL::Vector<Real> &unew,
                                                                   const ROL::Vector<Real> &z,
                                                                   const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);
  checkDirection(v, "applyInverseAdjointJacobian_un");

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x;
  Teuchos::Array<Real> t;
  dynamicConstraintStencil(uold, unew, ts, x, t);

  const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria = Teuchos::null;
  if (!usingAdjoint_) {
    solveTransposedJacobian(*(rtv_v.getVector()), *(rtv_jv.getVector()), x, t, *(rtv_z.getVector()));
  }
  else {
    stepperAdjoint_->computeStepSolver(*adjointJu_, x, t, *(rtv_z.getVector()), 0);
    rtv_jv.getVector()->assign(0.0);
    adjointJu_->solve(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
  }
} // applyInverseAdjointJacobian_un

template<class Real>
void ThyraProductME_ROL_DynamicConstraint<Real>::solveTransposedJacobian(
  const Thyra::VectorBase<Real>& v, Thyra::VectorBase<Real>& jv,
  const Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>>& x,
  const Teuchos::Array<Real>& t,
  const Thyra::VectorBase<Real>& z) const {

  const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria = Teuchos::null;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Real>> lows_factory = model_->get_W_factory();

  if (lows_factory.is_null()) {
    // No linear solver factory: use a transposed solve with the Jacobian solver of the model.
    stepper_->computeStepSolver(*Ju_, x, t, z, 0);
    TEUCHOS_TEST_FOR_EXCEPTION(!Ju_->solveSupports(Thyra::TRANS), std::logic_error,
      std::endl << "Piro::ThyraProductME_ROL_DynamicConstraint: without an adjoint model evaluator, the adjoint " <<
      "Jacobian is inverted with a transposed solve, which is not supported by the Jacobian solver of the model. " <<
      "Provide an adjoint model evaluator, or a model with a linear solver factory (get_W_factory())." << std::endl);
    jv.assign(0.0);
    Ju_->solve(Thyra::TRANS, v, Teuchos::ptrFromRef(jv), solveCriteria);
    return;
  }

  // Build a solver for the transposed Jacobian (and transposed preconditioner),
  // as done by Piro::ThyraProductME_Constraint_SimOpt for steady problems.
  int deriv_index = 0; // 0 = new state
  stepper_->computeStepJacobian(*Ju_op_, x, t, z, 0, deriv_index);
  Teuchos::RCP<const Thyra::LinearOpBase<Real>> Ju_op = Ju_op_;

  if (transposedJu_.is_null())
    transposedJu_ = lows_factory->createOp();

  Teuchos::RCP<Thyra::PreconditionerFactoryBase<Real>> prec_factory = lows_factory->getPreconditionerFactory();
  if (Teuchos::nonnull(prec_factory)) {
    if (transposedJuPrec_.is_null())
      transposedJuPrec_ = prec_factory->createPrec();
    prec_factory->initializePrec(Thyra::defaultLinearOpSource<Real>(Ju_op), transposedJuPrec_.get());
    Thyra::initializePreconditionedOp<Real>(*lows_factory,
        Thyra::transpose<Real>(Ju_op),
        Thyra::unspecifiedPrec<Real>(Thyra::transpose<Real>(transposedJuPrec_->getUnspecifiedPrecOp())),
        transposedJu_.ptr());
  }
  else {
    Thyra::initializeOp<Real>(*lows_factory, Thyra::transpose<Real>(Ju_op), transposedJu_.ptr());
  }

  jv.assign(0.0);
  transposedJu_->solve(Thyra::NOTRANS, v, Teuchos::ptrFromRef(jv), solveCriteria);
}

} // namespace Piro


#endif
