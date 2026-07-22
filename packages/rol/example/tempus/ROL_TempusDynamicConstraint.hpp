// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TEMPUS_DYNAMIC_CONSTRAINT_HPP
#define ROL_TEMPUS_DYNAMIC_CONSTRAINT_HPP


#include "Thyra_ModelEvaluator.hpp"

#include "Tempus_Integrator.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"

#include "ROL_DynamicConstraint.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"


namespace ROL{


template<class Real>
class TempusDynamicConstraint : public ROL::DynamicConstraint<Real> {
private:
  // Forward and (optional) adjoint model evaluator.
  ROL::Ptr<const Thyra::ModelEvaluator<Real>> model_;
  ROL::Ptr<const Thyra::ModelEvaluator<Real>> modelAdjoint_;
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

  // Flag indicating the use of the adjoint evaluator.
  bool usingAdjoint_;

public:

  TempusDynamicConstraint(const ROL::Ptr<Tempus::Integrator<Real>> & integrator);

  // The convention is that the adjoint model provides the application of
  // the adjoint Jacobian and its inverse.  All other operations are
  // provided by the forward model.
  TempusDynamicConstraint(const ROL::Ptr<Tempus::Integrator<Real>> & forward_integrator,
                          const ROL::Ptr<Tempus::Integrator<Real>> & adjoint_integrator);

  virtual ~TempusDynamicConstraint() {}

  void value(ROL::Vector<Real> &c,
             const ROL::Vector<Real> &uold,
             const ROL::Vector<Real> &unew,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts) const;

  void applyJacobian_uo(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const;

  void applyJacobian_un(ROL::Vector<Real> &jv,
                        const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &uold,
                        const ROL::Vector<Real> &unew,
                        const ROL::Vector<Real> &z,
                        const ROL::TimeStamp<Real> &ts) const;

  void applyJacobian_z(ROL::Vector<Real> &jv,
                       const ROL::Vector<Real> &v,
                       const ROL::Vector<Real> &uold,
                       const ROL::Vector<Real> &unew,
                       const ROL::Vector<Real> &z,
                       const ROL::TimeStamp<Real> &ts) const;

  void applyAdjointJacobian_uo(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const;

  void applyAdjointJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const;

  void applyAdjointJacobian_z(ROL::Vector<Real> &jv,
                              const ROL::Vector<Real> &v,
                              const ROL::Vector<Real> &uold,
                              const ROL::Vector<Real> &unew,
                              const ROL::Vector<Real> &z,
                              const ROL::TimeStamp<Real> &ts) const;

  void applyInverseJacobian_un(ROL::Vector<Real> &jv,
                               const ROL::Vector<Real> &v,
                               const ROL::Vector<Real> &uold,
                               const ROL::Vector<Real> &unew,
                               const ROL::Vector<Real> &z,
                               const ROL::TimeStamp<Real> &ts) const;

  void applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,
                                      const ROL::Vector<Real> &v,
                                      const ROL::Vector<Real> &uold,
                                      const ROL::Vector<Real> &unew,
                                      const ROL::Vector<Real> &z,
                                      const ROL::TimeStamp<Real> &ts) const;

};


//----------------------------------------------------------------------------
// Constructor with a single model.

template<class Real>
TempusDynamicConstraint<Real>::TempusDynamicConstraint(
  const ROL::Ptr<Tempus::Integrator<Real>> & integrator) {
  model_   = integrator->getStepper()->getModel();
  stepper_ = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(integrator->getStepper());
  Ju_      = model_->create_W();
  Ju_op_   = model_->create_W_op();
  Jz_op_   = model_->create_DfDp_op(0);
  // Set all other member variables to null.
  modelAdjoint_   = ROL::nullPtr;
  stepperAdjoint_ = ROL::nullPtr;
  adjointJu_      = ROL::nullPtr;
  adjointJu_op_   = ROL::nullPtr;
  // Set adjoint flag to false.
  usingAdjoint_ = false;
}


//----------------------------------------------------------------------------
// Constructor with forward and adjoint models.

template<class Real>
TempusDynamicConstraint<Real>::TempusDynamicConstraint(
  const ROL::Ptr<Tempus::Integrator<Real>> & forward_integrator,
  const ROL::Ptr<Tempus::Integrator<Real>> & adjoint_integrator) {
  model_          = forward_integrator->getStepper()->getModel();
  modelAdjoint_   = adjoint_integrator->getStepper()->getModel();
  stepper_        = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(forward_integrator->getStepper());
  stepperAdjoint_ = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(adjoint_integrator->getStepper());
  Ju_             = model_->create_W();
  Ju_op_          = model_->create_W_op();
  Jz_op_          = model_->create_DfDp_op(0);
  adjointJu_      = modelAdjoint_->create_W();
  adjointJu_op_   = modelAdjoint_->create_W_op();
  // Set adjoint flag to true.
  usingAdjoint_ = true;
}


//----------------------------------------------------------------------------
// Value

template<class Real>
void TempusDynamicConstraint<Real>::value(ROL::Vector<Real> &c,
                                          const ROL::Vector<Real> &uold,
                                          const ROL::Vector<Real> &unew,
                                          const ROL::Vector<Real> &z,
                                          const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&          rtv_c = dynamic_cast<ROL::ThyraVector<Real>& >(c);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  stepper_->computeStepResidual(*(rtv_c.getVector()), x, t, *(rtv_z.getVector()), 0);
} // value


//----------------------------------------------------------------------------
// Partial Jacobians

template<class Real>
void TempusDynamicConstraint<Real>::applyJacobian_uo(ROL::Vector<Real> &jv,
                                                     const ROL::Vector<Real> &v,
                                                     const ROL::Vector<Real> &uold,
                                                     const ROL::Vector<Real> &unew,
                                                     const ROL::Vector<Real> &z,
                                                     const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  int deriv_index = 1; // 1 = old state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_uo

template<class Real>
void TempusDynamicConstraint<Real>::applyJacobian_un(ROL::Vector<Real> &jv,
                                                     const ROL::Vector<Real> &v,
                                                     const ROL::Vector<Real> &uold,
                                                     const ROL::Vector<Real> &unew,
                                                     const ROL::Vector<Real> &z,
                                                     const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  int deriv_index = 0; // 0 = new state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_un

template<class Real>
void TempusDynamicConstraint<Real>::applyJacobian_z(ROL::Vector<Real> &jv,
                                                    const ROL::Vector<Real> &v,
                                                    const ROL::Vector<Real> &uold,
                                                    const ROL::Vector<Real> &unew,
                                                    const ROL::Vector<Real> &z,
                                                    const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  stepper_->computeStepParamDeriv(*Jz_op_, x, t, *(rtv_z.getVector()), 0);
  Jz_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyJacobian_z


//----------------------------------------------------------------------------
// Adjoint partial Jacobians

template<class Real>
void TempusDynamicConstraint<Real>::applyAdjointJacobian_uo(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  if (!usingAdjoint_) {
    int deriv_index = 1; // 1 = old state
    stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
  else {
    int deriv_index = 1; // 1 = old state
    stepperAdjoint_->computeStepJacobian(*adjointJu_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
    //adjointJu_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
} // applyAdjointJacobian_uo

template<class Real>
void TempusDynamicConstraint<Real>::applyAdjointJacobian_un(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  if (!usingAdjoint_) {
    int deriv_index = 0; // 0 = new state
    stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
  else {
    int deriv_index = 0; // 0 = new state
    stepperAdjoint_->computeStepJacobian(*adjointJu_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
    adjointJu_op_->apply(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
  }
} // applyAdjointJacobian_un

template<class Real>
void TempusDynamicConstraint<Real>::applyAdjointJacobian_z(ROL::Vector<Real> &jv,
                                                           const ROL::Vector<Real> &v,
                                                           const ROL::Vector<Real> &uold,
                                                           const ROL::Vector<Real> &unew,
                                                           const ROL::Vector<Real> &z,
                                                           const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  stepper_->computeStepParamDeriv(*Jz_op_, x, t, *(rtv_z.getVector()), 0);
  Jz_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
} // applyAdjointJacobian_z


//----------------------------------------------------------------------------
// Inverses

template<class Real>
void TempusDynamicConstraint<Real>::applyInverseJacobian_un(ROL::Vector<Real> &jv,
                                                            const ROL::Vector<Real> &v,
                                                            const ROL::Vector<Real> &uold,
                                                            const ROL::Vector<Real> &unew,
                                                            const ROL::Vector<Real> &z,
                                                            const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria = Teuchos::null;
  stepper_->computeStepSolver(*Ju_, x, t, *(rtv_z.getVector()), 0);
  Ju_->solve(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
} // applyInverseJacobian_un

template<class Real>
void TempusDynamicConstraint<Real>::applyInverseAdjointJacobian_un(ROL::Vector<Real> &jv,
                                                                   const ROL::Vector<Real> &v,
                                                                   const ROL::Vector<Real> &uold,
                                                                   const ROL::Vector<Real> &unew,
                                                                   const ROL::Vector<Real> &z,
                                                                   const ROL::TimeStamp<Real> &ts) const {
  ROL::ThyraVector<Real>&         rtv_jv = dynamic_cast<ROL::ThyraVector<Real>& >(jv);
  const ROL::ThyraVector<Real>&    rtv_v = dynamic_cast<const ROL::ThyraVector<Real>& >(v);
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

  const Teuchos::Ptr<const Thyra::SolveCriteria<Real>> solveCriteria = Teuchos::null;
  if (!usingAdjoint_) {
    stepper_->computeStepSolver(*Ju_, x, t, *(rtv_z.getVector()), 0);
    Ju_->solve(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
  }
  else {
    stepperAdjoint_->computeStepSolver(*adjointJu_, x, t, *(rtv_z.getVector()), 0);
    adjointJu_->solve(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
  }
} // applyInverseAdjointJacobian_un


} // namespace ROL


#endif
