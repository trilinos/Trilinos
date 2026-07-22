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


#include "Thyra_ModelEvaluator.hpp"

#include "Tempus_Integrator.hpp"
#include "Tempus_StepperOptimizationInterface.hpp"

#include "ROL_DynamicConstraint.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_ThyraVector.hpp"
#include "Piro_TempusIntegrator.hpp" 


namespace Piro{


template<class Real>
class ThyraProductME_ROL_DynamicConstraint : public ROL::DynamicConstraint<Real> {
private:
  Teuchos::RCP<Tempus::Integrator<Real>> integrator_;
  int num_responses_;
  Teuchos::ParameterList& optParams_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosityLevel_;
  Teuchos::RCP<ROL_ObserverBase<Real>> observer_;

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

  // Flag indicating the use of the adjoint evaluator.
  bool usingAdjoint_;

  using SMT=Teuchos::ScalarTraits<typename Thyra::ModelEvaluator<Real>::ScalarMag>;

public:

  ThyraProductME_ROL_DynamicConstraint(const Teuchos::RCP<Tempus::Integrator<Real>> & integrator,
    Teuchos::ParameterList& piroParams,
    Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  // The convention is that the adjoint model provides the application of
  // the adjoint Jacobian and its inverse.  All other operations are
  // provided by the forward model.
  ThyraProductME_ROL_DynamicConstraint(const Teuchos::RCP<Tempus::Integrator<Real>> & forward_integrator,
    const Teuchos::RCP<Tempus::Integrator<Real>> & adjoint_integrator,
    const Teuchos::RCP<Thyra::ModelEvaluator<Real>> & modelAdjoin,
    Teuchos::ParameterList& piroParams,
    Teuchos::EVerbosityLevel verbLevel= Teuchos::VERB_HIGH,
    Teuchos::RCP<ROL_ObserverBase<Real>> observer = Teuchos::null);

  virtual ~ThyraProductME_ROL_DynamicConstraint() {}

  void setNumResponses(int num_responses) {
    num_responses_ = num_responses;
  }

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
ThyraProductME_ROL_DynamicConstraint<Real>::ThyraProductME_ROL_DynamicConstraint(
  const Teuchos::RCP<Tempus::Integrator<Real>>& integrator,
  Teuchos::ParameterList& piroParams,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> observer) :
  integrator_(integrator),
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  observer_(observer) {
  model_   = integrator->getStepper()->getModel();
  stepper_ = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(integrator->getStepper());
  Ju_      = model_->create_W();
  Ju_op_   = model_->create_W_op();
  Jz_op_   = model_->create_DfDp_op(0);
  // Set all other member variables to null.
  modelAdjoint_   = integrator->getAdjointModel();
  stepperAdjoint_ = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(integrator->getStepper());
  adjointJu_      = modelAdjoint_->create_W();
  adjointJu_op_   = modelAdjoint_->create_W_op();
  // Set adjoint flag to false.
  usingAdjoint_ = false;
  num_responses_ = -1;
}


//----------------------------------------------------------------------------
// Constructor with forward and adjoint models.

template<class Real>
ThyraProductME_ROL_DynamicConstraint<Real>::ThyraProductME_ROL_DynamicConstraint(
  const Teuchos::RCP<Tempus::Integrator<Real>>& forward_integrator,
  const Teuchos::RCP<Tempus::Integrator<Real>>& adjoint_integrator,
  const Teuchos::RCP<Thyra::ModelEvaluator<Real>>& modelAdjoint,
  Teuchos::ParameterList& piroParams,
  Teuchos::EVerbosityLevel verbLevel,
  Teuchos::RCP<ROL_ObserverBase<Real>> observer) :
  optParams_(piroParams.sublist("Optimization Status")),
  out_(Teuchos::VerboseObjectBase::getDefaultOStream()),
  verbosityLevel_(verbLevel),
  observer_(observer) {
  model_          = forward_integrator->getStepper()->getModel();
  modelAdjoint_   = modelAdjoint;
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
void ThyraProductME_ROL_DynamicConstraint<Real>::value(ROL::Vector<Real> &c,
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
void ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_uo(ROL::Vector<Real> &jv,
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

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_uo, the norm of rtv_v is not a number." << std::endl);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_un, the norm of rtv_v is not a number." << std::endl);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyJacobian_z, the norm of rtv_v is not a number." << std::endl);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_uo, the norm of rtv_v is not a number." << std::endl);

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
    rtv_jv.getVector()->assign(0.0);
    //Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_un, the norm of rtv_v is not a number." << std::endl);

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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyAdjointJacobian_z, the norm of rtv_v is not a number." << std::endl);

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
void ThyraProductME_ROL_DynamicConstraint<Real>::applyInverseJacobian_un(ROL::Vector<Real> &jv,
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

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyInverseJacobian_un, the norm of rtv_v is not a number." << std::endl);

  Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Real>>> x(2);
  Teuchos::Array<Real> t(2);
  x[0] = rtv_unew.getVector(); t[0] = ts.t[1]; // note: x[0] is the NEW state
  x[1] = rtv_uold.getVector(); t[1] = ts.t[0]; // note: x[1] is the OLD state

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
  const ROL::ThyraVector<Real>& rtv_uold = dynamic_cast<const ROL::ThyraVector<Real>& >(uold);
  const ROL::ThyraVector<Real>& rtv_unew = dynamic_cast<const ROL::ThyraVector<Real>& >(unew);
  const ROL::ThyraVector<Real>&    rtv_z = dynamic_cast<const ROL::ThyraVector<Real>& >(z);

  TEUCHOS_TEST_FOR_EXCEPTION( SMT::isnaninf(rtv_v.getVector()->norm_2()), std::logic_error, 
    "ThyraProductME_ROL_DynamicConstraint<Real>::applyInverseAdjointJacobian_un, the norm of rtv_v is not a number." << std::endl);

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
    rtv_jv.getVector()->assign(0.0);
    adjointJu_->solve(Thyra::NOTRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
  }
} // applyInverseAdjointJacobian_un


} // namespace Piro


#endif
