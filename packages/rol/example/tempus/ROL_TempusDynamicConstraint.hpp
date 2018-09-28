// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
  ROL::Ptr<const Thyra::ModelEvaluator<Real>> model_;
  ROL::Ptr<Tempus::StepperOptimizationInterface<Real>> stepper_;
  ROL::Ptr<Thyra::LinearOpWithSolveBase<Real>> Ju_;
  ROL::Ptr<Thyra::LinearOpBase<Real>> Ju_op_;
  ROL::Ptr<Thyra::LinearOpBase<Real>> Jz_op_;

public:

  TempusDynamicConstraint(const ROL::Ptr<Tempus::Integrator<Real>> & integrator);

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


template<class Real>
TempusDynamicConstraint<Real>::TempusDynamicConstraint(
  const ROL::Ptr<Tempus::Integrator<Real>> & integrator) {
  model_   = integrator->getStepper()->getModel();
  stepper_ = ROL::dynamicPtrCast<Tempus::StepperOptimizationInterface<Real>>(integrator->getStepper());
  Ju_      = model_->create_W();
  Ju_op_   = model_->create_W_op();
  Jz_op_   = model_->create_DfDp_op(0);
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

  int deriv_index = 1; // 1 = old state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
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

  int deriv_index = 0; // 0 = new state
  stepper_->computeStepJacobian(*Ju_op_, x, t, *(rtv_z.getVector()), 0, deriv_index);
  Ju_op_->apply(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), 1.0, 0.0);
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
  stepper_->computeStepSolver(*Ju_, x, t, *(rtv_z.getVector()), 0);
  Ju_->solve(Thyra::TRANS, *(rtv_v.getVector()), rtv_jv.getVector().ptr(), solveCriteria);
} // applyInverseJacobian_un


} // namespace ROL


#endif
