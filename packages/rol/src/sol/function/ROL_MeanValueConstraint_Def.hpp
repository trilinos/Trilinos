// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVALUECONSTRAINT_DEF_HPP
#define ROL_MEANVALUECONSTRAINT_DEF_HPP

namespace ROL {

template<typename Real>
MeanValueConstraint<Real>::MeanValueConstraint( const Ptr<Constraint<Real>>      &con,
                                                const Ptr<SampleGenerator<Real>> &sampler)
  : con_(con) {
  std::vector<Real> param = computeSampleMean(sampler);
  con_->setParameter(param);
}

template<typename Real>
void MeanValueConstraint<Real>::update( const Vector<Real> &x, bool flag, int iter ) {
  con_->update(x,flag,iter);
}

template<typename Real>
void MeanValueConstraint<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  con_->update(x,type,iter);
}

template<typename Real>
void MeanValueConstraint<Real>::value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) {
  con_->value(c,x,tol);
}

template<typename Real>
void MeanValueConstraint<Real>::applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  con_->applyJacobian(jv,v,x,tol);
}

template<typename Real>
void MeanValueConstraint<Real>::applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  con_->applyAdjointJacobian(ajv,v,x,tol);
}

template<typename Real>
void MeanValueConstraint<Real>::applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
  con_->applyAdjointHessian(ahuv,u,v,x,tol);
}

template<typename Real>
std::vector<Real> MeanValueConstraint<Real>::computeSampleMean(const Ptr<SampleGenerator<Real>> &sampler) const {
  // Compute mean value of inputs and set parameter in constraint
  int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
  std::vector<Real> loc(dim), mean(dim), pt(dim);
  Real wt(0);
  for (int i = 0; i < nsamp; i++) {
    pt = sampler->getMyPoint(i);
    wt = sampler->getMyWeight(i);
    for (int j = 0; j < dim; j++) {
      loc[j] += wt*pt[j];
    }
  }
  sampler->sumAll(&loc[0],&mean[0],dim);
  return mean;
}

}

#endif
