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
