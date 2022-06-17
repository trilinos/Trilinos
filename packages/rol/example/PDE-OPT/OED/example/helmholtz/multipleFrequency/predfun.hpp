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

#ifndef PDEOPT_PREDFUN_H
#define PDEOPT_PREDFUN_H

#include "model.hpp"

template <class Real>
class PredFun : public ROL::Objective<Real> {
private:
  const ROL::Ptr<ROL::Constraint<Real>> model_;
  const ROL::Ptr<ROL::Vector<Real>> vec_;
  const ROL::Ptr<ROL::Vector<Real>> vdual_;

public:
  PredFun(const ROL::Ptr<ROL::Constraint<Real>> &model,
          const ROL::Ptr<ROL::Vector<Real>> &vec)
    : model_(model), vec_(vec), vdual_(vec_->dual().clone()) {}

  void setParameter(const std::vector<Real> &param) {
    ROL::Objective<Real>::setParameter(param);
    model_->setParameter(param);
  }

  void update(const ROL::Vector<Real> &x, bool flag = true, int iter = -1) {
    model_->update(x,flag,iter);
  }

  Real value(const ROL::Vector<Real> &x, Real &tol ) {
    model_->value(*vdual_,x,tol);
    return vdual_->apply(*vec_);
  }

  void gradient(ROL::Vector<Real> &g,
                const ROL::Vector<Real> &x,
                Real &tol) {
    model_->applyAdjointJacobian(g,*vec_,x,tol);
  }

}; // class PredFun

#endif
