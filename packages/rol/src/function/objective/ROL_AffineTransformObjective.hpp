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

#ifndef ROL_AFFINE_TRANSFORM_OBJECTIVE_H
#define ROL_AFFINE_TRANSFORM_OBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_SimController.hpp"

/** @ingroup func_group
    \class ROL::AffineTransformObjective
    \brief Compose an objective function with an affine transformation, i.e.,

    \f[ F(x) = f(Ax+b). \f]

*/


namespace ROL {

template <class Real>
class AffineTransformObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>>      obj_;
  const Ptr<LinearOperator<Real>> A_;
  const Ptr<Vector<Real>>         b_;

  Ptr<SimController<Real>> storage_;

  Ptr<Vector<Real>> primal_, dual_, Av_;

  const Ptr<const Vector<Real>> transform(const Vector<Real> &x) {
    bool isApplied = storage_->get(*primal_,Objective<Real>::getParameter());
    if (!isApplied) {
      Real tol = std::sqrt(ROL_EPSILON<Real>());
      A_->apply(*primal_,x,tol);
      primal_->plus(*b_);
      storage_->set(*primal_,Objective<Real>::getParameter());
    }
    return primal_;
  }

public:
  AffineTransformObjective(const Ptr<Objective<Real>>      &obj,
                           const Ptr<LinearOperator<Real>> &A,
                           const Ptr<Vector<Real>>         &b,
                           const Ptr<SimController<Real>>  &storage = nullPtr)
    : obj_(obj), A_(A), b_(b), storage_(storage) {
    primal_ = b->clone();
    Av_     = b->clone();
    dual_   = b->dual().clone();
    if (storage == nullPtr) {
      storage_ = makePtr<SimController<Real>>();
    }
  }

  virtual ~AffineTransformObjective() {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    storage_->objectiveUpdate(true);
    obj_->update(*transform(x),flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    return obj_->value(*transform(x),tol); 
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    obj_->gradient(*dual_,*transform(x),tol);
    A_->applyAdjoint(g,*dual_,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    A_->apply(*Av_,v,tol);
    obj_->hessVec(*dual_,*Av_,*transform(x),tol);
    A_->applyAdjoint(hv,*dual_,tol);
  }

}; // class AffineTransformObjective

} // namespace ROL

#endif // ROL_AFFINE_TRANSFORM_OBJECTIVE_H
