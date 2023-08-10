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

#ifndef ROL_TYPEBINDICATOROBJECTIVE_H
#define ROL_TYPEBINDICATOROBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"

/** @ingroup func_group
    \class ROL::TypeBIndicatorObjective
    \brief Provides the interface to evaluate the indicator function of linear constraints.

        ---
*/


namespace ROL {

template<typename Real>
class TypeBIndicatorObjective : public Objective<Real> {
private:
  const Ptr<PolyhedralProjection<Real>> proj_;
  const Ptr<Vector<Real>> res_;
  bool isInit_;
  Real tol_;
 
public:

  TypeBIndicatorObjective(const Ptr<BoundConstraint<Real>> &bnd)
    : proj_(makePtr<PolyhedralProjection<Real>>(bnd)),
      isInit_(true), tol_(0) {}

  TypeBIndicatorObjective(const Vector<Real>               &xprim,
                          const Vector<Real>               &xdual,
                          const Ptr<BoundConstraint<Real>> &bnd,
                          const Ptr<Constraint<Real>>      &con,
                          const Vector<Real>               &mul,
                          const Vector<Real>               &res,
                          ParameterList                    &list)
    : proj_(PolyhedralProjectionFactory<Real>(xprim,xdual,bnd,con,mul,res,list)),
      res_(res.clone()), isInit_(false) {}

  TypeBIndicatorObjective(const Ptr<PolyhedralProjection<Real>> &proj)
    : proj_(proj), res_(proj->getResidual()->clone()), isInit_(false) {}

  void initialize(const Vector<Real> &x) {
    if (!isInit_) {
      auto xz = x.clone(); xz->zero();
      Real tol(std::sqrt(ROL_EPSILON<Real>()));
      tol_ = static_cast<Real>(1e-2)*tol;
      proj_->getLinearConstraint()->value(*res_,*xz,tol);
      Real rnorm = res_->norm();
      if (rnorm > ROL_EPSILON<Real>()) tol_ *= rnorm;
      isInit_ = true;
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    initialize(x);
    const Real zero(0);
    bool isBndFeasible = proj_->getBoundConstraint()->isFeasible(x); 
    bool isConFeasible = true;
    if (res_ != nullPtr) {
      proj_->getLinearConstraint()->value(*res_,x,tol);
      if (res_->norm() > tol_) isConFeasible = false;
    }
    return (isBndFeasible && isConFeasible) ? zero : ROL_INF<Real>();
  }

  void prox( Vector<Real> &Pv, const Vector<Real> &v, Real t, Real &tol){
    Pv.set(v); proj_->project(Pv);
  }
}; // class TypeBIndicatorObjective

} // namespace ROL

#endif
