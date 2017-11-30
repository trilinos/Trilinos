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

#ifndef ROL_RISKLESSOBJECTIVE_HPP
#define ROL_RISKLESSOBJECTIVE_HPP

#include "ROL_RiskVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

template<class Real>
class RiskLessObjective : public Objective<Real> {
private:
  const ROL::Ptr<Objective<Real> > obj_;

public:
  RiskLessObjective(const ROL::Ptr<Objective<Real> > &obj) : obj_(obj) {}

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    ROL::Ptr<const Vector<Real> > x0
      = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    obj_->update(*x0,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const Vector<Real> > x0
      = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    return obj_->value(*x0,tol);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<Vector<Real> > g0
      = dynamic_cast<RiskVector<Real>&>(g).getVector();
    ROL::Ptr<const Vector<Real> > x0
      = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    obj_->gradient(*g0,*x0,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<Vector<Real> > hv0
      = dynamic_cast<RiskVector<Real>&>(hv).getVector();
    ROL::Ptr<const Vector<Real> > v0
      = dynamic_cast<const RiskVector<Real>&>(v).getVector();
    ROL::Ptr<const Vector<Real> > x0
      = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    obj_->hessVec(*hv0,*v0,*x0,tol);
  }

  void precond( Vector<Real> &Pv, const Vector<Real> &v,
                const Vector<Real> &x, Real &tol ) {
    ROL::Ptr<Vector<Real> > Pv0
      = dynamic_cast<RiskVector<Real>&>(Pv).getVector();
    ROL::Ptr<const Vector<Real> > v0
      = dynamic_cast<const RiskVector<Real>&>(v).getVector();
    ROL::Ptr<const Vector<Real> > x0
      = dynamic_cast<const RiskVector<Real>&>(x).getVector();
    obj_->precond(*Pv0,*v0,*x0,tol);
  }
};

}

#endif
