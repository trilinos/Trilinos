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

#ifndef ROL_SLACKLESSOBJECTIVE_HPP
#define ROL_SLACKLESSOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
 *  \class ROL::SlacklessObjective
 *  \brief This class strips out the slack variables from objective evaluations
 *         to create the new objective  \f$ F(x,s) = f(x) \f$
 */

namespace ROL {

template<class Real> 
class SlacklessObjective : public Objective<Real> {
private: 
  ROL::Ptr<Objective<Real> > obj_;

  ROL::Ptr<Vector<Real> > getOpt( Vector<Real> &xs ) {
    return dynamic_cast<PartitionedVector<Real>&>(xs).get(0);
  }

  ROL::Ptr<const Vector<Real> > getOpt( const Vector<Real> &xs ) {
    return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
  }

  void zeroSlack( Vector<Real> &x ) {
    PartitionedVector<Real> &xpv
      = dynamic_cast<PartitionedVector<Real>&>(x);
    const int nvec = static_cast<int>(xpv.numVectors());
    for (int i = 1; i < nvec; ++i) {
      xpv.get(i)->zero();
    }
  } 

public:
  SlacklessObjective( const ROL::Ptr<Objective<Real> > &obj ) : obj_(obj) {}
  ~SlacklessObjective() {}
 
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update( *getOpt(x), flag, iter );
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    return obj_->value( *getOpt(x), tol );
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    zeroSlack(g);
    obj_->gradient(*getOpt(g),*getOpt(x),tol);
  }

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    return obj_->dirDeriv(*getOpt(x),*getOpt(d),tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    zeroSlack(hv);
    obj_->hessVec(*getOpt(hv),*getOpt(v),*getOpt(x),tol);     
  }

  void invHessVec( Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    zeroSlack(ihv);
    obj_->invHessVec( *getOpt(ihv), *getOpt(v), *getOpt(x), tol );
  }

  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    zeroSlack(Pv);
    obj_->precond( *getOpt(Pv), *getOpt(v), *getOpt(x), tol );
  }

// Definitions for parametrized (stochastic) objective functions
public:
  void setParameter(const std::vector<Real> &param) {
    Objective<Real>::setParameter(param);
    obj_->setParameter(param);
  }
}; // class SlacklessObjective 

} // namespace ROL

#endif // ROL__SLACKLESSOBJECTIVE_HPP

