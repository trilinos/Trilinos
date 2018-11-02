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

#ifndef ROL_LINEAR_OBJECTIVE_SIMOPT_H
#define ROL_LINEAR_OBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::LinearObjective_SimOpt
    \brief Provides the interface to evaluate linear objective functions.

    This class implements the linear objective function
    \f[
       f(x) = \langle c, x\rangle_{\mathcal{X}^*,\mathcal{X}}
    \f]
    for fixed \f$c\in\mathcal{X}^*\f$.

    ---
*/


namespace ROL {

template <class Real>
class LinearObjective_SimOpt : public Objective_SimOpt<Real> {
private:
  const Ptr<const Vector<Real>> simcost_, optcost_;

public:
  LinearObjective_SimOpt(const Ptr<const Vector<Real>> &simcost = nullPtr,
                         const Ptr<const Vector<Real>> &optcost = nullPtr)
    : simcost_(simcost), optcost_(optcost) {}

  Real value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    Real valu(0), valz(0);
    if (simcost_ != nullPtr) {
      valu = u.dot(simcost_->dual());
    }
    if (optcost_ != nullPtr) {
      valz = z.dot(optcost_->dual());
    }
    return valu + valz;
  }

  void gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (simcost_ != nullPtr) {
      g.set(*simcost_);
    }
    else {
      g.zero();
    }
  }

  void gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    if (optcost_ != nullPtr) {
      g.set(*optcost_);
    }
    else {
      g.zero();
    }
  }

  void hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

  void hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
    hv.zero();
  }

}; // class LinearObjective_SimOpt

} // namespace ROL

#endif
