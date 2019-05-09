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

#ifndef ROL_STDOBJECTIVE_H
#define ROL_STDOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::StdObjective
    \brief Specializes the ROL::Objective interface for objective functions
    that operate on ROL::StdVector's.
*/


namespace ROL {

template <class Real>
class StdObjective : public virtual Objective<Real> {
public:
  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}

  using Objective<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    update(*(xs.getVector()),flag,true);
  }

  virtual Real value( const std::vector<Real> &x, Real &tol ) = 0;

  using Objective<Real>::value;
  Real value( const Vector<Real> &x, Real &tol ) {
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    return value(*(xs.getVector()),tol);
  }

  virtual void gradient( std::vector<Real> &g, const std::vector<Real> &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdObjective): gradient not implemented!");
  }

  using Objective<Real>::gradient;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    StdVector<Real> gs = dynamic_cast<StdVector<Real>&>(g);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    try {
      gradient(*(gs.getVector()),*(xs.getVector()),tol);
    }
    catch (std::exception &e) {
      Objective<Real>::gradient(g,x,tol);
    }
  }

  virtual Real dirDeriv( const std::vector<Real> &x, const std::vector<Real> &d, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdObjective): dirDeriv not implemented!");
  }

  using Objective<Real>::dirDeriv;
  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    const StdVector<Real> ds = dynamic_cast<const StdVector<Real>&>(d);
    try {
      return dirDeriv(*(xs.getVector()),*(ds.getVector()),tol);
    }
    catch (std::exception &e) {
      return Objective<Real>::dirDeriv(x,d,tol);
    }
  }

  virtual void hessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdObjective): hessVec not implemented!");
  }

  using Objective<Real>::hessVec;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    try {
      StdVector<Real> hvs = dynamic_cast<StdVector<Real>&>(hv);
      const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
      const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
      hessVec(*(hvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
    }
    catch (std::exception &e) {
      Objective<Real>::hessVec(hv,v,x,tol);
    }
  }

  virtual void invHessVec( std::vector<Real> &hv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      ">>> ERROR (ROL::StdObjective): invHessVec not implemented!");
  }

  using Objective<Real>::invHessVec;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    StdVector<Real> hvs = dynamic_cast<StdVector<Real>&>(hv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    invHessVec(*(hvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
  }

  virtual void precond( std::vector<Real> &Pv, const std::vector<Real> &v, const std::vector<Real> &x, Real &tol ) {
    Pv.assign(v.begin(),v.end());
  }

  using Objective<Real>::precond;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    StdVector<Real> Pvs = dynamic_cast<StdVector<Real>&>(Pv);
    const StdVector<Real> vs = dynamic_cast<const StdVector<Real>&>(v);
    const StdVector<Real> xs = dynamic_cast<const StdVector<Real>&>(x);
    precond(*(Pvs.getVector()),*(vs.getVector()),*(xs.getVector()),tol);
  }
};

} // namespace ROL

#endif
