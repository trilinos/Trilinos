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

#ifndef ROL_CONSTRAINT_PEBBL_H
#define ROL_CONSTRAINT_PEBBL_H

#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::Constraint_PEBBL
    \brief Defines the pebbl constraint operator interface.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is a ROL::StdVector.

    ---
*/


namespace ROL {

template <class Real>
class Constraint_PEBBL : public Constraint<Real> {
private:
  std::map<int,Real> map_;

  Ptr<std::vector<Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<StdVector<Real>&>(x).getVector();
  }

  Ptr<const std::vector<Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const StdVector<Real>&>(x).getVector();
  }

public:
  Constraint_PEBBL(void) {}

  Constraint_PEBBL(const Constraint_PEBBL &con) : map_(con.map_) {}

  void value(Vector<Real> &c,
       const Vector<Real> &x,
             Real &tol) {
    Ptr<std::vector<Real>> cval = getData(c);
    typename std::map<int,Real>::iterator it;
    int cnt = 0;
    for (it=map_.begin(); it!=map_.end(); ++it) {
      (*cval)[cnt] = x.dot(*x.basis(it->first))-it->second;
      cnt++;
    }
  }

  void applyJacobian(Vector<Real> &jv,
               const Vector<Real> &v,
               const Vector<Real> &x,
                     Real &tol) {
    Ptr<std::vector<Real>> jval = getData(jv);
    typename std::map<int,Real>::iterator it;
    int cnt = 0;
    for (it=map_.begin(); it!=map_.end(); ++it) {
      (*jval)[cnt] = v.dot(*x.basis(it->first));
      cnt++;
    }
  }

  void applyAdjointJacobian(Vector<Real> &ajv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                            Real &tol) {
    Ptr<const std::vector<Real>> vval = getConstData(v);
    typename std::map<int,Real>::iterator it;
    int cnt = 0;
    for (it=map_.begin(); it!=map_.end(); ++it) {
      ajv.axpy((*vval)[cnt],*x.basis(it->first));
      cnt++;
    }
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                     const Vector<Real> &u,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
  }

  Ptr<Vector<Real> > makeConstraintVector(void) const {
    return makePtr<StdVector<Real>>(makePtr<std::vector<Real>>(map_.size(),0));
  }

  bool isEmpty(void) const {
    return map_.empty();
  }

  void reset(void) {
    map_.clear();
  }

  void add(const std::map<int,Real> &in) {
    map_.insert(in.begin(),in.end());
  }

  void add(const std::pair<int,Real> &in) {
    map_.insert(in);
  }

}; // class Constraint_PEBBL

} // namespace ROL

#endif
