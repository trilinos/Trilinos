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

#ifndef ROL_TRANSFORM_PEBBL_H
#define ROL_TRANSFORM_PEBBL_H

#include "ROL_Constraint.hpp"

/** @ingroup func_group
    \class ROL::Transform_PEBBL
    \brief Defines the pebbl transform operator interface.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {

template <class Real>
class Transform_PEBBL : public Constraint<Real> {
protected:
  std::map<int,Real> map_;

public:
  Transform_PEBBL(void) {}

  Transform_PEBBL(const Transform_PEBBL &T) : map_(T.map_) {}

  void applyAdjointJacobian(Vector<Real> &ajv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                            Real &tol) {
    applyJacobian(ajv,v,x,tol);
  }

  void applyAdjointHessian(Vector<Real> &ahuv,
                     const Vector<Real> &u,
                     const Vector<Real> &v,
                     const Vector<Real> &x,
                           Real &tol) {
    ahuv.zero();
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

}; // class Transform_PEBBL

} // namespace ROL

#endif
