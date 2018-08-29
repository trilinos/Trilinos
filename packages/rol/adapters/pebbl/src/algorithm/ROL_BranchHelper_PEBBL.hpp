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

#ifndef ROL_BRANCHHELPER_PEBBL_H
#define ROL_BRANCHHELPER_PEBBL_H

#include "ROL_Vector.hpp"
#include "ROL_Transform_PEBBL.hpp"

/** @ingroup func_group
    \class ROL::BranchHelper_PEBBL
    \brief Defines the pebbl branch index interface.

    ---
*/


namespace ROL {

template <class Real>
class BranchHelper_PEBBL {
private:
  Ptr<const Vector<Real>> getVector(const Vector<Real> &xs ) const {
    try {
      return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return makePtrFromRef(xs);
    }
  }

public:
  virtual ~BranchHelper_PEBBL(void) {}

  BranchHelper_PEBBL(void) {}

  BranchHelper_PEBBL(const BranchHelper_PEBBL &con) {}

  virtual int getMyIndex(const Vector<Real> &x) const = 0;
  virtual void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                            const Vector<Real> &x) const = 0;

  int getIndex(const Vector<Real> &x) const {
    Ptr<const Vector<Real>> xp = getVector(x);
    return getMyIndex(*xp);
  }

  void getNumFrac(int &nfrac, Real &integralityMeasure,
                  const Vector<Real> &x) const {
    Ptr<const Vector<Real>> xp = getVector(x);
    getMyNumFrac(nfrac, integralityMeasure, *xp);
  }

  virtual Ptr<Transform_PEBBL<Real>> createTransform(void) const = 0;
}; // class BranchHelper_PEBBL

} // namespace ROL

#endif
