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

#ifndef ROL_PEBBL_BRANCHHELPER_H
#define ROL_PEBBL_BRANCHHELPER_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_PEBBL_IntegerTransformation.hpp"
#include "ROL_PEBBL_MixedVector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::BranchHelper
    \brief Defines the pebbl branch index interface.

    ---
*/

namespace ROL {
namespace PEBBL {

template <class Real>
class BranchHelper {
protected:
  Ptr<const Vector<Real>> getOptVector(const Vector<Real> &xs ) const {
    try {
      return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return makePtrFromRef(xs);
    }
  }

  Ptr<const Vector<Real>> getIntegerVector(const Vector<Real> &xs) const {
    try {
      return dynamic_cast<const MixedVector<Real>&>(*getOptVector(xs)).getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(xs);
    }
  }

public:
  virtual ~BranchHelper(void) {}
  BranchHelper(void) {}
  BranchHelper(const BranchHelper &con) {}

  virtual int getIndex(const Vector<Real> &x, const Vector<Real> &g) const = 0;
  virtual void getNumFrac(int &nfrac, Real &integralityMeasure, const Vector<Real> &x) const = 0;
  virtual Ptr<IntegerTransformation<Real>> createTransform(void) const = 0;

}; // class BranchHelper

} // namespace PEBBL
} // namespace ROL

#endif
