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

#ifndef ROL_TEUCHOSBRANCHHELPER_PEBBL_H
#define ROL_TEUCHOSBRANCHHELPER_PEBBL_H

#include "ROL_TeuchosVector.hpp"
#include "ROL_BranchHelper_PEBBL.hpp"
#include "ROL_TeuchosTransform_PEBBL.hpp"

/** @ingroup func_group
    \class ROL::TeuchosBranchHelper_PEBBL
    \brief Defines the pebbl branch index interface for TeuchosVectors.

    ---
*/


namespace ROL {

template <class Ordinal, class Real>
class TeuchosBranchHelper_PEBBL : public BranchHelper_PEBBL<Real> {
private:
  const Real tol_;
  const int method_;

  using BranchHelper_PEBBL<Real>::getVector;
  Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(*getVector(x)).getVector();
  }

  // Branching based on distance to integer
  int getMyIndex_D(const Vector<Real> &x, const Vector<Real> &lam,
                   Objective<Real> &obj, Constraint<Real> &con) const {
    // Get index closest to 0.5
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> xval = getConstData(x);
    int index = 0;
    Real minD(0), minX(ROL_INF<Real>()), half(0.5);
    int size = xval->length();
    for (int i = 0; i < size; ++i) {
      Real x  = (*xval)(i);
      Real fx = std::floor(x);
      Real cx = std::ceil(x);
      minD    = std::min(x-fx,cx-x);
      if (std::abs(minD-half) < minX) {
        minX = std::abs(minD-half);
        index = i;
      }
    }
    return index;
  }

  // Branching based on directional derivatives (similar to pseudo costs).
  int getMyIndex_PC(const Vector<Real> &x, const Vector<Real> &lam,
                    Objective<Real> &obj, Constraint<Real> &con) const {
    Real tol = static_cast<Real>(1e-8);
    Ptr<Vector<Real>> g = x.dual().clone();
    Ptr<Vector<Real>> J = x.dual().clone();
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> xval = getConstData(x);
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> gval = getConstData(*g);
    obj.gradient(*g,x,tol);
    con.applyAdjointJacobian(*J,lam,x,tol);
    g->plus(*J);
    Real maxD(ROL_NINF<Real>()), Li(0), Ui(0), mini(0);
    int index = 0, size = gval->length();
    for (int i = 0; i < size; ++i) {
      Li   = (*gval)[i] * (std::floor((*xval)[i]) - (*xval)[i]);
      Ui   = (*gval)[i] * (std::ceil((*xval)[i])  - (*xval)[i]);
      mini = std::min(std::abs(Li),std::abs(Ui));
      if (mini > maxD) {
        maxD  = mini;
        index = i;
      }
    }
    return index;
  }

public:
  TeuchosBranchHelper_PEBBL(const Real tol = 1e-6, const int method = 0)
    : tol_(tol), method_(method) {}

  TeuchosBranchHelper_PEBBL(const TeuchosBranchHelper_PEBBL &BH)
    : tol_(BH.tol_), method_(BH.method_) {}

  int getMyIndex(const Vector<Real> &x, const Vector<Real> &lam,
                 Objective<Real> &obj, Constraint<Real> &con) const {
    int index(0);
    if (method_ == 1) {
      index = getMyIndex_D(x,lam,obj,con);
    }
    else {
      index = getMyIndex_PC(x,lam,obj,con);
    }
    return index;
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const Vector<Real> &x) const {
    // Return number of fractional variables and the
    // sum of the distance to integer for the input vector
    Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> xval = getConstData(x);
    nfrac = 0;
    integralityMeasure = static_cast<Real>(0);
    Real minD(0);
    int size = xval->length();
    for (int i = 0; i < size; ++i) {
      Real x  = (*xval)(i);
      Real fx = std::floor(x);
      Real cx = std::ceil(x);
      minD    = std::min(x-fx,cx-x);
      integralityMeasure += minD;
      if (minD > tol_) {
        nfrac++;
      }
    }
  }

  Ptr<Transform_PEBBL<Real>> createTransform(void) const {
    return makePtr<TeuchosTransform_PEBBL<Ordinal,Real>>();
  }

}; // class StdBranchHelper_PEBBL

} // namespace ROL

#endif
