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

#ifndef ROL_PDEOPT_MULTIMAT_BRANCHHELPER_H
#define ROL_PDEOPT_MULTIMAT_BRANCHHELPER_H

#include "ROL_PEBBL_TpetraBranchHelper.hpp"
#include "transform.hpp"

template <class Real>
class TpetraMultiMatBranchHelper : public ROL::PEBBL::TpetraBranchHelper<Real> {
private:
  using ROL::PEBBL::TpetraBranchHelper<Real>::tol_;
  using ROL::PEBBL::TpetraBranchHelper<Real>::method_;
  using ROL::PEBBL::TpetraBranchHelper<Real>::getConstData;

  ROL::Ptr<const ROL::TpetraMultiVector<Real>> getData(const ROL::Vector<Real> &x) const {
    try {
      return ROL::dynamicPtrCast<const ROL::TpetraMultiVector<Real>>(dynamic_cast<const ROL::PartitionedVector<Real>&>(x).get(0));
    }
    catch (std::exception &e) {
      return ROL::dynamicPtrCast<const PDE_OptVector<Real>>(dynamic_cast<const ROL::PartitionedVector<Real>&>(x).get(0))->getField();
    }
  }

  int getIndex_nn(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    Teuchos::ArrayView<const Real> xview = (getConstData(x)->getData(0))();
    Teuchos::ArrayView<const Real> gview = (getConstData(g)->getData(0))();
    const Real one(1);
    Real maxD(ROL::ROL_NINF<Real>()), Li(0), Ui(0), mini(0), d(0), nd(0);
    int index = 0, size = xview.size(), cnt = 0;
    int nx = 2*std::sqrt(size/2); // nx = 2*ny
    int ny =   std::sqrt(size/2); // nc = nx*ny = 2*ny*ny
    // Compute number of neighbors
    std::multimap<int,int> map;
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        d = xview[i+j*nx];
        cnt = 0;
        if (d > tol_ && d < one - tol_) {
          for (int k = -1; k < 2; ++k) {
            for (int l = -1; l < 2; ++l) {
              if (k!=0 || l!=0) {
                if ( (i+k) + (j+l)*nx > 0 && (i+k) + (j+l)*nx < size ) {
                  nd = xview[(i+k) + (j+l)*nx];
                  if (nd == one) cnt++;
                }
              }
            }
          }
          map.insert({cnt,i+j*nx});
        }
      }
    }
    // Find the most fractional element with the most neighbors equal to one
    auto it = map.rbegin();
    //auto it = map.begin();
    cnt = map.count(it->first);
    index = it->second;
    if (cnt > 1) {
      for (int i = 0; i < cnt; ++i) {
        Li   = gview[it->second] * (std::floor(xview[it->second]) - xview[it->second]);
        Ui   = gview[it->second] * (std::ceil( xview[it->second]) - xview[it->second]);
        mini = std::min(std::abs(Li),std::abs(Ui));
        if (mini > maxD) {
          maxD  = mini;
          index = i;
        }
        it++;
      }
    }
    return index;
  }

public:
  TpetraMultiMatBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::TpetraBranchHelper<Real>(tol, method) {}

  TpetraMultiMatBranchHelper(const TpetraMultiMatBranchHelper &BH)
    : ROL::PEBBL::TpetraBranchHelper<Real>(BH) {}

  //int getMyIndex(const ROL::Vector<Real> &x) const {
  int getMyIndex(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    if (method_ < 2) // Use Std implementation
      return ROL::PEBBL::TpetraBranchHelper<Real>::getMyIndex(*getData(x),*getData(g));
    return getMyIndex_nn(*getData(x),*getData(g));
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const ROL::Vector<Real> &x) const {
    // Use Std implementation
    ROL::PEBBL::TpetraBranchHelper<Real>::getMyNumFrac(nfrac, integralityMeasure, *getData(x));
  }

  ROL::Ptr<ROL::PEBBL::IntegerTransformation<Real>> createTransform(void) const {
    return ROL::makePtr<TpetraMultiMatIntegerTransformation<Real>>();
  }

}; // class TpetraMultiMatBranchHelper

#endif
