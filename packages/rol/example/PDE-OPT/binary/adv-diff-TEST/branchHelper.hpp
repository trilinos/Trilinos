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

#ifndef ROL_PDEOPT_BRANCHHELPER_PEBBL_H
#define ROL_PDEOPT_BRANCHHELPER_PEBBL_H

#include "ROL_StdBranchHelper_PEBBL.hpp"
#include "transform.hpp"

template <class Real>
class PDEOPT_BranchHelper_PEBBL : public ROL::StdBranchHelper_PEBBL<Real> {
private:
  ROL::Ptr<const ROL::StdVector<Real>> getParameter(const ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<const ROL::StdVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
    }
  }

public:
  PDEOPT_BranchHelper_PEBBL(const Real tol = 1e-6, const int method = 0)
    : ROL::StdBranchHelper_PEBBL<Real>(tol, method) {}

  PDEOPT_BranchHelper_PEBBL(const PDEOPT_BranchHelper_PEBBL &BH)
    : ROL::StdBranchHelper_PEBBL<Real>(BH) {}

  //int getMyIndex(const ROL::Vector<Real> &x) const {
  int getMyIndex(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    // Use Std implementation
    return ROL::StdBranchHelper_PEBBL<Real>::getMyIndex(*getParameter(x),*getParameter(g));
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const ROL::Vector<Real> &x) const {
    // Use Std implementation
    ROL::StdBranchHelper_PEBBL<Real>::getMyNumFrac(nfrac, integralityMeasure, *getParameter(x));
  }

  ROL::Ptr<ROL::Transform_PEBBL<Real>> createTransform(void) const {
    return ROL::makePtr<PDEOPT_Transform_PEBBL<Real>>();
  }

}; // class PDEOPT_BranchHelper_PEBBL

#endif
