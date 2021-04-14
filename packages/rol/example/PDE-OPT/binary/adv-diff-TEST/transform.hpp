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

#ifndef ROL_ADVDIFFTEST_INTEGERTRANSFORMATION_H
#define ROL_ADVDIFFTEST_INTEGERTRANSFORMATION_H

#include "ROL_PEBBL_StdIntegerTransformation.hpp"
#include "ROL_PEBBL_TpetraIntegerTransformation.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class StdAdvDiffIntegerTransformation : public ROL::PEBBL::StdIntegerTransformation<Real> {
private:
  ROL::Ptr<ROL::StdVector<Real>> getParameter(ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<ROL::StdVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
    }
  }

public:
  StdAdvDiffIntegerTransformation(void)
    : ROL::PEBBL::StdIntegerTransformation<Real>() {}

  StdAdvDiffIntegerTransformation(const StdAdvDiffIntegerTransformation &T)
    : ROL::PEBBL::StdIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::StdIntegerTransformation<Real>::fixValues(*getParameter(c),zero);
  }

}; // class StdAdvDiffIntegerTransformation

template <class Real>
class TpetraAdvDiffIntegerTransformation : public ROL::PEBBL::TpetraIntegerTransformation<Real> {
private:
  ROL::Ptr<ROL::TpetraMultiVector<Real>> getData(ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<ROL::TpetraMultiVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<PDE_OptVector<Real>&>(x).getField();
    }
  }

public:
  TpetraAdvDiffIntegerTransformation(void)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>() {}

  TpetraAdvDiffIntegerTransformation(const TpetraAdvDiffIntegerTransformation &T)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::TpetraIntegerTransformation<Real>::fixValues(*getData(c),zero);
  }

}; // class TpetraAdvDiffIntegerTransformation

#endif
