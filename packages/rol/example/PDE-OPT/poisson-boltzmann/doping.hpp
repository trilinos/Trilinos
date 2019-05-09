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

#ifndef DOPING_HPP
#define DOPING_HPP

#include "../TOOLS/fe.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_MultiVector.hpp"
#include "ROL_Types.hpp"

template <class Real>
class Doping {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<Intrepid::FieldContainer<Real> > cellNodes_;
  const ROL::Ptr<Intrepid::FieldContainer<int> > cellDofs_;
  const Teuchos::Array<int> cellIds_;
  Real a_, b_, wx_, wy_;

public:
  Doping(const ROL::Ptr<FE<Real> > &fe,
         const ROL::Ptr<Intrepid::FieldContainer<Real> > & cellNodes,
         const ROL::Ptr<Intrepid::FieldContainer<int> > & cellDofs,
         const Teuchos::Array<int> & cellIds,
         Teuchos::ParameterList &parlist)
    : fe_(fe), cellNodes_(cellNodes), cellDofs_(cellDofs), cellIds_(cellIds) {
    // Get parameters from parameter list
    a_  = parlist.sublist("Problem").get("Desired Lower Doping Value", 0.3);
    b_  = parlist.sublist("Problem").get("Desired Upper Doping Value", 1.0);
    wx_ = parlist.sublist("Geometry").get("Width",0.6);
    wy_ = parlist.sublist("Geometry").get("Height",0.2);
  }

  Real evaluate(const std::vector<Real> &x) const {
    Real zero(0), dx(wx_/6.0), dy(wy_/4.0);
    Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
    bool inRegion1 = ((x[0] > zero-eps)   && (x[0] < dx+eps)  && (x[1] > wy_-dy-eps) && (x[1] < wy_+eps));
    bool inRegion2 = ((x[0] > wx_-dx-eps) && (x[0] < wx_+eps) && (x[1] > wy_-dy-eps) && (x[1] < wy_+eps));
    return (inRegion1 || inRegion2 ? b_ : a_);
  }

  void build(const ROL::Ptr<Tpetra::MultiVector<> > &dpVec) const {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
    ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,d);
    fe_->computeDofCoords(dofPoints, cellNodes_);
    
    std::vector<Real> coord(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < f; ++j) {
        int fidx = (*cellDofs_)(cellIds_[i],j);
        if (dpVec->getMap()->isNodeGlobalElement(fidx)) {
          for (int k = 0; k < d; ++k) {
            coord[k] = (*dofPoints)(i,j,k);
          }
          dpVec->replaceGlobalValue(fidx,
                                    0,
                                    evaluate(coord));
        }
      }
    }
  }
};

template <class Real>
class DopingBounds {
private:
  const ROL::Ptr<FE<Real> > fe_;
  const ROL::Ptr<Intrepid::FieldContainer<Real> > cellNodes_;
  const ROL::Ptr<Intrepid::FieldContainer<int> > cellDofs_;
  const Teuchos::Array<int> cellIds_;
  Real a_, b_, wx_, wy_, lo_, hi_;
  bool useConstant_;

public:
  DopingBounds(const ROL::Ptr<FE<Real> > &fe,
               const ROL::Ptr<Intrepid::FieldContainer<Real> > & cellNodes,
               const ROL::Ptr<Intrepid::FieldContainer<int> > & cellDofs,
               const Teuchos::Array<int> & cellIds,
               Teuchos::ParameterList &parlist)
    : fe_(fe), cellNodes_(cellNodes), cellDofs_(cellDofs), cellIds_(cellIds) {
    // Get parameters from parameter list
    a_  = parlist.sublist("Problem").get("Desired Lower Doping Value", 0.3);
    b_  = parlist.sublist("Problem").get("Desired Upper Doping Value", 1.0);
    wx_ = parlist.sublist("Geometry").get("Width",0.6);
    wy_ = parlist.sublist("Geometry").get("Height",0.2);
    lo_ = parlist.sublist("Problem").get("Lower Bound Value",  0.0);
    hi_ = parlist.sublist("Problem").get("Upper Bound Value", 10.0);
    useConstant_ = parlist.sublist("Problem").get("Use Constant Bounds",true);
  }

  Real evaluateLower(const std::vector<Real> &x) const {
    if (!useConstant_) {
      const Real zero(0), two(2), four(4), dx(wx_/6.0);
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      const bool inRegion1 = ((x[0] > zero-eps)   && (x[0] < dx+eps)      && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      const bool inRegion2 = ((x[0] > two*dx-eps) && (x[0] < four*dx+eps) && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      const bool inRegion3 = ((x[0] > wx_-dx-eps) && (x[0] < wx_+eps)     && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      return (inRegion1 || inRegion3 ? b_ : (inRegion2 ? a_ : lo_));
    }
    return lo_;
  }

  Real evaluateUpper(const std::vector<Real> &x) const {
    if (!useConstant_) {
      const Real zero(0), two(2), four(4), dx(wx_/6.0);
      const Real eps = std::sqrt(ROL::ROL_EPSILON<Real>());
      const bool inRegion1 = ((x[0] > zero-eps)   && (x[0] < dx+eps)      && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      const bool inRegion2 = ((x[0] > two*dx-eps) && (x[0] < four*dx+eps) && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      const bool inRegion3 = ((x[0] > wx_-dx-eps) && (x[0] < wx_+eps)     && (x[1] > wy_-eps) && (x[1] < wy_+eps));
      return (inRegion1 || inRegion3 ? b_ : (inRegion2 ? a_ : hi_));
    }
    return hi_;
  }

  void build(const ROL::Ptr<Tpetra::MultiVector<> > &loVec,
             const ROL::Ptr<Tpetra::MultiVector<> > &hiVec) const {
    int c = fe_->gradN()->dimension(0);
    int f = fe_->gradN()->dimension(1);
    int d = fe_->gradN()->dimension(3);
    ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints =
      ROL::makePtr<Intrepid::FieldContainer<Real>>(c,f,d);
    fe_->computeDofCoords(dofPoints, cellNodes_);
    
    std::vector<Real> coord(d);
    for (int i = 0; i < c; ++i) {
      for (int j = 0; j < f; ++j) {
        int fidx = (*cellDofs_)(cellIds_[i],j);
        if (loVec->getMap()->isNodeGlobalElement(fidx)) {
          for (int k = 0; k < d; ++k) {
            coord[k] = (*dofPoints)(i,j,k);
          }
          loVec->replaceGlobalValue(fidx,
                                    0,
                                    evaluateLower(coord));
        }
        if (hiVec->getMap()->isNodeGlobalElement(fidx)) {
          for (int k = 0; k < d; ++k) {
            coord[k] = (*dofPoints)(i,j,k);
          }
          hiVec->replaceGlobalValue(fidx,
                                    0,
                                    evaluateUpper(coord));
        }
      }
    }
  }
};

#endif
