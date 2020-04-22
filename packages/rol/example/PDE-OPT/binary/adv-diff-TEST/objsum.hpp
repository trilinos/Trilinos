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

#ifndef ROL_OBJSUM_H
#define ROL_OBJSUM_H

#include "misfit_robj.hpp"

template <class Real>
class Sum_Objective : public ROL::Objective<Real> {
private:
  const ROL::Ptr<ROL::Objective<Real>> misfit_, penalty_, binary_;

  ROL::Ptr<ROL::Vector<Real>> xdual_;
  bool initialized_;

  const ROL::Ptr<std::ostream> stream_;
  const bool printToStream_;

  Real misCost_, penCost_, binCost_;

public:
  Sum_Objective(const ROL::Ptr<FEMdata<Real>> &fem_,
                const ROL::Ptr<ROL::Objective<Real>> &pen,
                const ROL::Ptr<ROL::Objective<Real>> &bin,
                ROL::ParameterList &list,
                ROL::Ptr<std::ostream> &stream = ROL::nullPtr,
                bool printToStream = false)
    : misfit_(ROL::makePtr<Misfit_Objective<Real>>(fem_,list)),
      penalty_(pen), binary_(bin),
      xdual_(ROL::nullPtr), initialized_(false),
      stream_(stream), printToStream_(printToStream) {
    misCost_ = list.sublist("Problem").get("State Cost",       1.0);
    penCost_ = list.sublist("Problem").get("Control Cost",     1.0);
    binCost_ = list.sublist("Problem").get("Integrality Cost", 1.0);
  }

  void update(const ROL::Vector<Real> &x, bool flag = true, int iter = -1) {
    misfit_->update(x,flag,iter);
    penalty_->update(x,flag,iter);
    binary_->update(x,flag,iter);
  }

  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    Real misfit  = misfit_->value(x,tol);
    Real penalty = penalty_->value(x,tol);
    Real binary  = binary_->value(x,tol);
    if (printToStream_) {
      *stream_ << "Unscaled: "
               << "Misfit Value = "  << misfit  << "  "
               << "Penalty Value = " << penalty << "  "
               << "Binary Value = "  << binary  << std::endl;
      *stream_ << "Scaled:   "
               << "Misfit Value = "  << misCost_*misfit  << "  "
               << "Penalty Value = " << penCost_*penalty << "  "
               << "Binary Value = "  << binCost_*binary  << std::endl;
    }
    return misCost_*misfit + penCost_*penalty + binCost_*binary;
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol ) {
    if (!initialized_) {
      xdual_ = g.clone();
      initialized_ = true;
    }
    g.zero();
    misfit_->gradient(g,x,tol);
    g.scale(misCost_);
    penalty_->gradient(*xdual_,x,tol);
    g.axpy(penCost_,*xdual_);
    binary_->gradient(*xdual_,x,tol);
    g.axpy(binCost_,*xdual_);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol ) {
    if (!initialized_) {
      xdual_ = hv.clone();
      initialized_ = true;
    }
    hv.zero();
    misfit_->hessVec(hv,v,x,tol);
    hv.scale(misCost_);
    penalty_->hessVec(*xdual_,v,x,tol);
    hv.axpy(penCost_,*xdual_);
    binary_->hessVec(*xdual_,v,x,tol);
    hv.axpy(binCost_,*xdual_);
  }

}; // class Sum_Objective

#endif
