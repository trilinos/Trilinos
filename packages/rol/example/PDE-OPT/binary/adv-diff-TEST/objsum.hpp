// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  void update(const ROL::Vector<Real> &x, ROL::UpdateType type, int iter = -1) {
    misfit_->update(x,type,iter);
    penalty_->update(x,type,iter);
    binary_->update(x,type,iter);
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
