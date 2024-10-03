// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef PDEOPT_QOI_GINZBURGLANDAU_EX02_HPP
#define PDEOPT_QOI_GINZBURGLANDAU_EX02_HPP

#include "obj_ginzburg-landau.hpp"

template <class Real>
class QoI_GinzburgLandau_StateTracking_ex02 : public QoI_GinzburgLandau_StateTracking<Real> {
public:
  QoI_GinzburgLandau_StateTracking_ex02(const ROL::Ptr<FE<Real> > &fe,
                                        const ROL::Ptr<FieldHelper<Real> > &fieldHelper,
                                        Teuchos::ParameterList &parlist)
    : QoI_GinzburgLandau_StateTracking<Real>(fe,fieldHelper,parlist) {
    QoI_GinzburgLandau_StateTracking<Real>::computeTarget();
  }

  virtual Real evaluateRealTarget(const std::vector<Real> &x) const override {
    return static_cast<Real>(1)/std::sqrt(static_cast<Real>(2));
  }

  virtual Real evaluateImagTarget(const std::vector<Real> &x) const override {
    return static_cast<Real>(1)/std::sqrt(static_cast<Real>(2));
  }
}; // QoI_GinzburgLandau_StateTracking_ex02

#endif
