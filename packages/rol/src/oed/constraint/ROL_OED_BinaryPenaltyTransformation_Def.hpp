// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BINARYPENALTYTRANSFORMATION_DEF_HPP
#define ROL_OED_BINARYPENALTYTRANSFORMATION_DEF_HPP

namespace ROL::OED {

template<typename Real>
BinaryPenaltyTransformation<Real>::BinaryPenaltyTransformation(ParameterList& parlist)
  : ProfiledClass<Real,std::string>("OED::BinaryPenaltyTransformation"),
    pmin_(parlist.sublist("OED").sublist("Binary Penalty Transformation").get("Minimum Value",0.0)),
    ppow_(parlist.sublist("OED").sublist("Binary Penalty Transformation").get("Strength Value",3.0)),
    type_(parlist.sublist("OED").sublist("Binary Penalty Transformation").get("Type",0)) {
  pmin_ = std::max(pmin_,static_cast<Real>(0));
  ppow_ = std::max(ppow_,static_cast<Real>(2));
}

template<typename Real>
void BinaryPenaltyTransformation<Real>::value(
       std::vector<Real> &c,const std::vector<Real> &x,Real &tol) {
  startTimer("value");
  const Real one(1);
  const unsigned dim=x.size();
  Real val(0);
  for(unsigned i=0u; i<dim; ++i) {
    if (type_==1)      val = x[i]/(one+ppow_*(one-x[i]));
    else if (type_==2) val = (std::exp(ppow_*x[i])-one)/(std::exp(ppow_)-one);
    else               val = std::pow(x[i],ppow_);
    c[i] = pmin_ + (one-pmin_)*val;
  }
  stopTimer("value");
}

template<typename Real>
void BinaryPenaltyTransformation<Real>::applyJacobian(
       std::vector<Real> &jv,const std::vector<Real> &v,
       const std::vector<Real> &x,Real &tol) {
  startTimer("applyJacobian");
  const Real one(1), two(2);
  const unsigned dim=x.size();
  Real val(0);
  for(unsigned i=0u; i<dim; ++i) {
    if (type_==1)      val = (one+ppow_)/std::pow(one+ppow_*(one-x[i]),two);
    else if (type_==2) val = ppow_*std::exp(ppow_*x[i])/(std::exp(ppow_)-one);
    else               val = ppow_*std::pow(x[i],ppow_-one);
    jv[i] = (one-pmin_)*val*v[i];
  }
  stopTimer("applyJacobian");
}

template<typename Real>
void BinaryPenaltyTransformation<Real>::applyAdjointJacobian(
       std::vector<Real> &ajv,const std::vector<Real> &v,
       const std::vector<Real> &x,Real &tol) {
  startTimer("applyAdjointJacobian");
  const Real one(1), two(2);
  const unsigned dim=x.size();
  Real val(0);
  for(unsigned i=0u; i<dim; ++i) {
    if (type_==1)      val = (one+ppow_)/std::pow(one+ppow_*(one-x[i]),two);
    else if (type_==2) val = ppow_*std::exp(ppow_*x[i])/(std::exp(ppow_)-one);
    else               val = ppow_*std::pow(x[i],ppow_-one);
    ajv[i] = (one-pmin_)*val*v[i];
  }
  stopTimer("applyAdjointJacobian");
}

template<typename Real>
void BinaryPenaltyTransformation<Real>::applyAdjointHessian(
       std::vector<Real> &ahwv,const std::vector<Real> &w,
       const std::vector<Real> &v,const std::vector<Real> &x,Real &tol) {
  startTimer("applyAdjointHessian");
  const Real one(1), two(2), three(3);
  const unsigned dim=x.size();
  Real val(0);
  for(unsigned i=0u; i<dim; ++i) {
    if (type_==1)      val = two*ppow_*(one+ppow_)/std::pow(one+ppow_*(one-x[i]),three);
    else if (type_==2) val = ppow_*ppow_*std::exp(ppow_*x[i])/(std::exp(ppow_)-one);
    else               val = (ppow_-one)*ppow_*std::pow(x[i],ppow_-two);
    ahwv[i] = (one-pmin_)*val*w[i]*v[i];
  }
  stopTimer("applyAdjointHessian");
}

} // End ROL::OED Namespace

#endif
