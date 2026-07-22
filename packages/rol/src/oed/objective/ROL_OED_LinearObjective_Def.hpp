// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_LINEAR_OBJECTIVE_DEF_HPP
#define ROL_OED_LINEAR_OBJECTIVE_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void LinearObjective<Real>::computeG(Vector<Real> &g) {
  startTimer("computeG");
  if (type_ != "C" && type_ != "D") {
    std::vector<Real> param = Objective_SimOpt<Real>::getParameter();
    if (type_ == "I" || type_ == "R") {
      factors_->evaluate(g,param);
    }
    else if (type_ == "A") {
      traceSampler_->get(g,param);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OED::LinearObjective::computeG : Optimality type not implemented!");
    }
  }
  stopTimer("computeG");
}

// D Optimality
template<typename Real>
LinearObjective<Real>::LinearObjective() : ProfiledClass<Real,std::string>("OED::LinearObjective"), type_("D") {}

// I and R Optimality
template<typename Real>
LinearObjective<Real>::LinearObjective(const Ptr<Factors<Real>> &factors,
                const std::string             &type)
  : ProfiledClass<Real,std::string>("OED::LinearObjective"),
    factors_(factors), type_(type), g_(factors_->get(0)->clone()) {
  if (type_ != "I" && type_ != "R") {
    std::stringstream ss;
    ss << ">>> ROL::OED::LinearObjective : Wrong constructor for " << type_ << "-optimality!";
    throw Exception::NotImplemented(ss.str());
  }
}

// C optimality
template<typename Real>
LinearObjective<Real>::LinearObjective(const Ptr<Vector<Real>> &c)
  : ProfiledClass<Real,std::string>("OED::LinearObjective"),
    type_("C"), g_(c->clone()) {
  g_->set(*c);
}

// A Optimality
template<typename Real>
LinearObjective<Real>::LinearObjective(const Ptr<Vector<Real>>  &theta,
                const Ptr<TraceSampler<Real>> &traceSampler)
  : ProfiledClass<Real,std::string>("OED::LinearObjective"),
    traceSampler_(traceSampler), type_("A"), g_(theta->dual().clone()) {}

template<typename Real>
Real LinearObjective<Real>::value( const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("value");
  computeG(*g_);
  Real val = g_->apply(u);
  stopTimer("value");
  return val;
}

template<typename Real>
void LinearObjective<Real>::gradient_1( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("gradient_1");
  computeG(*g_);
  g.set(*g_);
  stopTimer("gradient_1");
}

template<typename Real>
void LinearObjective<Real>::gradient_2( Vector<Real> &g, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("gradient_2");
  g.zero();
  stopTimer("gradient_2");
}

template<typename Real>
void LinearObjective<Real>::hessVec_11( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_11");
  hv.zero();
  stopTimer("hessVec_11");
}

template<typename Real>
void LinearObjective<Real>::hessVec_12( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_12");
  hv.zero();
  stopTimer("hessVec_12");
}

template<typename Real>
void LinearObjective<Real>::hessVec_21( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_21");
  hv.zero();
  stopTimer("hessVec_21");
}

template<typename Real>
void LinearObjective<Real>::hessVec_22( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &u, const Vector<Real> &z, Real &tol ) {
  startTimer("hessVec_22");
  hv.zero();
  stopTimer("hessVec_22");
}

} // End OED Namespace
} // End ROL Namespace

#endif
