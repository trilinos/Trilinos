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
