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

/*! \file  obj.hpp
    \brief Provides the interface for local (cell-based) objective function computations.
*/

#ifndef BINARY_CON_STEFAN_BOLTZMANN_HPP
#define BINARY_CON_STEFAN_BOLTZMANN_HPP

#include "ROL_StdConstraint.hpp"
#include "ROL_Bounds.hpp"

template <class Real>
class BudgetConstraint : public ROL::StdConstraint<Real> {
private:
  Real budget_;

public:
  BudgetConstraint(ROL::ParameterList &pl) {
    budget_ = pl.sublist("Problem").get("Control Budget",8.0);
  }

  void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) override {
    c[0] = -budget_;
    for (const auto xv : x) c[0] += xv;
  }

  void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                      const std::vector<Real> &x, Real &tol ) override {
    jv[0] = static_cast<Real>(0);
    for (const auto vv : v) jv[0] += vv;
  }
  
  void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                             const std::vector<Real> &x, Real &tol ) override {
    ajv.assign(ajv.size(),v[0]);
  }

  void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                            const std::vector<Real> &v, const std::vector<Real> &x,
                            Real &tol ) {
    ahuv.assign(ahuv.size(),static_cast<Real>(0));
  }

  ROL::Ptr<ROL::Vector<Real>> createMultiplier(void) const {
    return ROL::makePtr<ROL::StdVector<Real>>(1,static_cast<Real>(0));
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> createBounds(void) const {
    ROL::Ptr<ROL::Vector<Real>> l = createMultiplier(); l->setScalar(-budget_);
    ROL::Ptr<ROL::Vector<Real>> u = createMultiplier(); u->zero();
    return ROL::makePtr<ROL::Bounds<Real>>(l,u);
  }
}; // BudgetConstraint

#endif
