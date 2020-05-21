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

#include "test_01.hpp"

// ROL_Types contains predefined constants and objects
#include "ROL_PEBBL_IntegerProblemFactory.hpp"

template<class Real>
class Test05Factory : public ROL::PEBBL::IntegerProblemFactory<Real> {
private:
  std::vector<Real> alpha_, beta_;
  const int N_;
  bool useIneq_;
  Real budget_;
  ROL::Ptr<ROL::Vector<Real>> xl_, xu_, ilo_, iup_;

public:
  Test05Factory(ROL::ParameterList &pl) : N_(10) {
    alpha_.resize(N_); beta_.resize(N_);
    for (int i = 0; i < N_; ++i) {
      alpha_[i] = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      beta_[i]  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    }
    std::cout << std::endl;
    useIneq_ = pl.get("Use Inequality", true);
    budget_  = static_cast<Real>(pl.get("Budget", 3));
    xl_  = ROL::makePtr<ROL::StdVector<Real>>(N_,0.0);
    xu_  = ROL::makePtr<ROL::StdVector<Real>>(N_,1.0);
    ilo_ = ROL::makePtr<ROL::StdVector<Real>>(1,-budget_);
    iup_ = ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
  }

  ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>> build(void) {
    ROL::Ptr<ROL::Objective<Real>>        obj = ROL::makePtr<Objective_SimpleBinary<Real>>(alpha_,beta_);
    ROL::Ptr<ROL::Vector<Real>>             x = ROL::makePtr<ROL::StdVector<Real>>(N_,0.0);
    ROL::Ptr<ROL::BoundConstraint<Real>>  bnd = ROL::makePtr<ROL::Bounds<Real>>(xl_,xu_);
    ROL::Ptr<ROL::Constraint<Real>>      icon = ROL::makePtr<Constraint_SimpleBinary<Real>>(static_cast<int>(budget_));
    ROL::Ptr<ROL::Vector<Real>>          imul = ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
    ROL::Ptr<ROL::BoundConstraint<Real>> ibnd = ROL::makePtr<ROL::Bounds<Real>>(ilo_,iup_);
    ROL::Ptr<ROL::PEBBL::IntegerProblem<Real>>
      problem = ROL::makePtr<ROL::PEBBL::IntegerProblem<Real>>(obj,x);
    problem->addBoundConstraint(bnd);
    if (useIneq_) problem->addLinearConstraint("Linear",icon,imul,ibnd);
    else          problem->addLinearConstraint("Linear",icon,imul);
    return problem;
  }
};
