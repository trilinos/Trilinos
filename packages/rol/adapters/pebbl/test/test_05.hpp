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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// ROL_Types contains predefined constants and objects
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_OptimizationSolver.hpp"
#include "ROL_Constraint_PEBBL.hpp"
#include "ROL_OptimizationProblemFactory.hpp"

template<class Real>
class Constraint_SimpleBinary : public ROL::Constraint<Real> {
private:
  const Real budget_;

  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  Constraint_SimpleBinary(int budget = 4) : budget_(budget) {}

  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &x, 
             Real &tol) {
    ROL::Ptr<std::vector<Real> >       cp = getVector(c);
    ROL::Ptr<const std::vector<Real> > xp = getConstVector(x);
    int dim(xp->size());
    (*cp)[0] = static_cast<Real>(0);
    for (int i = 0; i < dim; ++i) {
      (*cp)[0] += (*xp)[i];
    }
    (*cp)[0] -= static_cast<Real>(budget_);
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &x, 
                     Real &tol) {
    ROL::Ptr<std::vector<Real> >      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    (*jvp)[0] = static_cast<Real>(0);
    int dim(vp->size());
    for (int i = 0; i < dim; ++i) {
      (*jvp)[0] += (*vp)[i];
    }
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &x, 
                            Real &tol) {
    ROL::Ptr<std::vector<Real> >      jvp = getVector(ajv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    int dim(jvp->size());
    for (int i = 0; i < dim; ++i) {
      (*jvp)[i] = (*vp)[0];
    }
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                           Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class Objective_SimpleBinary : public ROL::Objective<Real> {
private:
  const std::vector<Real> alpha_, beta_;

  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  Objective_SimpleBinary(const std::vector<Real> &alpha, const std::vector<Real> &beta)
    : alpha_(alpha), beta_(beta) {}

  Real value(const ROL::Vector<Real> &x, Real &tol) {
    ROL::Ptr<const std::vector<Real> > xp = getConstVector(x);
    const Real half(0.5);
    int dim(xp->size());
    Real val(0);
    for (int i = 0; i < dim; ++i) {
      val += half * alpha_[i] * (*xp)[i] * (*xp)[i] - beta_[i] * (*xp)[i];
    }
    return val;
  }

  void gradient(ROL::Vector<Real> &g,
          const ROL::Vector<Real> &x,
                Real &tol) {
    ROL::Ptr<std::vector<Real> >       gp = getVector(g);
    ROL::Ptr<const std::vector<Real> > xp = getConstVector(x);
    int dim(xp->size());
    for (int i = 0; i < dim; ++i) {
      (*gp)[i] = alpha_[i] * (*xp)[i] - beta_[i];
    }
  }

  void hessVec(ROL::Vector<Real> &hv,
         const ROL::Vector<Real> &v,
         const ROL::Vector<Real> &x,
               Real &tol ) {
    ROL::Ptr<std::vector<Real> >      hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real> > vp = getConstVector(v);
    int dim(vp->size());
    for (int i = 0; i < dim; ++i) {
      (*hvp)[i] = alpha_[i] * (*vp)[i];
    }
  }
};

template<class Real>
class Test05Factory : public ROL::OptimizationProblemFactory<Real> {
private:
  std::vector<Real> alpha_, beta_;
  const int N_;
  ROL::ParameterList pl_;

public:
  Test05Factory(ROL::ParameterList &pl)
    : N_(10), pl_(pl) {
    alpha_.resize(N_); beta_.resize(N_);
    for (int i = 0; i < N_; ++i) {
      alpha_[i] = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
      beta_[i]  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    }
  }

  void update(void) {}

  ROL::Ptr<ROL::Objective<Real>> buildObjective(void) {
    return ROL::makePtr<Objective_SimpleBinary<Real>>(alpha_,beta_);
  }

  ROL::Ptr<ROL::Vector<Real>> buildSolutionVector(void) {
    return ROL::makePtr<ROL::StdVector<Real>>(N_,0.0);
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildBoundConstraint(void) {
    ROL::Ptr<ROL::Vector<Real>> xl = ROL::makePtr<ROL::StdVector<Real>>(N_,0.0);
    ROL::Ptr<ROL::Vector<Real>> xu = ROL::makePtr<ROL::StdVector<Real>>(N_,1.0);
    return ROL::makePtr<ROL::Bounds<Real>>(xl,xu);
  }

  ROL::Ptr<ROL::Constraint<Real>> buildEqualityConstraint(void) {
    bool useIneq = pl_.get("Use Inequality", true);
    if (!useIneq) {
      Real budget  = static_cast<Real>(pl_.get("Budget", 3));
      return ROL::makePtr<Constraint_SimpleBinary<Real>>(static_cast<int>(budget));
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Vector<Real>> buildEqualityMultiplier(void) {
    bool useIneq = pl_.get("Use Inequality", true);
    if (!useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Constraint<Real>> buildInequalityConstraint(void) {
    bool useIneq = pl_.get("Use Inequality", true);
    if (useIneq) {
      return ROL::makePtr<Constraint_SimpleBinary<Real>>(0);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::Vector<Real>> buildInequalityMultiplier(void) {
    bool useIneq = pl_.get("Use Inequality", true);
    if (useIneq) {
      return ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
    }
    return ROL::nullPtr;
  }

  ROL::Ptr<ROL::BoundConstraint<Real>> buildInequalityBoundConstraint(void) {
    bool useIneq = pl_.get("Use Inequality", true);
    if (useIneq) {
      Real budget  = static_cast<Real>(pl_.get("Budget", 3));
      ROL::Ptr<ROL::Vector<Real>> ilo = ROL::makePtr<ROL::StdVector<Real>>(1,0.0);
      ROL::Ptr<ROL::Vector<Real>> iup = ROL::makePtr<ROL::StdVector<Real>>(1,budget);
      return ROL::makePtr<ROL::Bounds<Real>>(ilo,iup);
    }
    return ROL::nullPtr;
  }
};
