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
//
#ifndef PH_OBJECTIVE_H
#define PH_OBJECTIVE_H

#include "ROL_PH_RiskObjective.hpp"
#include "ROL_PH_DeviationObjective.hpp"
#include "ROL_PH_RegretObjective.hpp"
#include "ROL_PH_ErrorObjective.hpp"
#include "ROL_PH_ProbObjective.hpp"
#include "ROL_PH_bPOEObjective.hpp"

/** @ingroup func_group
    \class ROL::PH_Objective
    \brief Provides the interface for the progressive hedging objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_Objective : public Objective<Real> {
private:
  Ptr<Objective<Real>> obj_;
  Ptr<const Vector<Real>> xbar_, w_;
  Ptr<Vector<Real>> xprimal_;
  Real penaltyParam_;

public:

  PH_Objective(const Ptr<Objective<Real>> &obj,
               const Ptr<Vector<Real>> &x,
               const Real penaltyParam,
               ParameterList &parlist) 
    : xbar_(nullPtr), w_(nullPtr), xprimal_(nullPtr),
      penaltyParam_(penaltyParam) {
    xprimal_ = x->clone();
    std::string type = parlist.sublist("SOL").get("Stochastic Component Type","Risk Neutral");
    if (type == "Risk Averse") {
      obj_ = makePtr<PH_RiskObjective<Real>>(obj,parlist);
    }
    else if (type == "Deviation") {
      obj_ = makePtr<PH_DeviationObjective<Real>>(obj,parlist);
    }
    else if (type == "Regret") {
      obj_ = makePtr<PH_RegretObjective<Real>>(obj,parlist);
    }
    else if (type == "Error") {
      obj_ = makePtr<PH_ErrorObjective<Real>>(obj,parlist);
    }
    else if (type == "Probability") {
      std::string prob = parlist.sublist("SOL").sublist("Probability").get("Name","bPOE");
      if (prob == "Smoothed POE") {
        obj_ = makePtr<PH_ProbObjective<Real>>(obj,parlist);
      }
      else if (prob == "bPOE") {
        obj_ = makePtr<PH_bPOEObjective<Real>>(obj,parlist);
      }
      else {
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               "Invalid probability type " << prob << "!");
      }
    }
    else if (type == "Risk Neutral") {
      obj_ = obj;
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                             "Invalid stochastic component type " << type << "!");
    }
  }

  void setData(const Ptr<const Vector<Real>> &xbar,
               const Ptr<const Vector<Real>> &w,
               const Real penaltyParam) {
    xbar_         = xbar;
    w_            = w;
    penaltyParam_ = penaltyParam;
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const Real half(0.5), one(1);
    Real val  = obj_->value(x,tol);
    Real wx   = x.dot(*w_);
    xprimal_->set(x);
    xprimal_->axpy(-one,*xbar_);
    Real xx   = xprimal_->dot(*xprimal_);
    return val + wx + half*penaltyParam_*xx;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    obj_->gradient(g,x,tol);
    xprimal_->set(*w_);
    xprimal_->axpy(penaltyParam_,x);
    xprimal_->axpy(-penaltyParam_,*xbar_);
    g.plus(xprimal_->dual());
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    obj_->hessVec(hv,v,x,tol);
    hv.axpy(penaltyParam_,v.dual());
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
