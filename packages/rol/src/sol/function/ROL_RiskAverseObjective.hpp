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

#ifndef ROL_RISKAVERSEOBJECTIVE_HPP
#define ROL_RISKAVERSEOBJECTIVE_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_ParametrizedObjective.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_RiskMeasure.hpp"

namespace ROL {

template<class Real>
class RiskAverseObjective : public Objective<Real> {
private:
  // Problem Data
  Teuchos::RCP<ParametrizedObjective<Real> > pObj_;
  Teuchos::RCP<SampleGenerator<Real> > vsampler_;
  Teuchos::RCP<SampleGenerator<Real> > gsampler_;
  Teuchos::RCP<RiskMeasure<Real> > rm_;

  // Storage Information
  bool storage_;
  std::map<std::vector<Real>,Real> value_storage_;
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > gradient_storage_;

public:
  virtual ~RiskAverseObjective() {}

  RiskAverseObjective( ParametrizedObjective<Real> &pObj, RiskMeasure<Real> &rm,  
                       SampleGenerator<Real> &vsampler, SampleGenerator<Real> &gsampler,
                       bool storage = true ) : storage_(storage) { 
    pObj_     = Teuchos::rcp(&pObj,false);     // Parametrized Objective Function Object
    vsampler_ = Teuchos::rcp(&vsampler,false); // Objective Function Value Sampler Object
    gsampler_ = Teuchos::rcp(&gsampler,false); // Gradient Sampler Object
    rm_       = Teuchos::rcp(&rm,false);       // Risk Measure Object
    
    value_storage_.clear();
    gradient_storage_.clear();
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    pObj_->update(x,flag,iter);
    vsampler_->update(x);
    if ( storage_ ) {
      value_storage_.clear();
    }
    if ( flag ) {
      gsampler_->update(x);
      if ( storage_ ) {
        gradient_storage_.clear();
      }
    }
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
    Real val = 0.0;
    Teuchos::RCP<Vector<Real> > x0;
    rm_->reset(x0,x);
    for ( int i = 0; i < vsampler_->numMySamples(); i++ ) {
      pObj_->setParameter(vsampler_->getMyPoint(i));
      if ( storage_ && value_storage_.count(vsampler_->getMyPoint(i)) ) {
        val = value_storage_[vsampler_->getMyPoint(i)];
      }
      else {
        val = pObj_->value(*x0,tol);
        if ( storage_ ) {
          value_storage_.insert(std::pair<std::vector<Real>,Real>(vsampler_->getMyPoint(i),val));
        }
      }
      rm_->update(val,vsampler_->getMyWeight(i));
    }
    return rm_->getValue(*vsampler_);
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    //gsampler_->refine(x,tol);
    g.zero();
    Teuchos::RCP<Vector<Real> > x0;
    rm_->reset(x0,x);
    Teuchos::RCP<Vector<Real> > g0 = x0->clone();
    Real val = 0.0;
    for ( int i = 0; i < gsampler_->numMySamples(); i++ ) {
      pObj_->setParameter(gsampler_->getMyPoint(i));
      if ( storage_ && value_storage_.count(gsampler_->getMyPoint(i)) ) {
        val = value_storage_[gsampler_->getMyPoint(i)];
      }
      else {
        val = pObj_->value(*x0,tol);
        if ( storage_ ) {
          value_storage_.insert(std::pair<std::vector<Real>,Real>(gsampler_->getMyPoint(i),val));
        }
      }
      if ( storage_ && gradient_storage_.count(gsampler_->getMyPoint(i)) ) {
        g0->set(*(gradient_storage_[gsampler_->getMyPoint(i)]));
      }
      else { 
        pObj_->gradient(*g0,*x0,tol);
        if ( storage_ ) {
          Teuchos::RCP<Vector<Real> > tmp = g0->clone();
          gradient_storage_.insert(std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(gsampler_->getMyPoint(i),tmp));
          gradient_storage_[gsampler_->getMyPoint(i)]->set(*g0);
        }
      }
      rm_->update(val,*g0,gsampler_->getMyWeight(i));
    }
    rm_->getGradient(g,*gsampler_);
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, 
                        const Vector<Real> &x, Real &tol ) {
    hv.zero();
    Teuchos::RCP<Vector<Real> > x0;
    Teuchos::RCP<Vector<Real> > v0;
    rm_->reset(x0,x,v0,v);
    Real val = 0.0;
    Real gv  = 0.0;
    Teuchos::RCP<Vector<Real> > g0 = x0->clone();
    Teuchos::RCP<Vector<Real> > h0 = x0->clone();
    for ( int i = 0; i < gsampler_->numMySamples(); i++ ) {
      pObj_->setParameter(gsampler_->getMyPoint(i));
      if ( storage_ && value_storage_.count(gsampler_->getMyPoint(i)) ) {
        val = value_storage_[gsampler_->getMyPoint(i)];
      }
      else {
        val = pObj_->value(*x0,tol);
        if ( storage_ ) {
          value_storage_.insert(std::pair<std::vector<Real>,Real>(gsampler_->getMyPoint(i),val));
        }
      }
      if ( storage_ && gradient_storage_.count(gsampler_->getMyPoint(i)) ) {
        g0->set(*(gradient_storage_[gsampler_->getMyPoint(i)]));
      }
      else { 
        pObj_->gradient(*g0,*x0,tol);
        if ( storage_ ) {
          Teuchos::RCP<Vector<Real> > tmp = g0->clone();
          gradient_storage_.insert(std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(gsampler_->getMyPoint(i),tmp));
          gradient_storage_[gsampler_->getMyPoint(i)]->set(*g0);
        }
      }
      gv  = g0->dot(*v0);
      pObj_->hessVec(*h0,*v0,*x0,tol);
      rm_->update(val,*g0,gv,*h0,gsampler_->getMyWeight(i));
    }
    rm_->getHessVec(hv,*gsampler_);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }
};

}

#endif
