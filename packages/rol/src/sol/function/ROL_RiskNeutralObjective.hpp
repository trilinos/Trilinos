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

#ifndef ROL_RISKNEUTRALOBJECTIVE_HPP
#define ROL_RISKNEUTRALOBJECTIVE_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_ParametrizedObjective.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class RiskNeutralObjective : public Objective<Real> {
private:
  Teuchos::RCP<ParametrizedObjective<Real> > pObj_;
  Teuchos::RCP<SampleGenerator<Real> > vsampler_;
  Teuchos::RCP<SampleGenerator<Real> > gsampler_;

  Real value_;
  Teuchos::RCP<Vector<Real> > gradient_;
 
  bool firstUpdate_;

public:
  virtual ~RiskNeutralObjective() {}

  RiskNeutralObjective( ParametrizedObjective<Real> &pObj, 
                        SampleGenerator<Real> &vsampler, 
                        SampleGenerator<Real> &gsampler ) {
    pObj_     = Teuchos::rcp(&pObj,false);     // Parametrized Objective Function Object
    vsampler_ = Teuchos::rcp(&vsampler,false); // Objective Function Value Sampler Object
    gsampler_ = Teuchos::rcp(&gsampler,false); // Gradient Sampler Object
    firstUpdate_ = true;
  }

  // Delegating constructors require C++11.
  //RiskNeutralObjective( ParametrizedObjective<Real> &pObj, SampleGenerator<Real> &sampler ) 
  //  : RiskNeutralObjective(pObj,sampler,sampler) {}
  RiskNeutralObjective( ParametrizedObjective<Real> &pObj, SampleGenerator<Real> &sampler ) {
    pObj_     = Teuchos::rcp(&pObj,false);     // Parametrized Objective Function Object
    vsampler_ = Teuchos::rcp(&sampler,false); // Objective Function Value Sampler Object
    gsampler_ = Teuchos::rcp(&sampler,false); // Gradient Sampler Object
    firstUpdate_ = true;
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( firstUpdate_ ) {
      gradient_ = (x.dual()).clone();
    }
    pObj_->update(x,flag,iter);
    vsampler_->update(x);
    value_ = 0.0;
    if ( flag ) {
      gsampler_->update(x);
      gradient_->zero();
    }
  }

  virtual Real value( const Vector<Real> &x, Real &tol ) {
//    std::cout << " Initial value = " << value_ << "\n";
    Real myval = 0.0, ptval = 0.0, val = 0.0, error = 2.0*tol + 1.0;
    std::vector<Real> ptvals;
    while ( error > tol ) {
      vsampler_->refine();
      for ( int i = vsampler_->start(); i < vsampler_->numMySamples(); i++ ) {
        pObj_->setParameter(vsampler_->getMyPoint(i));
        ptval  = pObj_->value(x,tol);
        myval += vsampler_->getMyWeight(i)*ptval;
        ptvals.push_back(ptval);
//        std::cout << "       ptval-" << i << " = " << ptval 
//                  << "  weight-"     << i << " = " << vsampler_->getMyWeight(i) << "\n";
      }
      error = vsampler_->computeError(ptvals);
//      std::cout << " v error = "     << error 
//                << "  tol = "        << tol 
//                << "  myval = "      << myval
//                << "  num points = " << vsampler_->numMySamples() 
//                << "  start = "      << vsampler_->start() << "\n";
      ptvals.clear();
    }
    vsampler_->sumAll(&myval,&val,1);
    value_ += val;
//    std::cout << " Final value = " << value_ << "\n";
    vsampler_->setSamples();
    return value_;
  }

  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
//    std::cout << " Initial norm(gradient) = " << gradient_->norm() << "\n";
    g.zero();
    Teuchos::RCP<Vector<Real> > ptg = g.clone(); ptg->zero();
    Teuchos::RCP<Vector<Real> > myg = g.clone(); myg->zero();
    std::vector<Teuchos::RCP<Vector<Real> > > ptgs;
    Real error = 2.0*tol + 1.0;
    while ( error > tol ) {
      gsampler_->refine();
      for ( unsigned i = gsampler_->start(); i < gsampler_->numMySamples(); i++ ) {
        pObj_->setParameter(gsampler_->getMyPoint(i));
        pObj_->gradient(*ptg,x,tol);
        myg->axpy(gsampler_->getMyWeight(i),*ptg); 
        ptgs.push_back(x.clone());
        (ptgs.back())->set(*ptg);
      }
      error = gsampler_->computeError(ptgs,x);
//      std::cout << " g error = "     << error 
//                << "  tol = "        << tol 
//                << "  norm(g) = "    << myg->norm() 
//                << "  num points = " << gsampler_->numMySamples() 
//                << "  start = "      << gsampler_->start() 
//                << "\n";
      ptgs.clear();
    }
    gsampler_->sumAll(*myg,g);
    gradient_->axpy(1.0,g);
    g.set(*(gradient_));
//    std::cout << " Final norm(gradient) = " << gradient_->norm() << "\n";
    gsampler_->setSamples();
  }

  virtual void hessVec( Vector<Real> &hv, const Vector<Real> &v, 
                        const Vector<Real> &x, Real &tol ) {
    hv.zero();
    Teuchos::RCP<Vector<Real> > pth = hv.clone(); pth->zero();
    Teuchos::RCP<Vector<Real> > myh = hv.clone(); myh->zero();
    for ( int i = 0; i < gsampler_->numMySamples(); i++ ) {
      pObj_->setParameter(gsampler_->getMyPoint(i));
      pObj_->hessVec(*pth,v,x,tol);
      myh->axpy(gsampler_->getMyWeight(i),*pth);
    }
    gsampler_->sumAll(*myh,hv);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x ) {
    Pv.set(v.dual());
  }
};

}

#endif

//  virtual Real value( const Vector<Real> &x, Real &tol ) {
//    Real myval = 0.0, ptval = 0.0, val = 0.0;
//    for ( unsigned i = 0; i < vsampler_->numMySamples(); i++ ) {
//      pObj_->setParameter(vsampler_->getMyPoint(i));
//      ptval  = pObj_->value(x,tol);
//      myval += vsampler_->getMyWeight(i)*ptval;
//    }
//    vsampler_->sumAll(&myval,&val,1);
//    return val;
//  }

//  virtual void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
//    g.zero();
//    Teuchos::RCP<Vector<Real> > ptg = x.clone(); ptg->zero();
//    Teuchos::RCP<Vector<Real> > myg = x.clone(); myg->zero();
//    for ( unsigned i = 0; i < gsampler_->numMySamples(); i++ ) {
//      pObj_->setParameter(gsampler_->getMyPoint(i));
//      pObj_->gradient(*ptg,x,tol);
//      myg->axpy(gsampler_->getMyWeight(i),*ptg); 
//    }
//    gsampler_->sumAll(*myg,g);
//  }
