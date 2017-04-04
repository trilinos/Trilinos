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

#ifndef ROL_BPOEOBJECTIVE_HPP
#define ROL_BPOEOBJECTIVE_HPP

#include "ROL_RiskAverseObjective.hpp"
#include "ROL_BPOE.hpp"

namespace ROL {

template<class Real>
class BPOEObjective : public Objective<Real> {
private:
  Teuchos::RCP<RiskMeasure<Real> > bpoe_;
  Teuchos::RCP<Objective<Real> > riskObj_;

public:
  BPOEObjective( const Teuchos::RCP<Objective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> > &vsampler, 
                 const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                 const Teuchos::RCP<SampleGenerator<Real> > &hsampler,
                 const bool storage = true ) {
    bpoe_    = Teuchos::rcp(new BPOE<Real>(threshold,order));
    riskObj_ = Teuchos::rcp(new RiskAverseObjective<Real>(pObj,bpoe_,vsampler,gsampler,hsampler,storage));
  }

  BPOEObjective( const Teuchos::RCP<Objective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> > &vsampler, 
                 const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                 const bool storage = true ) {
    bpoe_    = Teuchos::rcp(new BPOE<Real>(threshold,order));
    riskObj_ = Teuchos::rcp(new RiskAverseObjective<Real>(pObj,bpoe_,vsampler,gsampler,storage));
  }

  BPOEObjective( const Teuchos::RCP<Objective<Real> > &pObj,
                 const Real order, const Real threshold,
                 const Teuchos::RCP<SampleGenerator<Real> > &sampler,
                 const bool storage = true ) {
    bpoe_    = Teuchos::rcp(new BPOE<Real>(threshold,order));
    riskObj_ = Teuchos::rcp(new RiskAverseObjective<Real>(pObj,bpoe_,sampler,storage));
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    riskObj_->update(x,flag,iter);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    return riskObj_->value(x,tol);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    riskObj_->gradient(g,x,tol);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    riskObj_->hessVec(hv,v,x,tol);
  }

  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v,
                        const Vector<Real> &x, Real &tol ) {
    riskObj_->precond(Pv,v,x,tol);
  }
};

}

#endif
