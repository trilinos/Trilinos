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

#ifndef ROL_INTERIORPOINT_H
#define ROL_INTERIORPOINT_H

namespace ROL {


// Treat slack variable as Simulation variable
template <class Real>
class InteriorPointObjective : public Objective<Real> {
private:

  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type; 

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  Real mu_;
  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<Objective<Real> > barrier_;

public:
  InteriorPointObjective( Objective<Real> &obj, Objective<Real> &barrier, Real mu ) :
    mu_(mu) {

    obj_     = Teuchos::rcp(&obj, false);
    barrier_ = Teuchos::rcp(&barrier,false);
  }
 
  void updatePenalty( Real mu ) {
    mu_ = mu;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    
    PV xpv = Teuchos::dyn_cast<const PV>(x); 

    Teuchos::RCP<V> xopt  = xpv.get(OPT);    
    Teuchos::RCP<V> slack = xpv.get(SLACK);

    return obj_->value(*xopt,tol) + mu*barrier_(*slack,tol); 

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    PV gpv = Teuchos::dyn_cast<PV>(g); 
        
    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<V> gopt = gpv.get(OPT);
    Teuchos::RCP<V> gs   = gpv.get(SLACK);

    obj_->gradient(*gopt, *xopt, tol);
    barrier_->gradient(*gs, *s, tol); 
    gs->scale(mu)_;
   
  } 
   
  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                 const Vector<Real> &x, Real &tol ) {

    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    const PV vpv = Teuchos::dyn_cast<const PV>(v);
    PV hvpv = Teuchos::dyn_cast<PV>(hv); 
        
    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<const V> vopt = vpv.get(OPT);    
    Teuchos::RCP<const V> vs   = vpv.get(SLACK);

    Teuchos::RCP<V> hvopt = hvpv.get(OPT);
    Teuchos::RCP<V> hvs   = hvpv.get(SLACK);

    obj_->hessVec(*hvopt, *vopt, *xopt, tol);
    barrier_->hessVec(*hvs, *vs, *xs, tol);
    hvs->scale(mu_);
   
  }

}; // InteriorPointObjective




}

#endif
