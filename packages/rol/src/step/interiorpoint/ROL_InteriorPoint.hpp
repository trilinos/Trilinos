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


/** @ingroup func_group
 *  \class ROL::InteriorPointObjective
 *  \brief Adds barrier term to generic objective
 */
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
  InteriorPointObjective( const Teuchos::RCP<Objective<Real> > &obj, 
                          const Teuchos::RCP<Objective<Real> > &barrier, Real mu ) :
    mu_(mu) {

    obj_     = Teuchos::rcp(&obj, false);
    barrier_ = Teuchos::rcp(&barrier, false);
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

}; // class InteriorPointObjective


/** @ingroup func_group
 *  \class ROL::InteriorPointEqualityConstraint
 *  \brief Treat inequality constraint as equality with slack variable
 */

template<class Real> 
class InteriorPointEqualityConstraint : public EqualityConstraint<Real> {
private:
  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type; 

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  const static size_type EQUAL = 0;
  const static size_type INEQ  = 1;

  Teuchos::RCP<EqualityConstraint<Real> > eqcon_;
  Teuchos::RCP<EqualityConstraint<Real> > incon_;
   
public:

  InteriorPointEqualityConstraint( const Teuchos::RCP<EqualityConstraint<Real> > &eqcon,
                                   const Teuchos::RCP<EqualityConstraint<Real> > &incon ) :
    eqcon_(eqcon), incon_(incon) {}

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    PV cpv = Teuchos::dyn_cast<PV>(c); 

    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<V> ce = cpv.get(EQUAL);
    Teuchos::RCP<V> ci = cpv.get(INEQ);

    eqcon_->value(*ce, *xopt, tol);
    incon_->value(*ci, *xopt, tol);
    ci->axpy(-1.0,*s);
  }

  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) {
    
    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    const PV vpv = Teuchos::dyn_cast<const PV>(v);
    PV jvpv = Teuchos::dyn_cast<PV>(jv); 

    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<const V> vopt = vpv.get(OPT);    
    Teuchos::RCP<const V> vs   = vpv.get(SLACK);

    Teuchos::RCP<V> jvopt = jvpv.get(OPT);
    Teuchos::RCP<V> jvs   = jvpv.get(SLACK);

    eqcon_->applyJacobian(*jvopt, *vopt, *xopt, tol);  
    incon_->applyJacobian(*jvs, *vopt, *xopt, tol);
    jvs->axpy(-1.0,vopt);
  }

  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) {

    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    const PV vpv = Teuchos::dyn_cast<const PV>(v);
    PV ajvpv = Teuchos::dyn_cast<PV>(ajv); 

    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<const V> ve = vpv.get(EQUAL);    
    Teuchos::RCP<const V> vi = vpv.get(INEQ);

    Teuchos::RCP<V> ajve = ajvpv.get(EQUAL);
    Teuchos::RCP<V> ajvi = ajvpv.get(INEQ);

    eqcon_->applyAdjointJacobian(*ajve,*ve,*xopt,tol);
    incon_->applyAdjointJacobian(*ajvi,*vi,*xopt,tol);
    ajve->plus(*ajvi);
    ajvi->set(*vi);
    ajvi->scale(-1.0);
  }

  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) {
  
    const PV xpv = Teuchos::dyn_cast<const PV>(x);
    const PV upv = Teuchos::dyn_cast<const PV>(u);
    const PV vpv = Teuchos::dyn_cast<const PV>(v);
    PV ahuvpv = Teuchos::dyn_cast<PV>(ahuv); 

    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    Teuchos::RCP<const V> vopt = vpv.get(OPT);    
    
    Teuchos::RCP<const V> ue   = upv.get(EQUAL);    
    Teuchos::RCP<const V> ui   = upv.get(INEQ);
    
    Teuchos::RCP<V> ahuvopt = ahuvpv.get(OPT);
    Teuchos::RCP<V> ahuvs   = ahuvpv.get(SLACK);

    Teuchos::RCP<V> temp = ahuvopt->clone();
    
    eqcon_->applyAdjointHessian(*ahuvopt,*ue,*vopt,*xopt,tol);
    incon_->applyAdjointHessian(*temp,*ui,*vopt,*xopt,tol);

    ahuvopt->plus(temp);
    ahuvs->zero();

  }
   
}; // class InteriorPointEqualityConstraint


/** @ingroup func_group
 *  \class ROL::InteriorPointBoundConstraint
 *  \brief Require positivity of slack variables
 */

template<class Real>
class InteriorPointBoundConstraint : public BoundConstraint<Real> {
private:
  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type; 

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  Teuchos::RCP<BoundConstraint<Real> > bc_;

  Teuchos::RCP<Vector<Real> > lower_;
  Teuchos::RCP<Vector<Real> > upper_;  

public:

  InteriorPointBoundConstraint( const Vector<Real> &x ) {

    lower_ = x.clone();
    upper_ = x.clone();

    PV lowerpv = Teuchos::dyn_cast<PV>(lower_);
    PV upperpv = Teuchos::dyn_cast<PV>(upper_);
   
    Teuchos::RCP<V> lower_opt   = lowerpv.get(OPT);
    Teuchos::RCP<V> lower_slack = lowerpv.get(SLACK);

    Teuchos::RCP<V> upper_opt   = upperpv.get(OPT);
    Teuchos::RCP<V> upper_slack = upperpv.get(SLACK);

    Elementwise::Fill<Real> setToMax(std::numeric_limits<Real>::max());
    Elementwise::Fill<Real> setToMin(std::numeric_limits<Real>::lowest());
    
    lower_opt->applyUnary( setToMin );
    lower_slack->zero();

    upper_opt->applyUnary( setToMax );
    upper_slack->applyUnary( setToMax );

    bc_ = Teuchos::rcp( new BoundConstraint<Real>( lower_, upper_ ) );

  }

  void project( Vector<Real> &x ) {
    bc_->project(x);
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
    bc_->pruneUpperActive( v, x, eps );
  } 

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, 
                         const Vector<Real> &x, Real eps ) {
    bc_->pruneUpperActive( v, g, x, eps );
  }

   void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps ) {
    bc_->pruneLowerActive( v, x, eps );
  } 

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, 
                         const Vector<Real> &x, Real eps ) {
    bc_->pruneLowerActive( v, g, x, eps );
  }
 
  void setVectorToUpperBound( Vector<Real> &u ) {
    u.set(*upper_);
  }

  void setVectorToLowerBound( Vector<Real> &l ) {
    l.set(*lower_);
  }

  bool isFeasible( const Vector<Real> &v ) { 
    return bc_->isFeasible(v);
  }

}; // class InteriorPointBoundConstraint


}

#endif
