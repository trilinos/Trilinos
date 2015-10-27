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

#include "ROL_PartitionedVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_BoundConstraint.hpp"

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

  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<Objective<Real> > barrier_;

  Real mu_;
  int nfval_;
  int ngval_;

  // Downcast const ROL::Vector to const ROL::PartitionedVector
  const PV& partition(const V& x) {
    using Teuchos::dyn_cast;
    return dyn_cast<const PV>(x);
  }
  
  // Downcast ROL::Vector to ROL::PartitionedVector
  PV& partition(V &x) {
    using Teuchos::dyn_cast;
    return dyn_cast<PV>(x);
  }
  

public:
  InteriorPointObjective( const Teuchos::RCP<Objective<Real> > &obj, 
                          const Teuchos::RCP<Objective<Real> > &barrier, Real mu ) :
    obj_(obj),barrier_(barrier),mu_(mu),nfval_(0),ngval_(0) {
  }
 
  void updatePenalty( Real mu ) {
    mu_ = mu;
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  } 
 
  int getNumberGradientEvaluations(void) {
    return ngval_;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const PV &xpv = partition(x); 

    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);

    ++nfval_;
    return obj_->value(*xopt,tol) + mu_*barrier_->value(*s,tol); 

  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    const PV &xpv = partition(x);
    PV &gpv = partition(g); 
        
    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);    

    Teuchos::RCP<V> gopt = gpv.get(OPT);
    Teuchos::RCP<V> gs   = gpv.get(SLACK);

    obj_->gradient(*gopt, *xopt, tol);
    barrier_->gradient(*gs, *s, tol); 
    gs->scale(mu_);

    ++ngval_;   
  } 
   
  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                 const Vector<Real> &x, Real &tol ) {

    const PV &xpv = partition(x);
    const PV &vpv = partition(v);
    PV &hvpv = partition(hv); 
        
    Teuchos::RCP<const V> xopt = xpv.get(OPT);    
    Teuchos::RCP<const V> s    = xpv.get(SLACK);    

    Teuchos::RCP<const V> vopt = vpv.get(OPT);    
    Teuchos::RCP<const V> vs   = vpv.get(SLACK);

    Teuchos::RCP<V> hvopt = hvpv.get(OPT);
    Teuchos::RCP<V> hvs   = hvpv.get(SLACK);

    obj_->hessVec(*hvopt, *vopt, *xopt, tol);
    barrier_->hessVec(*hvs, *vs, *s, tol);
    hvs->scale(mu_);
   
  }

}; // class InteriorPointObjective



/** @ingroup func_group
 *  \class ROL::InteriorPointEqualityConstraint
 *  \brief Has both inequality and equality constraints. 
 *        Treat inequality constraint as equality with slack variable
 */

template<class Real> 
class InteriorPointEqualityConstraint : public EqualityConstraint<Real> {
private:
  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type; 

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  const static size_type INEQ  = 0;
  const static size_type EQUAL = 1;

  Teuchos::RCP<EqualityConstraint<Real> > incon_;
  Teuchos::RCP<EqualityConstraint<Real> > eqcon_;
   
  bool hasEquality_;         // True if an equality constraint is present
  int  ncval_;               // Number of constraint evaluations

  // Downcast const ROL::Vector to const ROL::PartitionedVector
  const PV& partition(const V& x) {
    using Teuchos::dyn_cast;
    return dyn_cast<const PV>(x);
  }
  
  // Downcast ROL::Vector to ROL::PartitionedVector
  PV& partition(V &x) {
    using Teuchos::dyn_cast;
    return dyn_cast<PV>(x);
  }
 
public:

  // Constructor with inequality and equality constraints
  InteriorPointEqualityConstraint( const Teuchos::RCP<EqualityConstraint<Real> > &incon, 
                            const Teuchos::RCP<EqualityConstraint<Real> > &eqcon ) :
     incon_(incon), eqcon_(eqcon),hasEquality_(true), ncval_(0) {}

  // Constructor with inequality constraint only
  InteriorPointEqualityConstraint( const Teuchos::RCP<EqualityConstraint<Real> > &incon ) :
    incon_(incon), hasEquality_(false), ncval_(0) {}

 
  int getNumberConstraintEvaluations(void) {
    return ncval_;
  }
 

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    using Teuchos::RCP;

    // Partition vectors and extract subvectors
    const PV &xpv = partition(x);

    RCP<const V> xopt = xpv.get(OPT);    
    RCP<const V> s    = xpv.get(SLACK);

    PV &cpv = partition(c); 

    RCP<V> ci = cpv.get(INEQ);
    incon_->value(*ci, *xopt, tol);
    ci->axpy(-1.0,*s);

    if(hasEquality_) {
      RCP<V> ce = cpv.get(EQUAL);
      eqcon_->value(*ce, *xopt, tol);
    }

    ++ncval_;
  }

  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) {
    
    using Teuchos::RCP;

    // Partition vectors and extract subvectors
    const PV &xpv = partition(x);
    const PV &vpv = partition(v);

    RCP<const V> xopt = xpv.get(OPT);    
    RCP<const V> s    = xpv.get(SLACK);

    RCP<const V> vopt = vpv.get(OPT);    
    RCP<const V> vs   = vpv.get(SLACK);

    PV &jvpv = partition(jv); 

    RCP<V> jvi = jvpv.get(INEQ);
    incon_->applyJacobian(*jvi, *vopt, *xopt, tol);
    jvi->axpy(-1.0,*vs);

    if(hasEquality_) {
      RCP<V> jve = jvpv.get(EQUAL);
      eqcon_->applyJacobian(*jve, *vopt, *xopt, tol);  
    }

  }

  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) {

    using Teuchos::RCP;

    // Partition vectors and extract subvectors
    const PV &xpv = partition(x);
    PV &ajvpv = partition(ajv); 

    RCP<const V> xopt = xpv.get(OPT);    
    RCP<const V> s    = xpv.get(SLACK);

    RCP<V> ajvopt = ajvpv.get(OPT);
    RCP<V> ajvs   = ajvpv.get(SLACK);

    const PV &vpv = partition(v);

    RCP<const V> vi = vpv.get(INEQ);

    incon_->applyAdjointJacobian(*ajvopt,*vi,*xopt,tol);


    ajvs->set(*vi);
    ajvs->scale(-1.0);
   
    if(hasEquality_) {

      RCP<const V> ve = vpv.get(EQUAL);    
      RCP<V> temp = ajvopt->clone();
      eqcon_->applyAdjointJacobian(*temp,*ve,*xopt,tol);
      ajvopt->plus(*temp);

    } 

  }

  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) {
  
    using Teuchos::RCP;

    const PV &xpv = partition(x);
    const PV &vpv = partition(v);
    PV &ahuvpv = partition(ahuv); 

    RCP<const V> xopt = xpv.get(OPT);    
    RCP<const V> s    = xpv.get(SLACK);

    RCP<const V> vopt = vpv.get(OPT);    
    
    RCP<V> ahuvopt = ahuvpv.get(OPT);
    RCP<V> ahuvs   = ahuvpv.get(SLACK);

    RCP<V> temp = ahuvopt->clone();
 
    const PV &upv = partition(u);

    RCP<const V> ui   = upv.get(INEQ);
 
    incon_->applyAdjointHessian(*ahuvopt,*ui,*vopt,*xopt,tol);
    ahuvs->zero();

    if(hasEquality_) {
      RCP<const V> ue   = upv.get(EQUAL);    
      eqcon_->applyAdjointHessian(*temp,*ue,*vopt,*xopt,tol);
      ahuvopt->plus(*temp);
    }

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
