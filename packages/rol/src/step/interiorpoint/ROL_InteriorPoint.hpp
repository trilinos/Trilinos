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

namespace InteriorPoint {

/** @ingroup func_group
 *  \class ROL::PenalizedObjective
 *  \brief Adds barrier term to generic objective
 */
template <class Real>
class PenalizedObjective : public ROL::Objective<Real> {
private:

  typedef Vector<Real>            V;
  typedef PartitionedVector<Real> PV;
  typedef typename PV::size_type  size_type; 

  const static size_type OPT   = 0;
  const static size_type SLACK = 1;

  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<Objective<Real> > barrier_;
  Teuchos::RCP<V> go_;
  Teuchos::RCP<V> gs_;


  Real mu_;
  Real fval_;     // Stored raw objective value
  Real pval_;     // Stored penalty objective value
  int nfval_;
  int ngval_;

  bool isValueComputed_;
  bool isGradientComputed_;

public:

  PenalizedObjective( Objective<Real> &obj, 
                      Objective<Real> &barrier, 
                      const Vector<Real> &x,
                      Real mu ) :
    obj_(Teuchos::rcp(&obj,false)),
    barrier_(Teuchos::rcp(&barrier,false)), 
    go_(Teuchos::null), gs_(Teuchos::null),
    mu_(mu), fval_(0.0), pval_(0.0), nfval_(0), ngval_(0),
    isValueComputed_(false), isGradientComputed_(false) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    RCP<V> g = x.dual().clone();
 
    PV &gpv = dyn_cast<PV>(*g);

    go_ = gpv.get(OPT);
    gs_ = gpv.get(SLACK);

  }
 
  void updatePenalty( Real mu ) {
    mu_ = mu;
  }

  Real getObjectiveValue( const Vector<Real> &x ) {
    Real tol = std::sqrt(ROL_EPSILON);
    if( !isValueComputed_ ) {
      // Evaluate (unpenalized) objective function value
      const PV &xpv = Teuchos::dyn_cast<const PV>(x);
      Teuchos::RCP<const Vector<Real> > xo = xpv.get(OPT);
      fval_ = obj_->value(*xo,tol);
      ++nfval_;
      isValueComputed_ = true;
    } 
    return fval_;
  }

  void getObjectiveGradient( Vector<Real> &g, const Vector<Real> &x ) {
    Real tol = std::sqrt(ROL_EPSILON);
    if( !isGradientComputed_ ) {
      // Evaluate (unpenalize) objective function gradient
      const PV &xpv = Teuchos::dyn_cast<const PV>(x);
      PV &gpv = Teuchos::dyn_cast<PV>(g);
      Teuchos::RCP<const Vector<Real> > xo = xpv.get(OPT);
      Teuchos::RCP<Vector<Real> > go = gpv.get(OPT);
      obj_->gradient(*go,*xo,tol);
      ++ngval_;
      isGradientComputed_ = true;
    }
  }

  int getNumberFunctionEvaluations(void) {
    return nfval_;
  } 
 
  int getNumberGradientEvaluations(void) {
    return ngval_;
  }

  void reset(void) { 
    nfval_ = 0.; nfval_ = 0.;
  } 

  /** \brief Update barrier penalized objective function

      This function updates the penalized objective function at new iterations. 
      @param[in]          x      is the new iterate. 
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
   
    const PV &xpv = Teuchos::dyn_cast<const PV>(x);

    Teuchos::RCP<const V> xo = xpv.get(OPT);
    Teuchos::RCP<const V> xs = xpv.get(SLACK);

    obj_->update(*xo,flag,iter);
    barrier_->update(*xs,flag,iter);

    if ( flag ) {
      isValueComputed_    = false;
      isGradientComputed_ = false;
    }     
  }

  /** \brief Compute value.

      This function returns the barrier objective value.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance. 
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    
    const PV &xpv = Teuchos::dyn_cast<const PV>(x); 

    Teuchos::RCP<const V> xo = xpv.get(OPT);    
    Teuchos::RCP<const V> xs = xpv.get(SLACK);

    if ( !isValueComputed_ ) {
      // Compute objective function value
      fval_ = obj_->value(*xo,tol);  
      pval_ = barrier_->value(*xs,tol);

      ++nfval_;
      isValueComputed_ = true;
    }

    Real val = fval_ + mu_*pval_; 

    return val; 
  }


  /** \brief Compute gradient.

      This function returns the barrier penalized objective gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance. 
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

    if ( !isGradientComputed_ ) {

      // Compute objective function gradient
      const PV &xpv = Teuchos::dyn_cast<const PV>(x);
        
      Teuchos::RCP<const V> xo = xpv.get(OPT);    
      Teuchos::RCP<const V> xs = xpv.get(SLACK);    

      obj_->gradient(*go_,*xo,tol);
      barrier_->gradient(*gs_,*xs,tol);
      ++ngval_;
      isGradientComputed_ = true;
    }

    PV &gpv = Teuchos::dyn_cast<PV>(g); 

    gpv.set(OPT,*go_); 
    gpv.set(SLACK,*gs_);

  } 

  /** \brief Apply Hessian approximation to vector.

      This function applies the Hessian of the barrier penalized objective 
      to the vector \f$v\f$.
      @param[out]         hv  is the the action of the Hessian on \f$v\f$.
      @param[in]          v   is the direction vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance.
  */   
  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                 const Vector<Real> &x, Real &tol ) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    const PV &xpv = dyn_cast<const PV>(x);
    const PV &vpv = dyn_cast<const PV>(v);
    PV &hvpv = dyn_cast<PV>(hv); 
        
    RCP<const V> xo = xpv.get(OPT);    
    RCP<const V> xs = xpv.get(SLACK);    

    RCP<const V> vo = vpv.get(OPT);    
    RCP<const V> vs = vpv.get(SLACK);

    RCP<V> hvo = hvpv.get(OPT);
    RCP<V> hvs = hvpv.get(SLACK);

    obj_->hessVec(*hvo, *vo, *xo, tol);
    barrier_->hessVec(*hvs, *vs, *xs, tol);
    hvs->scale(mu_);
   
  }

}; // class InteriorPointObjective



/** @ingroup func_group
 *  \class ROL::InteriorPoint::CompositeConstraint
 *  \brief Has both inequality and equality constraints. 
 *        Treat inequality constraint as equality with slack variable
 */

template<class Real> 
class CompositeConstraint : public EqualityConstraint<Real> {
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
   
  Teuchos::RCP<V> ci_;
  Teuchos::RCP<V> ce_;

  bool hasEquality_;         // True if an equality constraint is present
  int  ncval_;               // Number of constraint evaluations

  bool isConstraintComputed_;


public:

  // Constructor with inequality and equality constraints
  CompositeConstraint( EqualityConstraint<Real> &incon, 
                       EqualityConstraint<Real> &eqcon,
                       const Vector<Real> &c ) :
                       incon_(Teuchos::rcp(&incon,false)), 
                       eqcon_(Teuchos::rcp(&eqcon,false)), 
                       ci_(Teuchos::null), ce_(Teuchos::null), 
                       hasEquality_(true), ncval_(0), 
                       isConstraintComputed_(false) {

    const PV &cpv = Teuchos::dyn_cast<const PV>(c);

    ci_ = cpv.get(INEQ)->clone(); 
    ce_ = cpv.get(EQUAL)->clone();
  }

  // Constructor with inequality constraint only
  CompositeConstraint( EqualityConstraint<Real> &incon,
                       const Vector<Real> &c ) :
                       incon_(Teuchos::rcp(&incon,false)), 
                       eqcon_(Teuchos::null),
                       ci_(Teuchos::null), hasEquality_(false), 
                       ncval_(0), isConstraintComputed_(false) {

    const PV &cpv = Teuchos::dyn_cast<const PV>(c);

    ci_ = cpv.get(INEQ)->clone(); 
  }
 
  int getNumberConstraintEvaluations(void) {
    return ncval_;
  }
 
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {

    const PV &xpv = Teuchos::dyn_cast<const PV>(x);     
   
    Teuchos::RCP<const V> xo = xpv.get(OPT);
    Teuchos::RCP<const V> xs = xpv.get(SLACK);

    incon_->update(*xo,flag,iter);
 
    if( hasEquality_ ) {
      eqcon_->update(*xo,flag,iter);
    }

    if( flag ) {
      isConstraintComputed_ = false;
    }
 
  }

  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) {

    PV &cpv = Teuchos::dyn_cast<PV>(c); 

    if( !isConstraintComputed_ ) {
      
      // Compute constraint contributions
      const PV &xpv = Teuchos::dyn_cast<const PV>(x);
        

      Teuchos::RCP<const V> xo = xpv.get(OPT);    
      Teuchos::RCP<const V> xs = xpv.get(SLACK);
 
      incon_->value(*ci_, *xo, tol);
      ci_->axpy(-1.0,*xs);
       
      if(hasEquality_) {
        Teuchos::RCP<V> ce = cpv.get(EQUAL);
        eqcon_->value(*ce, *xo, tol);
      }

      ++ncval_;
      isConstraintComputed_ = true; 
    }

    cpv.set(INEQ,*ci_);
  
    if(hasEquality_) {
      cpv.set(EQUAL,*ce_);
    }
  }

  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) {
    
    using Teuchos::RCP;  using Teuchos::dyn_cast;

    // Partition vectors and extract subvectors
    const PV &xpv = dyn_cast<const PV>(x);
    const PV &vpv = dyn_cast<const PV>(v);

    RCP<const V> xo = xpv.get(OPT);    
    RCP<const V> xs = xpv.get(SLACK);

    RCP<const V> vo = vpv.get(OPT);    
    RCP<const V> vs = vpv.get(SLACK);

    PV &jvpv = dyn_cast<PV>(jv); 

    RCP<V> jvi = jvpv.get(INEQ);
    incon_->applyJacobian(*jvi, *vo, *xo, tol);
    jvi->axpy(-1.0,*vs);

    if(hasEquality_) {
      RCP<V> jve = jvpv.get(EQUAL);
      eqcon_->applyJacobian(*jve, *vo, *xo, tol);  
    }

  }

  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) {

    using Teuchos::RCP;  using Teuchos::dyn_cast;

    // Partition vectors and extract subvectors
    const PV &xpv = dyn_cast<const PV>(x);
    PV &ajvpv = dyn_cast<PV>(ajv); 

    RCP<const V> xo = xpv.get(OPT);    
    RCP<const V> xs = xpv.get(SLACK);

    RCP<V> ajvo = ajvpv.get(OPT);
    RCP<V> ajvs = ajvpv.get(SLACK);

    const PV &vpv = dyn_cast<const PV>(v);

    RCP<const V> vi = vpv.get(INEQ);

    incon_->applyAdjointJacobian(*ajvo,*vi,*xo,tol);

    ajvs->set(*vi);
    ajvs->scale(-1.0);
   
    if(hasEquality_) {

      RCP<const V> ve = vpv.get(EQUAL);    
      RCP<V> temp = ajvo->clone();
      eqcon_->applyAdjointJacobian(*temp,*ve,*xo,tol);
      ajvo->plus(*temp);

    } 

  }

  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) {
  
    using Teuchos::RCP;  using Teuchos::dyn_cast;

    const PV &xpv = dyn_cast<const PV>(x);
    const PV &vpv = dyn_cast<const PV>(v);
    PV &ahuvpv = dyn_cast<PV>(ahuv); 

    RCP<const V> xo = xpv.get(OPT);    
    RCP<const V> xs = xpv.get(SLACK);

    RCP<const V> vo = vpv.get(OPT);    
    
    RCP<V> ahuvo = ahuvpv.get(OPT);
    RCP<V> ahuvs = ahuvpv.get(SLACK);

    RCP<V> temp = ahuvo->clone();
 
    const PV &upv = dyn_cast<const PV>(u);

    RCP<const V> ui = upv.get(INEQ);
 
    incon_->applyAdjointHessian(*ahuvo,*ui,*vo,*xo,tol);
    ahuvs->zero();

    if(hasEquality_) {
      RCP<const V> ue   = upv.get(EQUAL);    
      eqcon_->applyAdjointHessian(*temp,*ue,*vo,*xo,tol);
      ahuvo->plus(*temp);
    }

  }
   
}; // class CompositeConstraint

} // namespace InteriorPoint
} // namespace ROL

#endif
