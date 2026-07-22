// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_INTERIORPOINTOBJECTIVE_H
#define ROL_INTERIORPOINTOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_ScalarController.hpp"

/** @ingroup func_group
    \class ROL::PrimalInteriorPointObjective
    \brief Provides the interface to evaluate the Interior Pointy
           log barrier penalty function with upper and lower bounds on
           some elements

    ---
*/

namespace ROL {

template<class Real>
class InteriorPointObjective : public Objective<Real> {

  typedef Elementwise::ValueSet<Real>   ValueSet;

private:
  
  const Ptr<Objective<Real>>       obj_;
  const Ptr<BoundConstraint<Real>> bnd_;
  const Ptr<const Vector<Real>>    lo_;
  const Ptr<const Vector<Real>>    up_;

  Ptr<Vector<Real>> maskL_;   // Elements are 1 when xl>-INF, zero for xl =-INF
  Ptr<Vector<Real>> maskU_;   // Elements are 1 when xu< INF, zero for xu = INF
  Ptr<Vector<Real>> maskL0_;  // Elements are 1 when xl>-INF and xu = INF, zero for xl =-INF
  Ptr<Vector<Real>> maskU0_;  // Elements are 1 when xu< INF and XL =-INF, zero for xu = INF
  Ptr<Vector<Real>> pwa_;     // Scratch vector

  bool useLinearDamping_;     // Add linear damping terms to the penalized objective
                              // to prevent the problems such as when the log barrier
                              // contribution is unbounded below on the feasible set
  Real kappaD_;               // Linear damping coefficient
  Real mu_;                   // Penalty parameter

  Ptr<ScalarController<Real,int>> fval_;
  Ptr<VectorController<Real,int>> gradient_;

  int nfval_;
  int ngrad_;

  // x <- f(x) = { log(x) if x >  0
  //             { -inf   if x <= 0
  class ModifiedLogarithm : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply( const Real &x ) const {
      const Real zero(0), NINF(ROL_NINF<Real>());
      return (x>zero) ? std::log(x) : NINF;
      //return std::log(x);
    }
  }; // class ModifiedLogarithm

  // x <- f(x) = { 1/x  if  x >  0
  //             { 0    if  x <= 0
  class ModifiedReciprocal : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply( const Real &x ) const {
      const Real zero(0), one(1);
      return (x>zero) ? one/x : zero;
      //return one/x;
    }

  }; // class ModifiedReciprocal

  // x <- f(x,y) = { y/x  if  x >  0
  //               { 0    if  x <= 0
  class ModifiedDivide : public Elementwise::BinaryFunction<Real> {
  public:
    Real apply( const Real &x, const Real &y ) const {
      const Real zero(0);
      return (x>zero) ? y/x : zero;
      //return y/x;
    }
  }; // class ModifiedDivide

  // x <- f(x,y) = { x  if  y != 0, complement == false
  //               { 0  if  y == 0, complement == false
  //               { 0  if  y != 0, complement == true
  //               { x  if  y == 0, complement == true
  class Mask : public Elementwise::BinaryFunction<Real> {
  private:
    bool complement_;
  public:
    Mask( bool complement ) : complement_(complement) {}
    Real apply( const Real &x, const Real &y ) const {
      const Real zero(0);
      return ( complement_ ^ (y != zero) ) ? zero : x;
    }
  }; // class Mask

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    const Real zero(0), one(1);

    fval_     = makePtr<ScalarController<Real,int>>();
    gradient_ = makePtr<VectorController<Real,int>>();

    // Determine the index sets where the
    ValueSet isBoundedBelow( ROL_NINF<Real>(), ValueSet::GREATER_THAN, one, zero );
    ValueSet isBoundedAbove( ROL_INF<Real>(),  ValueSet::LESS_THAN,    one, zero );

    maskL_ = x.clone(); maskL_->applyBinary(isBoundedBelow,*lo_);
    maskU_ = x.clone(); maskU_->applyBinary(isBoundedAbove,*up_);

    pwa_ = x.clone();

    if( useLinearDamping_ ) {
      maskL0_ = x.clone();
      maskL0_->set(*maskL_);                     // c_i = { 1 if l_i > NINF
                                                 //       { 0 otherwise
      maskL0_->applyBinary(Mask(true),*maskU_);  // c_i = { 1 if l_i > NINF and u_i = INF
                                                 //       { 0 otherwise
      maskU0_ = x.clone();
      maskU0_->set(*maskU_);                     // c_i = { 1 if u_i < INF
                                                 //       { 0 otherwise
      maskU0_->applyBinary(Mask(true),*maskL_);  // c_i = { 1 if u_i < INF and l_i = NINF
                                                 //       { 0 otherwise
    }
  }

public:

  InteriorPointObjective( const Ptr<Objective<Real>>       &obj,
                          const Ptr<BoundConstraint<Real>> &bnd, 
                          const Vector<Real>               &x,
                          const Vector<Real>               &g,
                          const bool useLinearDamping,
                          const Real kappaD,
                          const Real mu )
    : obj_(obj), bnd_(bnd), lo_(bnd->getLowerBound()), up_(bnd->getUpperBound()),
      useLinearDamping_(useLinearDamping), kappaD_(kappaD), mu_(mu),
      nfval_(0), ngrad_(0) {
    initialize(x,g);
  }

  InteriorPointObjective( const Ptr<Objective<Real>>       &obj,
                          const Ptr<BoundConstraint<Real>> &bnd, 
                          const Vector<Real>               &x,
                          const Vector<Real>               &g,
                          ParameterList                    &parlist )
    : obj_(obj), bnd_(bnd), lo_(bnd->getLowerBound()), up_(bnd->getUpperBound()),
      nfval_(0), ngrad_(0) {
    ParameterList &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");
    ParameterList &lblist = iplist.sublist("Barrier Objective");

    useLinearDamping_ = lblist.get("Use Linear Damping",         true);
    kappaD_           = lblist.get("Linear Damping Coefficient", 1.e-4);
    mu_               = lblist.get("Initial Barrier Parameter",  0.1);

    initialize(x,g);
  }

  Real getObjectiveValue(const Vector<Real> &x, Real &tol) {
    int key(0);
    Real val(0);
    bool isComputed = fval_->get(val,key);
    if (!isComputed) {
      val = obj_->value(x,tol); nfval_++;
      fval_->set(val,key);
    }
    return val;
  }

  const Ptr<const Vector<Real>> getObjectiveGradient(const Vector<Real> &x, Real &tol) {
    int key(0);
    if (!gradient_->isComputed(key)) {
      if (gradient_->isNull(key)) gradient_->allocate(x.dual(),key);
      obj_->gradient(*gradient_->set(key),x,tol); ngrad_++;
    }
    return gradient_->get(key);
  }

  int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) const {
    return ngrad_;
  }

  void updatePenalty(const Real mu) {
    mu_ = mu;
  }

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) {
    obj_->update(x,type,iter);
    fval_->objectiveUpdate(type);
    gradient_->objectiveUpdate(type);
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    const Real zero(0), one(1);
    Real linearTerm = zero;
    // Compute the unpenalized objective value
    Real fval = getObjectiveValue(x,tol);
    // Compute log barrier
    ModifiedLogarithm                mlog;
    Elementwise::ReductionSum<Real>  sum;
    Elementwise::Multiply<Real>      mult;

    pwa_->set(x);                             // pwa = x
    pwa_->axpy(-one,*lo_);                    // pwa = x-l
    if( useLinearDamping_ ) {
      // Penalizes large positive x_i when only a lower bound exists
      linearTerm += maskL0_->dot(*pwa_);
    }
    pwa_->applyUnary(mlog);                   // pwa = mlog(x-l)
    Real aval = pwa_->dot(*maskL_);

    pwa_->set(*up_);                          // pwa = u
    pwa_->axpy(-one,x);                       // pwa = u-x
    if( useLinearDamping_ ) {
      // Penalizes large negative x_i when only an upper bound exists
      linearTerm += maskU0_->dot(*pwa_);
    }
    pwa_->applyUnary(mlog);                   // pwa = mlog(u-x)
    Real bval = pwa_->dot(*maskU_);

    fval -= mu_*(aval+bval);
    fval += kappaD_*mu_*linearTerm;
    return fval;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    const Real one(1);
    // Compute gradient of objective function
    g.set(*getObjectiveGradient(x,tol));

    // Add gradient of the log barrier penalty
    ModifiedReciprocal mrec;

    pwa_->set(x);                             // pwa = x
    pwa_->axpy(-one,*lo_);                    // pwa = x-l
    pwa_->applyUnary(mrec);                   // pwa_i = 1/(x_i-l_i) for i s.t. x_i > l_i, 0 otherwise
    pwa_->applyBinary(Mask(true),*maskL_);    // zero elements where l = NINF
    g.axpy(-mu_,pwa_->dual());
    if( useLinearDamping_ ) {
      g.axpy(-mu_*kappaD_,maskL0_->dual());
    }

    pwa_->set(*up_);                          // pwa = u
    pwa_->axpy(-one,x);                       // pwa = u-x
    pwa_->applyUnary(mrec);                   // pwa_i = 1/(u_i-x_i) for i s.t. u_i > x_i, 0 otherwise
    pwa_->applyBinary(Mask(true),*maskU_);    // zero elements where u = INF
    g.axpy( mu_,pwa_->dual());
    if( useLinearDamping_ ) {
      g.axpy( mu_*kappaD_,maskU0_->dual());
    }
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    const Real one(1), two(2);
    // Evaluate objective hessian
    obj_->hessVec(hv,v,x,tol);

    // Evaluate log barrier hessian
    ModifiedReciprocal mrec;
    Elementwise::Multiply<Real> mult;
    Elementwise::Power<Real> square(two);

    pwa_->set(x);                             // pwa = x
    pwa_->axpy(-one,*lo_);                    // pwa = x-l
    pwa_->applyUnary(mrec);                   // pwa_i = 1/(x_i-l_i) for i s.t. x_i > l_i, 0 otherwise
    pwa_->applyBinary(Mask(true),*maskL_);    // zero elements where l = NINF
    pwa_->applyUnary(square);                 // pwa_i = { (x_i-l_i)^(-2)  if l_i > NINF
                                              //         { 0               if l_i = NINF
    pwa_->applyBinary(mult,v);
    hv.axpy(mu_,pwa_->dual());

    pwa_->set(*up_);                          // pwa = u
    pwa_->axpy(-one,x);                       // pwa = u-x
    pwa_->applyUnary(mrec);                   // pwa_i = 1/(u_i-x_i) for i s.t. u_i > x_i, 0 otherwise
    pwa_->applyBinary(Mask(true),*maskU_);    // zero elements where u = INF
    pwa_->applyUnary(square);                 // pwa_i = { (u_i-x_i)^(-2)  if u_i < INF
                                              //         { 0               if u_i = INF
    pwa_->applyBinary(mult,v);
    hv.axpy(mu_,pwa_->dual());
 }

}; // class InteriorPointObjective

}

#endif // ROL_INTERIORPOINTOBJECTIVE_H
