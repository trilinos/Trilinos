// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_INTERIORPOINTPENALTY_H
#define ROL_INTERIORPOINTPENALTY_H

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_ParameterList.hpp"

/** @ingroup func_group
    \class ROL::InteriorPointPenalty
    \brief Provides the interface to evaluate the Interior Pointy
           log barrier penalty function with upper and lower bounds on
           some elements

    ---
*/

namespace ROL {

template<class Real>
class InteriorPointPenalty : public Objective<Real> {

  typedef Vector<Real>             V;
  typedef Objective<Real>          OBJ;
  typedef BoundConstraint<Real>    BND;

  typedef Elementwise::ValueSet<Real>   ValueSet;

private:
  
  const ROL::Ptr<OBJ>       obj_;
  const ROL::Ptr<BND>       bnd_;
  const ROL::Ptr<const V>   lo_;
  const ROL::Ptr<const V>   up_;

  ROL::Ptr<V>   g_;       // Gradient of penalized objective

  ROL::Ptr<V>   maskL_;   // Elements are 1 when xl>-INF, zero for xl = -INF
  ROL::Ptr<V>   maskU_;   // Elements are 1 when xu< INF, zero for xu =  INF

  ROL::Ptr<V>   a_;       // Scratch vector
  ROL::Ptr<V>   b_;       // Scratch vector
  ROL::Ptr<V>   c_;       // Scratch vector

  bool useLinearDamping_;     // Add linear damping terms to the penalized objective
                              // to prevent the problems such as when the log barrier
                              // contribution is unbounded below on the feasible set

  Real mu_;                   // Penalty parameter
  Real kappaD_;               // Linear damping coefficient
  Real fval_;                 // Unpenalized objective value

  int nfval_;
  int ngval_;



  // x <- f(x) = { log(x) if x >  0
  //             { 0      if x <= 0
  class ModifiedLogarithm : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply( const Real &x ) const {
      return (x>0) ? std::log(x) : Real(0.0);
    }
  }; // class ModifiedLogarithm

  // x <- f(x) = { 1/x  if  x >  0
  //             { 0    if  x <= 0
  class ModifiedReciprocal : public Elementwise::UnaryFunction<Real> {
  public:
    Real apply( const Real &x ) const {
      return (x>0) ? 1.0/x : Real(0.0);
    }

  }; // class ModifiedReciprocal

  // x <- f(x,y) = { y/x  if  x >  0
  //               { 0    if  x <= 0
  class ModifiedDivide : public Elementwise::BinaryFunction<Real> {
  public:
    Real apply( const Real &x, const Real &y ) const {
      return (x>0) ? y/x : Real(0.0);
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
      return ( complement_ ^ (y !=0) ) ? 0 : x;
    }
  }; // class Mask


public:

  ~InteriorPointPenalty() {}

  InteriorPointPenalty( const ROL::Ptr<Objective<Real> > &obj,
                        const ROL::Ptr<BoundConstraint<Real> > &con, 
                        ROL::ParameterList &parlist ) :
    obj_(obj), bnd_(con), lo_( con->getLowerBound() ), up_( con->getUpperBound() ) {

    Real one(1);
    Real zero(0);

    // Determine the index sets where the
    ValueSet isBoundedBelow( ROL_NINF<Real>(), ValueSet::GREATER_THAN, one, zero );
    ValueSet isBoundedAbove( ROL_INF<Real>(),  ValueSet::LESS_THAN,    one, zero );

    maskL_ = lo_->clone();
    maskU_ = up_->clone();

    maskL_->applyBinary(isBoundedBelow,*lo_);
    maskU_->applyBinary(isBoundedAbove,*up_);

    ROL::ParameterList &iplist = parlist.sublist("Step").sublist("Primal Dual Interior Point");
    ROL::ParameterList &lblist = iplist.sublist("Barrier Objective");

    useLinearDamping_ = lblist.get("Use Linear Damping",true);
    kappaD_ = lblist.get("Linear Damping Coefficient",1.e-4);
    mu_ = lblist.get("Initial Barrier Parameter",0.1);


    a_ = lo_->clone();
    b_ = up_->clone();
    g_ = lo_->dual().clone();

    if( useLinearDamping_ ) {
      c_ = lo_->clone();
    }
  }

  Real getObjectiveValue(void) const {
    return fval_;
  }

  ROL::Ptr<Vector<Real> > getGradient(void) const {
    return g_;
  }

  int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) const {
    return ngval_;
  }

  ROL::Ptr<const Vector<Real> > getLowerMask(void) const {
    return maskL_;
  }

  ROL::Ptr<const Vector<Real> > getUpperMask(void) const {
    return maskU_;
  }

  /** \brief Update the interior point penalized objective.

      This function updates the log barrier penalized function at new iterations.
      @param[in]          x      is the new iterate.
      @param[in]          flag   is true if the iterate has changed.
      @param[in]          iter   is the outer algorithm iterations count.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
  }

  /** \brief Compute value.

      This function returns the log barrier penalized objective value.

      \f[ \varphi_\mu(x) = f(x) - \mu \sum\limits_{i\int I_L} \ln(x_i-l_i)
                                - \mu \sum\limits_{i\in I_Y} \ln(u_i-x_i) \f]
      Where \f$ I_L=\{i:l_i>-\infty\} \f$ and \f$ I_U = \{i:u_i<\infty\}\f$

      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for interior point penalty computation.
  */
  Real value( const Vector<Real> &x, Real &tol ) {

    ModifiedLogarithm                mlog;
    Elementwise::ReductionSum<Real>  sum;
    Elementwise::Multiply<Real>      mult;

    // Compute the unpenalized objective value
    fval_ = obj_->value(x,tol);
    nfval_++;

    Real fval = fval_;
    Real linearTerm = 0.0;   // Linear damping contribution

    a_->set(x);                             // a_i = x_i
    a_->axpy(-1.0,*lo_);                    // a_i = x_i-l_i

    if( useLinearDamping_ ) {

      c_->set(*maskL_);                     // c_i = { 1 if l_i > NINF
                                            //       { 0 otherwise
      c_->applyBinary(Mask(true),*maskU_);  // c_i = { 1 if l_i > NINF and u_i = INF
                                            //       { 0 otherwise
      c_->applyBinary(mult,*a_);            // c_i = { x_i-l_i if l_i > NINF and u_i = INF

      // Penalizes large positive x_i when only a lower bound exists
      linearTerm += c_->reduce(sum);
    }

    a_->applyUnary(mlog);                   // a_i = mlog(x_i-l_i)

    Real aval = a_->dot(*maskL_);

    b_->set(*up_);                          // b_i = u_i
    b_->axpy(-1.0,x);                       // b_i = u_i-x_i

    if( useLinearDamping_ ) {

      c_->set(*maskU_);                     // c_i = { 1 if u_i < INF
                                            //       { 0 otherwise
      c_->applyBinary(Mask(true),*maskL_);  // c_i = { 1 if u_i < INF and l_i = NINF
                                            //       { 0 otherwise
      c_->applyBinary(mult,*b_);            // c_i = { u_i-x_i if u_i < INF and l_i = NINF
                                            //       { 0 otherwise

      // Penalizes large negative x_i when only an upper bound exists
      linearTerm += c_->reduce(sum);

    }

    b_->applyUnary(mlog);                   // b_i = mlog(u_i-x_i)

    Real bval = b_->dot(*maskU_);


    fval -= mu_*(aval+bval);
    fval += kappaD_*mu_*linearTerm;

    return fval;

  }

  /** \brief Compute gradient.

      This function returns the log barrier penalized gradient.
      @param[out]         g   is the gradient.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact log barrier penalty computation.
  */

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Compute gradient of objective function
    obj_->gradient(*g_,x,tol);
    ngval_++;
    g.set(*g_);

    // Add gradient of the log barrier penalty
    ModifiedReciprocal mrec;

    a_->set(x);                             // a = x
    a_->axpy(-1.0,*lo_);                    // a = x-l

    a_->applyUnary(mrec);                   // a_i = 1/(x_i-l_i) for i s.t. x_i > l_i, 0 otherwise
    a_->applyBinary(Mask(true),*maskL_);   // zero elements where l = NINF

    b_->set(*up_);                          // b = u
    b_->axpy(-1.0,x);                       // b = u-x
    b_->applyUnary(mrec);                   // b_i = 1/(u_i-x_i) for i s.t. u_i > x_i, 0 otherwise
    b_->applyBinary(Mask(true),*maskU_);   // zero elements where u = INF

    g.axpy(-mu_,*a_);
    g.axpy(mu_,*b_);

    if( useLinearDamping_ ) {

      a_->set(*maskL_);
      a_->applyBinary(Mask(true),*maskU_);  // a_i = { 1 if l_i > NINF and u_i = INF
                                            //       { 0 otherwise
      b_->set(*maskU_);
      b_->applyBinary(Mask(true),*maskL_);  // b_i = { 1 if u_i < INF and l_i = NINF
                                            //       { 0 otherwise
      g.axpy(-mu_*kappaD_,*a_);
      g.axpy( mu_*kappaD_,*b_);

    }
  }

  /** \brief Compute action of Hessian on vector.

      This function returns the log barrier penalized Hessian acting on a given vector.
      @param[out]         hv  is the Hessian-vector product.
      @param[in]          v   is the given vector.
      @param[in]          x   is the current iterate.
      @param[in]          tol is a tolerance for inexact log barrier penalty computation.
  */

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

    ModifiedReciprocal mrec;
    Elementwise::Multiply<Real> mult;
    Elementwise::Power<Real> square(2.0);

    obj_->hessVec(hv,v,x,tol);

    a_->set(x);                             // a = x
    a_->axpy(-1.0,*lo_);                    // a = x-l
    a_->applyUnary(mrec);                   // a_i = 1/(x_i-l_i) for i s.t. x_i > l_i, 0 otherwise
    a_->applyBinary(Mask(true),*maskL_);    // zero elements where l = NINF
    a_->applyUnary(square);                 // a_i = { (x_i-l_i)^(-2)  if l_i > NINF
                                            //       { 0               if l_i = NINF
    a_->applyBinary(mult,v);

    b_->set(*up_);                          // b = u
    b_->axpy(-1.0,x);                       // b = u-x
    b_->applyUnary(mrec);                   // b_i = 1/(u_i-x_i) for i s.t. u_i > x_i, 0 otherwise
    b_->applyBinary(Mask(true),*maskU_);    // zero elements where u = INF
    b_->applyUnary(square);                 // b_i = { (u_i-x_i)^(-2)  if u_i < INF
                                            //       { 0               if u_i = INF
    b_->applyBinary(mult,v);

    hv.axpy(mu_,*a_);
    hv.axpy(mu_,*b_);
 }

  // Return the unpenalized objective
  const ROL::Ptr<OBJ> getObjective( void ) { return obj_; }

  // Return the bound constraint
  const ROL::Ptr<BND> getBoundConstraint( void ) { return bnd_; }

}; // class InteriorPointPenalty

}


#endif // ROL_INTERIORPOINTPENALTY_H
