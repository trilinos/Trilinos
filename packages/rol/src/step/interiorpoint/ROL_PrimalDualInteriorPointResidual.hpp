// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PRIMALDUALINTERIORPOINTRESIDUAL_H
#define ROL_PRIMALDUALINTERIORPOINTRESIDUAL_H

#include "ROL_BoundConstraint.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Objective.hpp"
#include "ROL_PartitionedVector.hpp"
#include "ROL_RandomVector.hpp"

#include <iostream>

/** @ingroup func_group
    \class   ROL::PrimalDualInteriorPointResidual
    \brief   Symmetrized form of the KKT operator for the Type-EB problem
             with equality and bound multipliers

    The system is symmetrized by multiplying through by

    S = [  I        0        0       0      ]
        [  0        I        0       0      ]
        [  0        0    -inv(Zl)    0      ]
        [  0        0        0     -inv(Zu) ] 
    
    Where Zl and Zu are diagonal matrices containing the lower and upper 
    bound multipliers respectively

    Infinite bounds have identically zero-valued lagrange multipliers.
    
    --- 
*/

namespace ROL { 

template<class Real> 
class PrimalDualInteriorPointResidual : public Constraint<Real> {

  typedef ROL::ParameterList     PL;

  typedef Vector<Real>               V;
  typedef PartitionedVector<Real>    PV;
  typedef Objective<Real>            OBJ;
  typedef Constraint<Real>           CON;
  typedef BoundConstraint<Real>      BND;

  typedef Elementwise::ValueSet<Real>   ValueSet;

  typedef typename PV::size_type size_type;

private: 

  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 2;
  static const size_type UPPER = 3;

  const ROL::Ptr<OBJ> obj_;
  const ROL::Ptr<CON> con_;
  const ROL::Ptr<BND> bnd_;

  ROL::Ptr<const V>   x_;            // Optimization vector
  ROL::Ptr<const V>   l_;            //  constraint multiplier
  ROL::Ptr<const V>   zl_;           // Lower bound multiplier
  ROL::Ptr<const V>   zu_;           // Upper bound multiplier

  ROL::Ptr<const V>   xl_;           // Lower bound
  ROL::Ptr<const V>   xu_;           // Upper bound 

  const ROL::Ptr<const V> maskL_;   
  const ROL::Ptr<const V> maskU_;

  ROL::Ptr<V> scratch_;              // Scratch vector the same dimension as x

  Real mu_;

  bool symmetrize_;

  Real one_;
  Real zero_;
  
  int nfval_;
  int ngrad_;
  int ncval_;

  Elementwise::Multiply<Real>  mult_;

  class SafeDivide : public Elementwise::BinaryFunction<Real> {
  public:
    Real apply( const Real &x, const Real &y ) const {
      return y != 0 ? x/y : 0;
    }
  }; 
  
  SafeDivide div_;

  class SetZeros : public Elementwise::BinaryFunction<Real> {
  public:
    Real apply( const Real &x, const Real &y ) const {
      return y==1.0 ? 0 : x;
    }
  };

  SetZeros setZeros_;  

  // Fill in zeros of x with corresponding values of y
  class InFill : public Elementwise::BinaryFunction<Real> {
  public:
    Real apply( const Real &x, const Real &y ) const {
      return x == 0 ? y : x;
    }
  };
   
  InFill inFill_;
 
  // Extract the optimization and lagrange multiplier
  ROL::Ptr<V> getOptMult( V &vec ) {
    PV &vec_pv = dynamic_cast<PV&>(vec);
 
    return CreatePartitioned(vec_pv.get(OPT),vec_pv.get(EQUAL));
  }

  // Extract the optimization and lagrange multiplier
  ROL::Ptr<const V> getOptMult( const V &vec ) {
    const PV &vec_pv = dynamic_cast<const PV&>(vec);
 
    return CreatePartitioned(vec_pv.get(OPT),vec_pv.get(EQUAL));
  }




public:

  PrimalDualInteriorPointResidual( const ROL::Ptr<OBJ> &obj, 
                                   const ROL::Ptr<CON> &con, 
                                   const ROL::Ptr<BND> &bnd,
                                   const V &x,
                                   const ROL::Ptr<const V> &maskL,
                                   const ROL::Ptr<const V> &maskU,
                                         ROL::Ptr<V> &scratch,
                                         Real mu, bool symmetrize ) :
    obj_(obj), con_(con), bnd_(bnd), xl_(bnd->getLowerBound()),
    xu_(bnd->getUpperBound()), maskL_(maskL), maskU_(maskU), scratch_(scratch),
    mu_(mu), symmetrize_(symmetrize), one_(1.0), zero_(0.0), nfval_(0),
    ngrad_(0), ncval_(0) {

    // Get access to the four components
    const PV &x_pv = dynamic_cast<const PV&>(x);
    
    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
    zl_ = x_pv.get(LOWER);
    zu_ = x_pv.get(UPPER); 
  
  }
 

  void update( const Vector<Real> &x, bool flag = true, int iter = -1  ) {

    // Get access to the four components
    const PV &x_pv = dynamic_cast<const PV&>(x);
    
    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
    zl_ = x_pv.get(LOWER);
    zu_ = x_pv.get(UPPER);  

    obj_->update(*x_,flag,iter);
    con_->update(*x_,flag,iter);

  }


  // Evaluate the gradient of the Lagrangian
  void value( V &c, const V &x, Real &tol ) {

    

    Elementwise::Shift<Real> subtract_mu(-mu_);
    Elementwise::Fill<Real>  fill_minus_mu(-mu_);

    const PV &x_pv = dynamic_cast<const PV&>(x);
    PV &c_pv = dynamic_cast<PV&>(c);
  
    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
    zl_ = x_pv.get(LOWER);
    zu_ = x_pv.get(UPPER);  

    ROL::Ptr<V> cx  = c_pv.get(OPT);
    ROL::Ptr<V> cl  = c_pv.get(EQUAL);
    ROL::Ptr<V> czl = c_pv.get(LOWER);
    ROL::Ptr<V> czu = c_pv.get(UPPER);  
    
    /********************************************************************************/
    /* Optimization Components                                                      */
    /********************************************************************************/
    obj_->gradient(*cx,*x_,tol);
    ngrad_++;
   
    con_->applyAdjointJacobian(*scratch_,*l_,*x_,tol);

    cx->plus(*scratch_);

    cx->axpy(-one_,*zl_);
    cx->plus(*zu_);            // cx = g+J'l-zl+zu

    /********************************************************************************/
    /*  Constraint Components                                               */
    /********************************************************************************/

    con_->value(*cl,*x_,tol);     
    ncval_++;

    /********************************************************************************/
    /* Lower Bound Components                                                       */
    /********************************************************************************/
    if( symmetrize_ ) {  // -(x-l)+mu/zl

      czl->applyUnary(fill_minus_mu);
      czl->applyBinary(div_,*zl_);

      scratch_->set(*x_);
      scratch_->axpy(-1.0,*xl_);
      
      czl->plus(*scratch_);  
      czl->scale(-1.0);
    }
    else {  // czl = zl*(x-l)-mu*e
      czl->set(*x_);                   // czl = x
      czl->axpy(-1.0,*xl_);            // czl = x-l
      czl->applyBinary(mult_,*zl_);    // czl = zl*(x-l)
      czl->applyUnary(subtract_mu);    // czl = zl*(x-l)-mu*e
    }

    // Zero out elements corresponding to infinite lower bounds
    czl->applyBinary(mult_,*maskL_);

    /********************************************************************************/
    /* Upper Bound Components                                                       */
    /********************************************************************************/
    if( symmetrize_ ) { // -(u-x)+mu/zu
     
      czu->applyUnary(fill_minus_mu);
      czu->applyBinary(div_,*zu_);

      scratch_->set(*xu_);
      scratch_->axpy(-1.0, *x_);

      czu->plus(*scratch_);
      czu->scale(-1.0);
    }
    else { // zu*(u-x)-mu*e
      czu->set(*xu_);                  // czu = u
      czu->axpy(-1.0,*x_);             // czu = u-x
      czu->applyBinary(mult_,*zu_);    // czu = zu*(u-x)
      czu->applyUnary(subtract_mu);    // czu = zu*(u-x)-mu*e
    }

    // Zero out elements corresponding to infinite upper bounds
    czu->applyBinary(mult_,*maskU_);

  }

  // Evaluate the action of the Hessian of the Lagrangian on a vector
  //
  // [ J11 J12 J13 J14 ] [ vx  ]   [ jvx  ]   [ J11*vx + J12*vl + J13*vzl + J14*vzu ]
  // [ J21  0   0   0  ] [ vl  ] = [ jvl  ] = [ J21*vx                              ]
  // [ J31  0  J33  0  ] [ vzl ]   [ jvzl ]   [ J31*vx          + J33*vzl           ]
  // [ J41  0   0  J44 ] [ vzu ]   [ jvzu ]   [ J41*vx                    + J44*vzu ]
  //
  void applyJacobian( V &jv, const V &v, const V &x, Real &tol ) {

    

    PV &jv_pv = dynamic_cast<PV&>(jv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
    const PV &x_pv = dynamic_cast<const PV&>(x); 

    // output vector components
    ROL::Ptr<V> jvx  = jv_pv.get(OPT);
    ROL::Ptr<V> jvl  = jv_pv.get(EQUAL);
    ROL::Ptr<V> jvzl = jv_pv.get(LOWER);
    ROL::Ptr<V> jvzu = jv_pv.get(UPPER);

    // input vector components
    ROL::Ptr<const V> vx  = v_pv.get(OPT);
    ROL::Ptr<const V> vl  = v_pv.get(EQUAL); 
    ROL::Ptr<const V> vzl = v_pv.get(LOWER);
    ROL::Ptr<const V> vzu = v_pv.get(UPPER);

    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
    zl_ = x_pv.get(LOWER);
    zu_ = x_pv.get(UPPER);  


    /********************************************************************************/
    /* Optimization Components                                                      */
    /********************************************************************************/

    obj_->hessVec(*jvx,*vx,*x_,tol);
    con_->applyAdjointHessian(*scratch_,*l_,*vx,*x_,tol);
    jvx->plus(*scratch_); 
    con_->applyAdjointJacobian(*scratch_,*vl,*x_,tol);
    jvx->plus(*scratch_);

    // H_13 = -I for l_i > -infty
    scratch_->set(*vzl);
    scratch_->applyBinary(mult_,*maskL_);
    jvx->axpy(-1.0,*scratch_);
    
    // H_14 = I for u_i < infty
    scratch_->set(*vzu);
    scratch_->applyBinary(mult_,*maskU_);
    jvx->plus(*scratch_);

    /********************************************************************************/
    /*  Constraint Components                                               */
    /********************************************************************************/

    con_->applyJacobian(*jvl,*vx,*x_,tol);

    /********************************************************************************/
    /* Lower Bound Components                                                       */
    /********************************************************************************/

    if( symmetrize_ ) {
      // czl  = x-l-mu/zl
      // jvzl = -vx - inv(Zl)*(X-L)*vzl
      
      jvzl->set(*x_);
      jvzl->axpy(-1.0,*xl_);
      jvzl->applyBinary(mult_,*vzl);
      jvzl->applyBinary(div_,*zl_);
          
      jvzl->plus(*vx);
      jvzl->scale(-1.0);
      
    }

    else {
      // czl  = zl*(x-l)-mu*e 
      // jvzl = Zl*vx + (X-L)*vzl 

      // H_31 = Zl
      jvzl->set(*vx);
      jvzl->applyBinary(mult_,*zl_);

      // H_33 = X-L
      scratch_->set(*x_);
      scratch_->axpy(-1.0,*xl_);
      scratch_->applyBinary(mult_,*vzl);

      jvzl->plus(*scratch_);

    }

    // jvzl[i] = vzl[i] if l[i] = -inf
    jvzl->applyBinary(mult_,*maskL_);
    jvzl->applyBinary(inFill_,*vzl);

    /********************************************************************************/
    /* Upper Bound Components                                                       */
    /********************************************************************************/
  
    if( symmetrize_ ) {
      // czu  = u-x-mu/zu
      // jvzu = vx - inv(Zu)*(U-X)*vzu
     
      jvzu->set(*xu_);
      jvzu->axpy(-1.0,*x_);
      jvzu->applyBinary(mult_,*vzu);
      jvzu->applyBinary(div_,*zu_);
      jvzu->scale(-1.0); 
      jvzu->plus(*vx);

    }
    else {
      // czu = zu*(u-x)-mu*e
      // jvzu = -Zu*vx + (U-X)*vzu 
   
      // H_41 = -Zu
      scratch_->set(*zu_);
      scratch_->applyBinary(mult_,*vx);

      // H_44 = U-X
      jvzu->set(*xu_);
      jvzu->axpy(-1.0,*x_);
      jvzu->applyBinary(mult_,*vzu);
 
      jvzu->axpy(-1.0,*scratch_);   
  
    }

    // jvzu[i] = vzu[i] if u[i] = inf 
    jvzu->applyBinary(mult_,*maskU_);        
    jvzu->applyBinary(inFill_,*vzu);         

  }

  // Call this whenever mu changes
  void reset( const Real mu ) {
    mu_ = mu;
    nfval_ = 0;
    ngrad_ = 0;
    ncval_ = 0;
  }
  
  int getNumberFunctionEvaluations(void) const {
    return nfval_;
  }

  int getNumberGradientEvaluations(void) const {
    return ngrad_;
  }

  int getNumberConstraintEvaluations(void) const {
    return ncval_;
  }

}; // class PrimalDualInteriorPointResidual

} // namespace ROL

#endif // ROL_PRIMALDUALINTERIORPOINTRESIDUAL_H
