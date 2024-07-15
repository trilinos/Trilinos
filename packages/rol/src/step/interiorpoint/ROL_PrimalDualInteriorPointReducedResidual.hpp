// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PRIMALDUALINTERIORPOINTREDUCEDRESIDUAL_H
#define ROL_PRIMALDUALINTERIORPOINTREDUCEDRESIDUAL_H

#include "ROL_InteriorPointPenalty.hpp"
#include "ROL_Constraint.hpp"

#include <iostream>

/** @ingroup func_group
    \class   ROL::PrimalDualInteriorPointReducedResidual
    \brief   Reduced form of the Primal Dual Interior Point
             residual and the action of its Jacobian

             The orginal system has the form

             \f[
             \begin{pmatrix}
              W   & J^\ast & -I   & I   \\
              J   & 0      &  0   & 0   \\
              Z_l & 0      &  X-L & 0   \\
             -Z_u & 0      &  0   & U-X
             \end{pmatrix}
             \begin{pmatrix}
             d_x \\ d_\lambda \\ d_{z_l} \\ d_{z_u} 
             \end{pmatrix}
             =
             -\begin{pmatrix}
             \nabla f+J^\ast \lambda -z_l + z_u \\
             c \\
             (x-l)z_l - \mu e \\
             (u-x)z_u - \mu e
             \end{pmatrix} = -\begin{pmatrix} g_x \\ g_\lambda \\ g_{z_l} \\ g_{z_u}
             \f]

             Using the last two equations, we have

             \f[ d_{z_l} = -(X-L)^{-1} g_{z_l} - (X-L)^{-1}Z_l d_x \f]
             \f[ d_{z_u} = -(U-X)^{-1} g_{z_u} + (U-X)^{-1}Z_u d_x \f] 
 
             Substituting into the first equation, we get 

             \f[ [W+(X-L)^{-1}Z_l+(U-X)^{-1}Z_u]d_x + J^\ast d_\lambda = 
                  -g_x + (U-X)^{-1}g_{z_u} - (L-X)^{-1}g_{z_l} \f] 

             This leads to the reduced system

             \f[
             \begin{pmatrix}
             W+\Sigma_l+\Sigma_u & J^\ast \\
             J                   & 0
             \end{pmatrix} 
             \begin{pmatrix}
             d_x \\ d_\lambda
             \end{pmatrix} 
             = -
             \begin{pmatrix} 
             \nabla \varphi_{\mu}(x) + J^\ast \lambda \\ c 
             \end{pmatrix} 
             \f]
             Where \f$\varphi_\mu(x)\f$ is the barrier penalty objective
    
    --- 
*/

namespace ROL { 

template<class Real> 
class PrimalDualInteriorPointResidual : public Constraint<Real> {

  typedef ROL::ParameterList     PL;

  typedef Vector<Real>               V;
  typedef PartitionedVector<Real>    PV;
  typedef Objective<Real>            OBJ;
  typedef Constraint<Real>   CON;
  typedef BoundConstraint<Real>      BND;
  typedef LinearOperator<Real>       LOP;
  typedef InteriorPointPenalty<Real> PENALTY;

  typedef Elementwise::ValueSet<Real>   ValueSet;

  typedef typename PV::size_type size_type;

private: 

  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 2;
  static const size_type UPPER = 3;

  ROL::Ptr<const V>   x_;            // Optimization vector
  ROL::Ptr<const V>   l_;            //  constraint multiplier
  ROL::Ptr<const V>   zl_;           // Lower bound multiplier
  ROL::Ptr<const V>   zu_;           // Upper bound multiplier

  ROL::Ptr<const V>   xl_;           // Lower bound
  ROL::Ptr<const V>   xu_;           // Upper bound 

  const ROL::Ptr<const V> maskL_;   
  const ROL::Ptr<const V> maskU_;

  ROL::Ptr<V> scratch_;

  const ROL::Ptr<PENALTY> penalty_;
  const ROL::Ptr<OBJ>     obj_;
  const ROL::Ptr<CON>     con_;


public:

  PrimalDualInteriorPointResidual( const ROL::Ptr<PENALTY> &penalty, 
                                   const ROL::Ptr<CON> &con,
                                   const V &x,
                                         ROL::Ptr<V> &scratch ) :
    penalty_(penalty), con_(con), scratch_(scratch) {

    obj_   = penalty_->getObjective();
    maskL_ = penalty_->getLowerMask();
    maskU_ = penalty_->getUpperMask();

    ROL::Ptr<BND> bnd = penalty_->getBoundConstraint();
    xl_ = bnd->getLowerBound();
    xu_ = bnd->getUpperBound();


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
    

    PV &c_pv = dynamic_cast<PV&>(c);
    const PV &x_pv = dynamic_cast<const PV&>(x); 

    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
    zl_ = x_pv.get(LOWER);
    zu_ = x_pv.get(UPPER); 

    ROL::Ptr<V> cx  = c_pv.get(OPT);
    ROL::Ptr<V> cl  = c_pv.get(EQUAL);
   
    // TODO: Add check as to whether we really need to recompute these
    penalty_->gradient(*cx,*x_,tol);
    con_->applyAdjointJacobian(*scratch_,*l_,*x_,tol);
   
    // \f$ c_x = \nabla \phi_mu(x) + J^\ast \lambda \f$
    cx_->plus(*scratch_);      
    
    con_->value(*cl_,*x_,tol);
    
  }

  void applyJacobian( V &jv, const V &v, const V &x, Real &tol ) {

    

    PV &jv_pv = dynamic_cast<PV&>(jv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
    const PV &x_pv = dynamic_cast<const PV&>(x); 

    // output vector components
    ROL::Ptr<V> jvx  = jv_pv.get(OPT);
    ROL::Ptr<V> jvl  = jv_pv.get(EQUAL);

    // input vector components
    ROL::Ptr<const V> vx  = v_pv.get(OPT);
    ROL::Ptr<const V> vl  = v_pv.get(EQUAL); 

    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
 
    // \f$ 
    obj_->hessVec(*jvx,*vx,*x_,tol);
    con_->applyAdjointHessian(*scratch_,*l_,*vx,*x_,tol);
    jvx->plus(*scratch_);   

    // \f$J^\ast d_\lambda \f$
    con_->applyAdjointJacobian(*scratch_,*vl,*x_,tol);
    jvx->plus(*scratch_);


    Elementwise::DivideAndInvert<Real> divinv;
    Elementwise::Multiply<Real> mult;

    // Note that indices corresponding to infinite bounds should automatically lead to 
    // zero diagonal contributions to the Sigma operators
    scratch_->set(*x_);
    scratch_->axpy(-1.0,*xl_);           // x-l
    scratch_->applyBinary(divinv,*zl_);  // zl/(x-l)
    scratch_->applyBinary(mult,*vx);     // zl*vx/(x-l) = Sigma_l*vx

    jvx->plu(*scratch_);
    
    scratch_->set(*xu_);
    scratch_->axpy(-1.0,*x_);           // u-x
    scratch_->applyBinary(divinv,*zu_); // zu/(u-x)
    scratch_->applyBinary(mult,*vx);    // zu*vx/(u-x) = Sigma_u*vx

    jvx->plus(*scratch_);  
 
   
     
    



//            \f[ [W+(X-L)^{-1}+(U-X)^{-1}]d_x + J^\ast d_\lambda = 
//                  -g_x + (U-X)^{-1}g_{z_u} - (L-X)^{-1}g_{z_l} \f] 



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

#endif // ROL_PRIMALDUALKKTOPERATOR_H
