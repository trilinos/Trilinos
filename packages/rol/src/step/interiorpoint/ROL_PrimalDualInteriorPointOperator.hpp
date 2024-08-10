// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PRIMALDUALINTERIORPOINTOPERATOR_H
#define ROL_PRIMALDUALINTERIORPOINTOPERATOR_H

#include "ROL_LinearOperator.hpp"


namespace ROL { 



template<class Real> 
class PrimalDualInteriorPointBlock11 : public LinearOperator<Real> {

  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;
  typedef Objective<Real>          OBJ;
  typedef Constraint<Real> CON;      

  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 0;
  static const size_type UPPER = 1;
 
  ROL::Ptr<const V>   x_;            // Optimization vector
  ROL::Ptr<const V>   l_;            //  constraint multiplier

  ROL::Ptr<V> scratch_;

  Real delta_; // Initial correction


public: 

  PrimalDualInteriorPointBlock11( ROL::Ptr<OBJ> &obj, ROL::Ptr<CON> &con, 
                           const V &x, ROL::Ptr<V> & scratch, 
                           Real delta=0 ) : 
    obj_(obj), con_(con), scratch_(scratch), delta_(delta) { 

    const PV &x_pv = dynamic_cast<const PV&>(x);

    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {

    const PV &x_pv = dynamic_cast<const PV&>(x);
    
    x_  = x_pv.get(OPT);
    l_  = x_pv.get(EQUAL);

    obj_->update(*x_,flag,true);
    con_->update(*x_,flag,true);  
  }
 
  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {

    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);

    // output vector components
    ROL::Ptr<V> Hvx  = Hv_pv.get(OPT);
    ROL::Ptr<V> Hvl  = Hv_pv.get(EQUAL);

    // input vector components
    ROL::Ptr<const V> vx  = v_pv.get(OPT);
    ROL::Ptr<const V> vl  = v_pv.get(EQUAL); 

    obj_->hessVec(*jvx,*vx,*x_,tol);
    con_->applyAdjointHessian(*scratch_,*l_,*vx,*x_,tol);
    jvx->plus(*scratch_); 
    con_->applyAdjointJacobian(*scratch_,*vl,*x_,tol);
    jvx->plus(*scratch_);

    // Inertial correction
    if( delta_ != 0 ) {
      jvx->axpy(delta_,*vx);
    }

    con_->applyJacobian(*jvl,*vx,*x_,tol);

  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
     ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_PrimalDualInteriorPointBlock11, applyInverse): "
                                "Not implemented.");    
  }
  
  void setInertia( Real delta ) { 
    delta_ = delta;
  }


};  // class PrimalDualInteriorPointBlock11


template<class Real> 
class PrimalDualInteriorPointBlock12 : public LinearOperator<Real> {

  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;

  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 0;
  static const size_type UPPER = 1;
 
public:
 
  void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);

    // output vector components
    ROL::Ptr<V> Hvx  = Hv_pv.get(OPT);
    ROL::Ptr<V> Hvl  = Hv_pv.get(EQUAL);

    // input vector components
    ROL::Ptr<const V> vzl  = v_pv.get(LOWER);
    ROL::Ptr<const V> vzu  = v_pv.get(UPPER); 

    Hvx->set(*vzu);
    Hvx->axpy(-1.0,*vzl);
    Hvl->zero();
   
  }

  void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
     ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_PrimalDualInteriorPointBlock12, applyInverse): "
                                "Not implemented.");    
  }


};  // class PrimalDualInteriorPointBlock12


template<class Real> 
class PrimalDualInteriorPointBlock21 : public LinearOperator<Real> {

  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;

  static const size_type OPT   = 0;
  static const size_type EQUAL = 1;
  static const size_type LOWER = 0;
  static const size_type UPPER = 1;
 
  ROL::Ptr<const V> zl_;
  ROL::Ptr<const V> zu_;

public:
 
  PrimalDualInteriorPointBlock21( const V &z ) {
    const PV &z_pv = dynamic_cast<const PV&>(z);
    zl_ = z_pv.get(LOWER);
    zu_ = z_pv.get(UPPER);  
  }

  void update( const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    const PV &z_pv = dynamic_cast<const PV&>(z);
    zl_ = z_pv.get(LOWER);
    zu_ = z_pv.get(UPPER);  
  } 

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);

    // output vector components
    ROL::Ptr<V> Hvzl  = Hv_pv.get(LOWER);
    ROL::Ptr<V> Hvzu  = Hv_pv.get(UPPER);

    // input vector components
    ROL::Ptr<const V> vx  = v_pv.get(OPT);
    ROL::Ptr<const V> vl  = v_pv.get(EQUAL); 

    Hvzl->set(*vx);
    Hvzl->applyBinary(mult_,*zl_);

    Hvzu->set(*vx);
    Hvzu->applyBinary(mult_,*zu_);
    Hvzu->scale(-1.0); 
  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
     ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_PrimalDualInteriorPointBlock21, applyInverse): "
                                "Not implemented.");    
  }

};  // class PrimalDualInteriorPointBlock21


template<class Real> 
class PrimalDualInteriorPointBlock22 : public LinearOperator<Real> {

  typedef Vector<Real>             V;
  typedef PartitionedVector<Real>  PV;
  typedef BoundConstraint<Real>    BND;      

  static const size_type OPT   = 0;
  static const size_type LOWER = 0;
  static const size_type UPPER = 1;
 
  ROL::Ptr<const V> x_;
  ROL::Ptr<const V> xl_;
  ROL::Ptr<const V> xu_;

  Elementwise::Multiply<Real> mult_;
  Elementwise::Multiply<Real> divinv_;


public:
 
  PrimalDualInteriorPointBlock22( const ROL::Ptr<BND> &bnd, const Vector<Real> &x ) {

    const PV &x_pv = dynamic_cast<const PV&>(x);
    
    x_  = x_pv.get(OPT);
    xl_ = bnd.getLowerBound();
    xu_ = bnd.getUpperBound();

  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {

    const PV &x_pv = dynamic_cast<const PV&>(x);
    x_  = x_pv.get(OPT);
  }
 
  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);

    // output vector components
    ROL::Ptr<V> Hvzl  = Hv_pv.get(LOWER);
    ROL::Ptr<V> Hvzu  = Hv_pv.get(UPPER);

    // input vector components
    ROL::Ptr<const V> vzl = v_pv.get(LOWER);
    ROL::Ptr<const V> vzu = v_pv.get(UPPER); 

    Hvzl->set(*x_);
    Hvzl->axpy(-1.0,*xl_);
    Hvzl->applyBinary(mult_,*vzl);

    Hvzu->set(*xu_);
    Hvzu->axpy(-1.0,*x_);
    Hvzu->applyBinary(mult_,*vzu);

  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {

    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);

    // output vector components
    ROL::Ptr<V> Hvzl  = Hv_pv.get(LOWER);
    ROL::Ptr<V> Hvzu  = Hv_pv.get(UPPER);

    // input vector components
    ROL::Ptr<const V> vzl = v_pv.get(LOWER);
    ROL::Ptr<const V> vzu = v_pv.get(UPPER); 

    Hvzl->set(*x_);
    Hvzl->axpy(-1.0,*xl_);
    Hvzl->applyBinary(divinv_,*vzl); // Hvzl = vzl/(x-xl)

    Hvzu->set(*xu_);
    Hvzu->axpy(-1.0,*x_);
    Hvzu->applyBinary(divinv_,*vzu); // Hvzu = vzu/(xu-x)
  }

};  // class PrimalDualInteriorPointBlock22


} // namespace ROL



#endif  // ROL_PRIMALDUALINTERIORPOINTOPERATOR_H

