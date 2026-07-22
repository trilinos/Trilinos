// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLOCKOPERATOR2DETERMINANT_H
#define ROL_BLOCKOPERATOR2DETERMINANT_H

#include "ROL_BlockOperator2.hpp"

/** @ingroup func_group
    \class ROL::BlockOperator2Determinant
    \brief Provides the interface to the block determinant of a 2x2 block operator
    ---
*/


namespace ROL {

template <class Real>
class BlockOperator2Determinant : public LinearOperator<Real> {

  typedef Vector<Real>         V;
  typedef LinearOperator<Real> OP;
 
private:

  ROL::Ptr<OP> A_, B_, C_, D_;
  ROL::Ptr<V> scratch_;  

public:

  BlockOperator2Determinant( ROL::Ptr<OP> &A, 
                             ROL::Ptr<OP> &B, 
                             ROL::Ptr<OP> &C, 
                             ROL::Ptr<OP> &D,
                             ROL::Ptr<V>  &scratch ) : 
    A_(A), B_(B), C_(C), D_(D), scratch_(scratch) {}


  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    A_->update(x,flag,true);
    B_->update(x,flag,true);
    C_->update(x,flag,true);
    D_->update(x,flag,true);
  }

  // Apply the determinant \f$(A-BD^{-1}B)\f$
  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
    B_->apply(*scratch_,v,tol);
    D_->applyInverse(Hv,*scratch_,tol);
    C_->apply(*scratch_,Hv,tol);
    A_->apply(Hv,v,tol);
    Hv.axpy(-1.0,*scratch_);  
  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {

    ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_BlockOperator2Determinant, applyInverse): "
                                "Not implemented."); 
  }

}; // class BlockOperator2Determinant

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2DETERMINANT_H
