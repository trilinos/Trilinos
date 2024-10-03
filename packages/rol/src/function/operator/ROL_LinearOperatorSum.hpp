// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAROPERATORSUM_H
#define ROL_LINEAROPERATORSUM_H

#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::LinearOperatorSum
    \brief Provides the interface to sum of linear operators applied to a vector
    ---
*/


namespace ROL {

template <class Real>
class LinearOperatorSum : public LinearOperator<Real>  {

  typedef Vector<Real>         V;
  typedef LinearOperator<Real> OP;
 
  typedef typename std::vector<ROL::Ptr<OP> >::size_type size_type;

private:

  ROL::Ptr<std::vector<ROL::Ptr<OP> > > ops_;
  ROL::Ptr<V> scratch_;

public:

  LinearOperatorSum( ROL::Ptr<OP> &A, 
                     ROL::Ptr<OP> &B, 
                     ROL::Ptr<V> & scratch ) :
    scratch_(scratch) {
    ops_ = ROL::makePtr<std::vector<OP> >>();
    ops_->push_back(A);
    ops_->push_back(B);
  }

  LinearOperatorSum( ROL::Ptr<OP> &A, 
                     ROL::Ptr<OP> &B, 
                     ROL::Ptr<OP> &C, 
                     ROL::Ptr<V> & scratch ) :
    scratch_(scratch) {
    ops_ = ROL::makePtr<std::vector<OP> >>();
    ops_->push_back(A);
    ops_->push_back(B);
    ops_->push_back(C);
  }

  // TODO: implementation for arbitrary sum

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    for( size_type i=0; i<ops_->size(); ++i ) {
      (*ops_)[i]->update(x,flag,true);
    }
  }

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    (*ops_)[0]->apply(Hv,v,tol);
    for( size_type i=1; i<ops_->size(); ++i ) {
      (*ops_)[i]->apply(*scratch_,v,tol);
      Hv.plus(*scratch_);
    }
  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_LinearOperatorSum, applyInverse): "
                                  "Inverse is not defined for general sum of operators.");     
  }

}; // class LinearOperatorSum

} // namespace ROL

#endif // ROL_LINEAROPERATOR_PRODUCT_H
