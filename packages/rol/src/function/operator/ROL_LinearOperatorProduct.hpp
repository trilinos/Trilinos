// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEAROPERATORPRODUCT_H
#define ROL_LINEAROPERATORPRODUCT_H

#include "ROL_LinearOperator.hpp"

/** @ingroup func_group
    \class ROL::LinearOperatorProduct
    \brief Provides the interface to the sequential application of linear operators.
    ---
*/


namespace ROL {

template <class Real>
class LinearOperatorProduct : public LinearOperator<Real> {

  typedef Vector<Real>         V;
  typedef LinearOperator<Real> OP;
 
  typedef typename std::vector<ROL::Ptr<OP> >::size_type size_type;

private:

  ROL::Ptr<std::vector<ROL::Ptr<OP> > > ops_;

public:

  LinearOperatorSum( ROL::Ptr<OP> &A, 
                     ROL::Ptr<OP> &B) {
    ops_ = ROL::makePtr<std::vector<OP> >>();
    ops_->push_back(A);
    ops_->push_back(B);
  }

  LinearOperatorSum( ROL::Ptr<OP> &A, 
                     ROL::Ptr<OP> &B, 
                     ROL::Ptr<OP> &C) {
    ops_ = ROL::makePtr<std::vector<OP> >>();
    ops_->push_back(A);
    ops_->push_back(B);
    ops_->push_back(C);
  }

  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    for( size_type i=0; i<ops_->size(); ++i ) {
      (*ops_)[i]->update(x,flag,true);
    }
  }

  virtual void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const { 
     Hv.set(v);   
     for( size_type i=0; i<ops_->size(); ++i ) {
      (*ops_)[i]->apply(Hv,Hv,tol);
    }
  }

  virtual void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
    Hv.set(v);
     for( size_type i=ops_->size()-1; i>=0; --i ) {
      (*ops_)[i]->applyInverse(Hv,Hv,tol);
    }
  }

}; // class LinearOperatorProduct

} // namespace ROL

#endif // ROL_LINEAROPERATORPRODUCT_H
