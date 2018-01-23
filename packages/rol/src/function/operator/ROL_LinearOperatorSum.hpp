// @HEADER
// ************************************************************************
//
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
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_LinearOperatorSum, applyInverse): "
                                  "Inverse is not defined for general sum of operators.");     
  }

}; // class LinearOperatorSum

} // namespace ROL

#endif // ROL_LINEAROPERATOR_PRODUCT_H
