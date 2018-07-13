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
