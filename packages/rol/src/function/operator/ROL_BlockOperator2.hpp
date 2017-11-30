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

#ifndef ROL_BLOCKOPERATOR2_H
#define ROL_BLOCKOPERATOR2_H

#include "ROL_BlockOperator.hpp"

/** @ingroup func_group
    \class ROL::BlockOperator2
    \brief Provides the interface to apply a 2x2 block operator to a 
           partitioned vector.

    ---
*/

namespace ROL {

template<class Real>
class BlockOperator2 : public LinearOperator<Real> {

  typedef Vector<Real>               V;    // ROL Vector
  typedef PartitionedVector<Real>    PV;   // ROL PartitionedVector
  typedef LinearOperator<Real>       OP;   // Linear Operator

private:

  ROL::Ptr<OP> bkop_;
  ROL::Ptr<std::vector<ROL::Ptr<OP> > > ops_;

public:   

  BlockOperator2( ROL::Ptr<OP> &a11, ROL::Ptr<OP> &a12,
                  ROL::Ptr<OP> &a21, ROL::Ptr<OP> &a22 ) {

    using std::vector;
    
    

    ops_ = ROL::makePtr<vector<ROL::Ptr<OP> >>();
    ops_->push_back(a11);
    ops_->push_back(a12);
    ops_->push_back(a21);
    ops_->push_back(a22);
            
    bkop_ = ROL::makePtr<BlockOperator<Real>>(ops_);
 
  }


  void apply( V &Hv, const V &v, Real &tol ) const {
    bkop_->apply(Hv,v,tol);  
  }


  void applyInverse( V &Hv, const V &v, Real &tol ) const {

    TEUCHOS_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_BlockOperator2, applyInverse): "
                                "Not implemented."); 

  } 

  ROL::Ptr<LinearOperator<Real> > getOperator( int row, int col ) const {
    int dex = 2*row+col;
    if( 0<=dex && dex<=3 ) {
      return (*ops_)[dex];
    } 
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_BlockOperator2, getOperator): "
                                  "invalid block indices."); 
    }
    
  }


}; // class BlockOperator2

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2_H
