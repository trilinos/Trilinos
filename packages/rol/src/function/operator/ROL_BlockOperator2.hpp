// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_BlockOperator2, applyInverse): "
                                "Not implemented."); 

  } 

  ROL::Ptr<LinearOperator<Real> > getOperator( int row, int col ) const {
    int dex = 2*row+col;
    if( 0<=dex && dex<=3 ) {
      return (*ops_)[dex];
    } 
    else {
      ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_BlockOperator2, getOperator): "
                                  "invalid block indices."); 
    }
    
  }


}; // class BlockOperator2

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2_H
