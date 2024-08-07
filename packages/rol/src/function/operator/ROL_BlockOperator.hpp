// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLOCKOPERATOR_H
#define ROL_BLOCKOPERATOR_H

#include "ROL_LinearOperator.hpp"
#include "ROL_PartitionedVector.hpp"

/** @ingroup func_group
    \class ROL::BlockOperator
    \brief Provides the interface to apply a block operator to a 
           partitioned vector.

    ---
*/

namespace ROL {

template<class Real>
class BlockOperator : public LinearOperator<Real> {

  typedef Vector<Real>               V;    // ROL Vector
  typedef PartitionedVector<Real>    PV;   // Partitioned Vector
  typedef LinearOperator<Real>       OP;   // Linear Operator

  typedef std::vector<ROL::Ptr<OP> > OpVec;  // Vector (column-stacked matrix) of pointers to operators
  typedef typename OpVec::size_type      uint;   // index type

private:
  
  ROL::Ptr<OpVec> blocks_;

public:
  BlockOperator() {}
  BlockOperator( const ROL::Ptr<OpVec> &blocks ) : blocks_(blocks) {}
  
  virtual void apply( V &Hv, const V &v, Real &tol ) const {
 
    // Downcast to Partitioned Vectors
    PV &Hv_part = dynamic_cast<PV&>(Hv);
    const PV &v_part = dynamic_cast<const PV&>(v);
    
    uint nvec1 = v_part.numVectors();
    uint nvec2 = Hv_part.numVectors();
    uint nblks = blocks_->size();

    ROL_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                                ">>> ERROR (ROL_BlockOperator, apply): "
                                "Mismatch between input and output number of subvectors.");

    ROL_TEST_FOR_EXCEPTION( (nblks != nvec1*nvec2 ) , std::invalid_argument, 
                                ">>> ERROR (ROL_BlockOperator, apply): "
                                "Block operator dimension mismatch."); 

    for( uint i=0; i<nvec1; ++i ) {
      
      ROL::Ptr<V> Hvi = Hv_part.get(i);
      ROL::Ptr<V> u = Hvi->clone();

      u->zero(); 

      for( uint j=0; j<nvec2; ++j ) {
        uint k = i+nvec1*j;
          (*blocks_)[k]->apply(*u,*v_part.get(j),tol);
          Hvi->plus(*u);
      }
    }
  }  


}; // class BlockOperator

} // namespace ROL

#endif // ROL_BLOCKOPERATOR_H
