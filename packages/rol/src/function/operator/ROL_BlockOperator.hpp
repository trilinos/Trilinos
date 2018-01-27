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

    TEUCHOS_TEST_FOR_EXCEPTION( (nvec1 != nvec2), std::invalid_argument,
                                ">>> ERROR (ROL_BlockOperator, apply): "
                                "Mismatch between input and output number of subvectors.");

    TEUCHOS_TEST_FOR_EXCEPTION( (nblks != nvec1*nvec2 ) , std::invalid_argument, 
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
