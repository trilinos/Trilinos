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

#ifndef ROL_BLOCKOPERATOR2UNITLOWER_H
#define ROL_BLOCKOPERATOR2UNITLOWER_H

#include "ROL_BlockOperator2.hpp"

/** @ingroup func_group
    \class ROL::BlockOperator2UnitLower
    \brief Provides the interface to apply a 2x2 block unit lower operator to a 
           partitioned vector.

           \f[ U = \begin{pmatrix} I & 0 \\ C & I \end{pmatrix} \f]
           \f[ U^{-1} = \begin{pmatrix} I & 0 \\ -C & I \end{pmatrix} \f]

    ---
*/

namespace ROL {

template<class Real>
class BlockOperator2UnitLower : public LinearOperator<Real> {

  typedef Vector<Real>               V;    // ROL Vector
  typedef PartitionedVector<Real>    PV;   // ROL PartitionedVector
  typedef LinearOperator<Real>       OP;   // Linear Operator

private:

  ROL::Ptr<OP> C_;


public:   

  BlockOperator2UnitLower( ROL::Ptr<OP> &C ) : C_(C) {}
  }

  void apply( V &Hv, const V &v, Real &tol ) const {
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
      
    ROL::Ptr<V> Hv1 = Hv_pv.get(0);
    ROL::Ptr<V> Hv2 = Hv_pv.get(1);
    ROL::Ptr<const V> v1 = v_pv.get(0);
    ROL::Ptr<const V> v2 = v_pv.get(1);

    Hv1->set(*v1);
    C_->apply(*Hv2,*v1,tol);
    Hv2->plus(*v2);
  }


  void applyInverse( V &Hv, const V &v Real &tol ) const {
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
      
    ROL::Ptr<V> Hv1 = Hv_pv.get(0);
    ROL::Ptr<V> Hv2 = Hv_pv.get(1);
    ROL::Ptr<const V> v1 = v_pv.get(0);
    ROL::Ptr<const V> v2 = v_pv.get(1);
 
    Hv1->set(*v1);
    C_->apply(*Hv2,*v1,tol);
    Hv2->scale(-1.0); 
    Hv2->plus(*v2);

  } 

  ROL::Ptr<LinearOperator<Real> > getOperator( int row, int col ) const {
    if( row == 1 && col == 0 ) {
      return C_;
    } 
    else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_BlockOperator2UnitLower, getOperator): "
                                  "invalid block indices."); 
    }
    
  }


}; // class BlockOperator2UnitLower

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2UNITLOWER_H
