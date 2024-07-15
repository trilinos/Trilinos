// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_BlockOperator2UnitLower, getOperator): "
                                  "invalid block indices."); 
    }
    
  }


}; // class BlockOperator2UnitLower

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2UNITLOWER_H
