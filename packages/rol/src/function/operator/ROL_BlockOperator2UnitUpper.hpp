// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLOCKOPERATOR2UNITUPPER_H
#define ROL_BLOCKOPERATOR2UNITUPPER_H

#include "ROL_BlockOperator2.hpp"

/** @ingroup func_group
    \class ROL::BlockOperator2UnitUpper
    \brief Provides the interface to apply a 2x2 block unit upper operator to a 
           partitioned vector.

           \f[ U = \begin{pmatrix} I & B \\ 0 & I \end{pmatrix} \f]
           \f[ U^{-1} = \begin{pmatrix} I & -B \\ 0 & I \end{pmatrix} \f]

    ---
*/

namespace ROL {

template<class Real>
class BlockOperator2UnitUpper : public LinearOperator<Real> {

  typedef Vector<Real>               V;    // ROL Vector
  typedef PartitionedVector<Real>    PV;   // ROL PartitionedVector
  typedef LinearOperator<Real>       OP;   // Linear Operator

private:

  ROL::Ptr<OP> B_;


public:   

  BlockOperator2UnitUpper( ROL::Ptr<OP> &B ) : B_(B) {}
  }

  void apply( V &Hv, const V &v, Real &tol ) const {
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
      
    ROL::Ptr<V> Hv1 = Hv_pv.get(0);
    ROL::Ptr<V> Hv2 = Hv_pv.get(1);
    ROL::Ptr<const V> v1 = v_pv.get(0);
    ROL::Ptr<const V> v2 = v_pv.get(1);

    B_->apply(*Hv1,*v2,tol);
    Hv1->plus(*v1);
    Hv2->set(*v2); 
  }


  void applyInverse( V &Hv, const V &v Real &tol ) const {
    

    PV &Hv_pv = dynamic_cast<PV&>(Hv);
    const PV &v_pv = dynamic_cast<const PV&>(v);
      
    ROL::Ptr<V> Hv1 = Hv_pv.get(0);
    ROL::Ptr<V> Hv2 = Hv_pv.get(1);
    ROL::Ptr<const V> v1 = v_pv.get(0);
    ROL::Ptr<const V> v2 = v_pv.get(1);
 
    B_->apply(*Hv1,*v2,tol);
    Hv1->scale(-1.0);
    Hv1->plus(*v1);
    Hv2->set(*v2); 

  } 

  ROL::Ptr<LinearOperator<Real> > getOperator( int row, int col ) const {
    if( row == 0 && col == 1 ) {
      return B_;
    } 
    else {
      ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument, 
                                  ">>> ERROR (ROL_BlockOperator2UnitUpper, getOperator): "
                                  "invalid block indices."); 
    }
    
  }


}; // class BlockOperator2UnitUpper

} // namespace ROL

#endif // ROL_BLOCKOPERATOR2UNITUPPER_H
