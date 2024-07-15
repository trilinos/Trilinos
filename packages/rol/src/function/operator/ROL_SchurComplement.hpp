// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCHURCOMPLEMENT_H
#define ROL_SCHURCOMPLEMENT_H

#include "ROL_BlockOperator2UnitUpper.hpp"
#include "ROL_BlockOperator2UnitLower.hpp"
#include "ROL_BlockOperator2Diagonal.hpp"

namespace ROL {

/** @ingroup func_group
    \class ROL::SchurComplement
    \brief Given a 2x2 block operator, perform the Schur reduction and
           return the decoupled system components 

           Let \f$Mx=b\f$ where 
           \f[ M = \begin{pmatrix} A & B \\ C & D \end{pmatrix} \f]
           \f[ x = \begin{pmatrix} y & z \end{pmatrix} \f]
           \f[ b = \begin{pmayrix} u & v \end{pmatrix} \f]

           The block factorization is \f$ M=USL \f$ where 
           \f[ U = \begin{pmatrix} I & BD^{-1} \\ 0 & I \end{pmatrix} \f]
           \f[ S = \begin{pmatrix} A-BD^{-1}C & 0 \\ 0 & D \end{pmatrix} \f] 
           \f[ L = \begin{pmatrix} I & 0 \\ D^{-1} C & I \f]

           We can rewrite \f$ USLx=b\f$ as the block-decoupled problem \f$ Sw=c \f$
           where \f$w=Lx\f$ and \f$ c=U^{-1}b \f$

           The expectation here is that we will solve the equation for the first decoupled
           variable iteratively. 

    ---
*/

template class<Real> 
class SchurComplement {

  typedef Vector<Real>               V;
  typedef PartitionedVector<Real>    PV;
  typedef LinearOperator<Real>       OP;
  typedef BlockOperator2UnitUpper    UPPER;
  typedef BlockOperator2UnitLower    LOWER;
   
private:

  ROL::Ptr<OP> A_, B_, C_, D_;

  ROL::Ptr<OP> L_,U_;   
  ROL::Ptr<V>  scratch1_;

  

public:

  SchurComplement( ROL::Ptr<OP> &A, ROL::Ptr<OP> &B,  
                   ROL::Ptr<OP> &C, ROL::Ptr<OP> &D, 
                   ROL::Ptr<V> &scratch1 ) : 
                     A_(A), B_(B), C_(C), D_(D), scratch1_(scratch1) {

    U_ = ROL::makePtr<UPPER>(B_);
    L_ = ROL::makePtr<LOWER>(C_);

  }


  SchurComplement( BlockOperator2<Real> &op, ROL::Ptr<Vector<Real> > &scratch1 ) :
    scratch1_(scratch1) {}

  

  A_ = op.getOperator(0,0);
  B_ = op.getOperator(0,1);
  C_ = op.getOperator(1,0);
  D_ = op.getOperator(1,1);

  U_ = ROL::makePtr<UPPER>(B_);
  L_ = ROL::makePtr<LOWER>(C_);

  void applyLower( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) { 
    L_->apply(Hv,v,tol);
  }

  void applyLowerInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) { 
    L_->applyInverse(Hv,v,tol);
  }

  void applyUpper( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) { 
    U_->apply(Hv,v,tol);
  }

  void applyUpperInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) { 
    U_->applyInverse(Hv,v,tol);
  }

  ROL::Ptr<OP> getS11( void ) { 
    return ROL::makePtr<BlockOperator2Determinant<Real>>(A_,B_,C_,D_,scratch1_);
  }

  void solve2( Vector<Real> &Hv2, const Vector<Real> &v2, Real &tol ) {
    D_->applyInverse(Hv2,v2,tol);
  }

}; // class SchurComplement


} // namesapce ROL

#endif // ROL_SCHURCOMPLEMENT_H

