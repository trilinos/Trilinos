//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_HELPERS_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_HELPERS_HPP

#include "Rythmos_Types.hpp"

#include "Rythmos_RKButcherTableauBase.hpp"
#include "Teuchos_Assert.hpp"
#include "Thyra_ProductVectorBase.hpp"

namespace Rythmos {

/* \brief . */
template<class Scalar>
void assembleIRKState(
  const int stageIndex,
  const Teuchos::SerialDenseMatrix<int,Scalar> &A_in,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages_in = A_in.numRows();
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( stageIndex, 0, numStages_in );
  TEUCHOS_ASSERT_EQUALITY( A_in.numRows(), numStages_in );
  TEUCHOS_ASSERT_EQUALITY( A_in.numCols(), numStages_in );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages_in );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages_in; ++j ) {
    Vp_StV( outArg(x_out), dt * A_in(stageIndex,j), *x_stage_bar.getVectorBlock(j) );
  }

}


/* \brief . */
template<class Scalar>
void assembleIRKSolution(
  const Teuchos::SerialDenseVector<int,Scalar> &b_in,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages_in = b_in.length();
  TEUCHOS_ASSERT_EQUALITY( b_in.length(), numStages_in );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages_in );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages_in; ++j ) {
    Vp_StV( outArg(x_out), dt * b_in(j), *x_stage_bar.getVectorBlock(j) );
  }

}

/* \brief . */
template<class Scalar>
void assembleERKState(
  const int stageIndex,
  const Teuchos::SerialDenseMatrix<int,Scalar> &A_in,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::VectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{
  TEST_FOR_EXCEPT(true);
}

/* \brief . */
template<class Scalar>
void assembleERKSolution(
  const Teuchos::SerialDenseVector<int,Scalar> &b_in,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::VectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
bool isEmptyRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt ) {
  typedef ScalarTraits<Scalar> ST;
  
  // Check that numStages > 0
  if (rkbt.numStages() == 0) {
    return true;
  }

  // Check that the b vector has _some_ non-zero entry
  int numNonZero = 0;
  int numStages_local = rkbt.numStages();
  const Teuchos::SerialDenseVector<int,Scalar> b_local = rkbt.b();
  for (int i=0 ; i<numStages_local ; ++i) {
    if (b_local(i) != ST::zero()) {
      numNonZero++;
    }
  }
  if (numNonZero == 0) {
    return true;
  }
  // There is no reason to check A and c because they can be zero and you're
  // producing an explicit method as long as b has something in it.
  return false;
}


template<class Scalar>
void assertNonEmptyRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  TEST_FOR_EXCEPTION( isEmptyRKButcherTableau(rkbt), std::logic_error,
      "Error, this RKButcherTableau is either empty or the b vector is all zeros!\n"
      );
}

template<class Scalar>
bool isDIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  if (isEmptyRKButcherTableau(rkbt)) {
    return false;
  }
  typedef ScalarTraits<Scalar> ST;
  int numStages_local = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A_local = rkbt.A();
  for (int i=0 ; i<numStages_local ; ++i) {
    for (int j=0 ; j<numStages_local ; ++j) {
      if ((j>i) && (A_local(i,j) != ST::zero())) {
        return false;
      }
    }
  }
  return true;
}

template<class Scalar>
bool isIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt ) 
{
  if (isEmptyRKButcherTableau(rkbt)) {
    return false;
  }
  return true;
}

template<class Scalar>
void validateIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  TEST_FOR_EXCEPTION( !isIRKButcherTableau(rkbt), std::logic_error,
    "Error!  This implicit RK Butcher Tableau is empty!\n"
    );
}

template<class Scalar>
void validateDIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  TEST_FOR_EXCEPTION( !isDIRKButcherTableau(rkbt), std::logic_error,
      "Error!  This Diagonal Implicit RK Butcher Tableau has non-zeros in the upper triangular part!\n" 
      );
}

template<class Scalar>
bool isSDIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt ) 
{
  if (isEmptyRKButcherTableau(rkbt)) {
    return false;
  }
  if (!isDIRKButcherTableau(rkbt)) {
    return false;
  }
  // Verify the diagonal entries are all equal.
  typedef ScalarTraits<Scalar> ST;
  int numStages_local = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A_local = rkbt.A();
  Scalar val = A_local(0,0);
  for (int i=0 ; i<numStages_local ; ++i) {
    if (A_local(i,i) != val) {
      return false;
    }
  }
  return true;
}

template<class Scalar>
void validateSDIRKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  TEST_FOR_EXCEPTION( !isSDIRKButcherTableau(rkbt), std::logic_error,
      "Error!  This Singly Diagonal Implicit RK Butcher Tableau does not have equal diagonal entries!\n"
      );
}

template<class Scalar>
bool isERKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt) 
{
  if (isEmptyRKButcherTableau(rkbt)) {
    return false;
  }
  // Verify the diagonal is zero and the upper triangular part is zero
  typedef ScalarTraits<Scalar> ST;
  int numStages_local = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A_local = rkbt.A();
  for (int i=0 ; i<numStages_local ; ++i) {
    for (int j=0 ; j<numStages_local ; ++j) {
      if ((j>=i) && ((A_local(i,j) != ST::zero()))) {
        return false;
      }
    }
  }
  const Teuchos::SerialDenseVector<int,Scalar> c_local = rkbt.c();
  if( c_local(0) != ST::zero() ) {
    return false;
  }
  // 08/13/08 tscoffe:  I'm not sure what else I can check for b & c...
  return true;
}



template<class Scalar>
void validateERKButcherTableau( const RKButcherTableauBase<Scalar>& rkbt )
{
  TEST_FOR_EXCEPTION( !isERKButcherTableau(rkbt), std::logic_error,
      "Error!  This ERK Butcher Tableau is not lower triangular or c(0) is not zero!\n" 
      );
}

/*
template<class Scalar>
void validateERKOrder( RKButcherTableauBase<Scalar> rkbt, int order_in )
{
  typedef ScalarTraits<Scalar> ST;
  Teuchos::SerialDenseMatrix<int,Scalar> A_local = rkbt.A();
  Teuchos::SerialDenseVector<int,Scalar> b_local = rkbt.b();
  Teuchos::SerialDenseVector<int,Scalar> c_local = rkbt.c();
  int N = rkbt.numStages();
  TEST_FOR_EXCEPT(N == 0);

  if (order_in == 3) {
    Scalar sum1 = ST::zero();
    Scalar sum2 = ST::zero();
    Scalar sum3 = ST::zero();
    Scalar sum4 = ST::zero();
    for (int j=0 ; j<N ; ++j) {
      sum1 += b_local(j);
      for (int k=0 ; k<N ; ++k) {
        sum2 += 2*b_local(j)*A_local(j,k);
        for (int l=0 ; l<N ; ++l) {
          sum3 += 3*b_local(j)*A_local(j,k)*A_local(j,l);
          sum4 += 6*b_local(j)*A_local(j,k)*A_local(k,l);
        }
      }
    }
    TEST_FOR_EXCEPTION(
        (
         ( sum1 != ST::one() ) || 
         ( sum2 != ST::one() ) ||
         ( sum3 != ST::one() ) ||
         ( sum4 != ST::one() )
        ), 
        std::logic_error,
        "Error!, this RK Butcher Tableau does not meet the order conditions for 3rd order\n"
        );
  } else {
    TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error!  this function is only defined for order 3\n"
        );
  }
}

template<class Scalar>
void validateIRKOrder( RKButcherTableauBase<Scalar> rkbt, int order_in )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void validateDIRKOrder( RKButcherTableauBase<Scalar> rkbt, int order_in )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void validateSDIRKOrder( RKButcherTableauBase<Scalar> rkbt, int order_in )
{
  TEST_FOR_EXCEPT(true);
}
*/
 
enum E_RKButcherTableauTypes {
  RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_INVALID,
  RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_ERK,
  RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK,
  RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK,
  RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK
};

template<class Scalar>
E_RKButcherTableauTypes determineRKBTType(const RKButcherTableauBase<Scalar>& rkbt) {
  if (isEmptyRKButcherTableau(rkbt)) {
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_INVALID;
  }
  if (isERKButcherTableau(rkbt)) {
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_ERK;
  }
  if (rkbt.numStages() == 1) { 
    // In this case, its just an IRK method.
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK;
  }
  if (isSDIRKButcherTableau(rkbt)) {
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK;
  }
  if (isDIRKButcherTableau(rkbt)) {
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK;
  }
  if (isIRKButcherTableau(rkbt)) {
    return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK;
  }
  return RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_INVALID;
}



} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_HELPERS_HPP
