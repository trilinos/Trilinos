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


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_HPP


#include "Rythmos_Types.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

// 08/13/08 tscoffe:  These are namespaced so they are only availabe in this file.
namespace {
  using Teuchos::RCP;
  using Teuchos::rcp;

  const std::string ForwardEuler_name = "Forward Euler";
  const std::string BackwardEuler_name = "Backward Euler";
  const std::string Explicit4Stage_name = "Explicit 4 Stage";
  const std::string Explicit3_8Rule_name = "Explicit 3/8 Rule";

  const std::string Explicit2Stage2ndOrderRunge_name = "Explicit 2 Stage 2nd order by Runge";
  const std::string Explicit3Stage3rdOrderHeun_name = "Explicit 3 Stage 3rd order by Heun";
  const std::string Explicit3Stage3rdOrder_name = "Explicit 3 Stage 3rd order";
  const std::string Explicit4Stage3rdOrderRunge_name = "Explicit 4 Stage 3rd order by Runge";
  const std::string Explicit4Stage3rdOrder_name = "Explicit 4 stage 4rd order";

  const std::string Implicit1Stage2ndOrderGauss_name = "Implicit 1 Stage 2nd order Gauss";
  const std::string Implicit2Stage4thOrderGauss_name = "Implicit 2 Stage 4th order Gauss";
  const std::string Implicit3Stage6thOrderGauss_name = "Implicit 3 Stage 6th order Gauss";

  const std::string SelectionTypeByName_name = "Method by name";
  const std::string SelectionTypeByName_default = BackwardEuler_name;
  const std::string SelectionTypeExplicitByOrder_name = "Explicit method by order";
  const int SelectionTypeExplicitByOrder_default = 1;
  const std::string SelectionTypeImplicitByOrder_name = "Implicit method by order";
  const int SelectionTypeImplicitByOrder_default = 1;
  const std::string SelectionTypeDIRKByOrder_name = "DIRK method by order";
  const int SelectionTypeDIRKByOrder_default = 1;
  const std::string SelectionTypeSDIRKByOrder_name = "SDIRK method by order";
  const int SelectionTypeSDIRKByOrder_default = 1;

  const std::string SelectionType_name = "Selection Type";
  const std::string SelectionType_default = SelectionTypeByName_name;

  enum E_RKButcherTableauSelectionTypes {
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_INVALID,
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_BY_NAME,
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_EXPLICIT_BY_ORDER,
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_IMPLICIT_BY_ORDER,
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_DIRK_BY_ORDER,
    RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_SDIRK_BY_ORDER
  };
  Teuchos::Array<std::string>
    S_RKButcherTableauSelectionTypes = Teuchos::tuple<std::string>(
        SelectionTypeByName_name,
        SelectionTypeExplicitByOrder_name,
        SelectionTypeImplicitByOrder_name,
        SelectionTypeDIRKByOrder_name,
        SelectionTypeSDIRKByOrder_name
        );
  const RCP<Teuchos::StringToIntegralParameterEntryValidator<E_RKButcherTableauSelectionTypes> >
    rkButcherTableauSelectionTypeValidator = rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<E_RKButcherTableauSelectionTypes>(
          S_RKButcherTableauSelectionTypes,
          Teuchos::tuple<E_RKButcherTableauSelectionTypes>(
            RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_BY_NAME,
            RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_EXPLICIT_BY_ORDER,
            RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_IMPLICIT_BY_ORDER,
            RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_DIRK_BY_ORDER,
            RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_SDIRK_BY_ORDER
            ),
          SelectionType_name
          )
        );


  enum E_RKButcherTableauSelectionMethodNames {
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_INVALID,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_FORWARD_EULER,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_BACKWARD_EULER,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_8_RULE,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_2_STAGE_2_ORDER_RUNGE,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER_HEUN,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER_RUNGE,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_1_STAGE_2_ORDER_GAUSS,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_2_STAGE_4_ORDER_GAUSS,
    RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_3_STAGE_6_ORDER_GAUSS
  };

  Teuchos::Array<std::string>
    S_RKButcherTableauMethodNames = Teuchos::tuple<std::string>(
      ForwardEuler_name,
      BackwardEuler_name,
      Explicit4Stage_name,
      Explicit3_8Rule_name,
      Explicit2Stage2ndOrderRunge_name,
      Explicit3Stage3rdOrderHeun_name,
      Explicit3Stage3rdOrder_name,
      Explicit4Stage3rdOrderRunge_name,
      Explicit4Stage3rdOrder_name,
      Implicit1Stage2ndOrderGauss_name,
      Implicit2Stage4thOrderGauss_name,
      Implicit3Stage6thOrderGauss_name
      );

  const RCP<Teuchos::StringToIntegralParameterEntryValidator<E_RKButcherTableauSelectionMethodNames> >
    rkButcherTableauMethodNameValidator = rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<E_RKButcherTableauSelectionMethodNames>(
          S_RKButcherTableauMethodNames,
          Teuchos::tuple<E_RKButcherTableauSelectionMethodNames>(
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_FORWARD_EULER,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_BACKWARD_EULER,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_8_RULE,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_2_STAGE_2_ORDER_RUNGE,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER_HEUN,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER_RUNGE,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_1_STAGE_2_ORDER_GAUSS,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_2_STAGE_4_ORDER_GAUSS,
            RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_IMPLICIT_3_STAGE_6_ORDER_GAUSS
            ),
          SelectionTypeByName_name
          )
        );

} // namespace


namespace Rythmos {


/* \brief . */
template<class Scalar>
class RKButcherTableau {
public:
  /** \brief. */
  RKButcherTableau()
    {}
  /** \brief. */
  RKButcherTableau(
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c
    )
    : A_(A), b_(b), c_(c)
    {
      const int numStages = A.numRows();
      TEUCHOS_ASSERT_EQUALITY( A.numRows(), numStages );
      TEUCHOS_ASSERT_EQUALITY( A.numCols(), numStages );
      TEUCHOS_ASSERT_EQUALITY( b.length(), numStages );
      TEUCHOS_ASSERT_EQUALITY( c.length(), numStages );
    }
  /** \brief . */
  int numStages() const { return A_.numRows(); }
  /** \brief . */
  const Teuchos::SerialDenseMatrix<int,Scalar>& A() const { return A_; }
  /** \brief . */
  const Teuchos::SerialDenseVector<int,Scalar> b() const { return b_; }
  /** \brief . */
  const Teuchos::SerialDenseVector<int,Scalar> c() const { return c_; }
private:
  Teuchos::SerialDenseMatrix<int,Scalar> A_;
  Teuchos::SerialDenseVector<int,Scalar> b_;
  Teuchos::SerialDenseVector<int,Scalar> c_;
};


/* \brief . */
template<class Scalar>
void assembleIRKState(
  const int stageIndex,
  const Teuchos::SerialDenseMatrix<int,Scalar> &A,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages = A.numRows();
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( stageIndex, 0, numStages );
  TEUCHOS_ASSERT_EQUALITY( A.numRows(), numStages );
  TEUCHOS_ASSERT_EQUALITY( A.numCols(), numStages );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages; ++j ) {
    Vp_StV( outArg(x_out), dt * A(stageIndex,j), *x_stage_bar.getVectorBlock(j) );
  }

}


/* \brief . */
template<class Scalar>
void assembleIRKSolution(
  const Teuchos::SerialDenseVector<int,Scalar> &b,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::ProductVectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{

  typedef ScalarTraits<Scalar> ST;

  const int numStages = b.length();
  TEUCHOS_ASSERT_EQUALITY( b.length(), numStages );
  TEUCHOS_ASSERT_EQUALITY( x_stage_bar.productSpace()->numBlocks(), numStages );
  Thyra::VectorBase<Scalar>& x_out = *x_out_ptr;

  V_V( outArg(x_out), x_base );
  for ( int j = 0; j < numStages; ++j ) {
    Vp_StV( outArg(x_out), dt * b(j), *x_stage_bar.getVectorBlock(j) );
  }

}

/* \brief . */
template<class Scalar>
void assembleERKState(
  const int stageIndex,
  const Teuchos::SerialDenseMatrix<int,Scalar> &A,
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
  const Teuchos::SerialDenseVector<int,Scalar> &b,
  const Scalar dt,
  const Thyra::VectorBase<Scalar> &x_base,
  const Thyra::VectorBase<Scalar> &x_stage_bar,
  Teuchos::Ptr<Thyra::VectorBase<Scalar> > x_out_ptr
  )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void assertNonEmptyRKButcherTableau( RKButcherTableau<Scalar> rkbt )
{
  typedef ScalarTraits<Scalar> ST;
  
  // Check that numStages > 0
  TEST_FOR_EXCEPTION( rkbt.numStages() == 0, std::logic_error,
     "Error!  This RKButcherTableau has no stages!\n"
     );

  // Check that the b vector has _some_ non-zero entry
  int numNonZero = 0;
  int numStages = rkbt.numStages();
  const Teuchos::SerialDenseVector<int,Scalar> b = rkbt.b();
  for (int i=0 ; i<numStages ; ++i) {
    if (b(i) != ST::zero()) {
      numNonZero++;
    }
  }
  TEST_FOR_EXCEPTION( numNonZero == 0, std::logic_error,
      "Error!  This RKButcherTableau's b vector is all zeros!\n"
      );

  // There is no reason to check A and c because they can be zero and you're
  // producing an explicit method as long as b has something in it.
}

template<class Scalar>
void validateIRKButcherTableau( RKButcherTableau<Scalar> rkbt )
{
  // For now, all I can assert is that the tableau is non-empty
  assertNonEmptyRKButcherTableau(rkbt);
}

template<class Scalar>
void validateDIRKButcherTableau( RKButcherTableau<Scalar> rkbt )
{
  assertNonEmptyRKButcherTableau(rkbt);
  // Verify the upper triangular part is zero
  typedef ScalarTraits<Scalar> ST;
  int numStages = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A = rkbt.A();
  for (int i=0 ; i<numStages ; ++i) {
    for (int j=0 ; j<numStages ; ++j) {
      if (j>i) {
        TEST_FOR_EXCEPTION( A(i,j) != ST::zero(), std::logic_error,
           "Error!  This Diagonal Implicit RK Butcher Tableau has non-zeros in the upper triangular part!\n" 
           );
      }
    }
  }
}

template<class Scalar>
void validateSDIRKButcherTableau( RKButcherTableau<Scalar> rkbt )
{
  // For now, all I can assert is that the tableau is non-empty
  assertNonEmptyRKButcherTableau(rkbt);
  validateDIRKButcherTableau( rkbt );
  // Verify the diagonal entries are all equal.
  typedef ScalarTraits<Scalar> ST;
  int numStages = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A = rkbt.A();
  Scalar val = A(0,0);
  for (int i=0 ; i<numStages ; ++i) {
    TEST_FOR_EXCEPTION( A(i,i) != val, std::logic_error,
        "Error!  This Singly Diagonal Implicit RK Butcher Tableau does not have equal diagonal entries!\n"
        );
  }
}


template<class Scalar>
void validateERKButcherTableau( RKButcherTableau<Scalar> rkbt )
{
  assertNonEmptyRKButcherTableau(rkbt);
  // Verify the diagonal is zero and the upper triangular part is zero
  typedef ScalarTraits<Scalar> ST;
  int numStages = rkbt.numStages();
  const Teuchos::SerialDenseMatrix<int,Scalar> A = rkbt.A();
  for (int i=0 ; i<numStages ; ++i) {
    for (int j=0 ; j<numStages ; ++j) {
      if (j>=i) {
        TEST_FOR_EXCEPTION( A(i,j) != ST::zero(), std::logic_error,
           "Error!  This ERK Butcher Tableau is not lower triangular!\n" 
           );
      }
    }
  }
  const Teuchos::SerialDenseVector<int,Scalar> c = rkbt.c();
  TEST_FOR_EXCEPTION( c(0) != ST::zero(), std::logic_error,
      "Error!  c(0) must be zero for an explicit RK method!\n"
      );
  // 08/13/08 tscoffe:  I'm not sure what else I can check for b & c...
}

template<class Scalar>
void validateERKOrder( RKButcherTableau<Scalar> rkbt, int order )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void validateIRKOrder( RKButcherTableau<Scalar> rkbt, int order )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void validateDIRKOrder( RKButcherTableau<Scalar> rkbt, int order )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
void validateSDIRKOrder( RKButcherTableau<Scalar> rkbt, int order )
{
  TEST_FOR_EXCEPT(true);
}

template<class Scalar>
class RKButcherTableauFactory 
{
  public:
    RKButcherTableau<Scalar> create(Teuchos::ParameterList& paramList) const = 0;
    RCP<const Teuchos::ParameterList> getValidParameters() const = 0;
};

template<class Scalar>
class DefaultRKButcherTableauFactory : virtual public RKButcherTableauFactory<Scalar>
{
  public:
    DefaultRKButcherTableauFactory();
    ~DefaultRKButcherTableauFactory();
    RKButcherTableau<Scalar> create(Teuchos::ParameterList& paramList) const; 
    RCP<const Teuchos::ParameterList> getValidParameters() const;
};


template<class Scalar>
RCP<const Teuchos::ParameterList> DefaultRKButcherTableauFactory<Scalar>::getValidParameters() const
{
  static RCP<Teuchos::ParameterList> validPL;

  if (Teuchos::is_null(validPL)) {
    RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    
    pl->set( SelectionType_name, SelectionType_default,
        "", rkButcherTableauSelectionTypeValidator );

    pl->set( SelectionTypeByName_name, SelectionTypeByName_default,
        "", rkButcherTableauMethodNameValidator );

    pl->set( SelectionTypeExplicitByOrder_name, SelectionTypeExplicitByOrder_default );
    pl->set( SelectionTypeImplicitByOrder_name, SelectionTypeImplicitByOrder_default );
    pl->set( SelectionTypeDIRKByOrder_name, SelectionTypeDIRKByOrder_default );
    pl->set( SelectionTypeSDIRKByOrder_name, SelectionTypeSDIRKByOrder_default );

    validPL = pl;
  }
  return (validPL);
}

template<class Scalar>
RKButcherTableau<Scalar> createBackwardEulerRKBT() 
{
  // Backward Euler Method:
  // c = [ 1 ]'
  // A = [ 1 ]
  // b = [ 1 ]'
  // This method is 1st order
  typedef ScalarTraits<Scalar> ST;
  Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
  A(0,0) = ST::one();
  Teuchos::SerialDenseVector<int,Scalar> b(1);
  b(0) = ST::one();
  Teuchos::SerialDenseVector<int,Scalar> c(1);
  c(0) = ST::one();
  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createForwardEulerRKBT() 
{
  // Forward Euler Method:
  // c = [ 0 ]'
  // A = [ 0 ]
  // b = [ 1 ]'
  // This method is 1st order
  typedef ScalarTraits<Scalar> ST;
  Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
  Teuchos::SerialDenseVector<int,Scalar> b(1);
  b(0) = ST::one();
  Teuchos::SerialDenseVector<int,Scalar> c(1);
  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit4StageRKBT()
{
  // "The" Runge-Kutta Method: 
  // "Solving Ordinary Differential Equations I:  Nonstiff Problems", 2nd Edition
  // E. Hairer, S.P. Norsett, G. Wanner
  // Table 1.2, pg 138
  // c = [  0  1/2 1/2  1  ]'
  // A = [  0              ] 
  //     [ 1/2  0          ]
  //     [  0  1/2  0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 1/3 1/3 1/6 ]'
  // This method is 4th order
  typedef ScalarTraits<Scalar> ST;
  Scalar one = ST::one();
  Scalar zero = ST::zero();
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar onesixth = ST::one()/(6*ST::one());
  Scalar onethird = ST::one()/(3*ST::one());

  int numStages = 4;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;
  A(0,3) = zero;

  A(1,0) = onehalf;
  A(1,1) = zero;
  A(1,2) = zero;
  A(1,3) = zero;

  A(2,0) = zero;
  A(2,1) = onehalf;
  A(2,2) = zero;
  A(2,3) = zero;

  A(3,0) = zero;
  A(3,1) = zero;
  A(3,2) = one;
  A(3,3) = zero;

  // Fill b:
  b(0) = onesixth;
  b(1) = onethird;
  b(2) = onethird;
  b(3) = onesixth;
  
  // fill b_c_
  c(0) = zero;
  c(1) = onehalf;
  c(2) = onehalf;
  c(3) = one;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit3_8RuleRKBT()
{
  // 3/8 Rule Runge-Kutta Method: 
  // "Solving Ordinary Differential Equations I:  Nonstiff Problems", 2nd Edition
  // E. Hairer, S.P. Norsett, G. Wanner
  // Table 1.2, pg 138
  // c = [  0  1/3 2/3  1  ]'
  // A = [  0              ]
  //     [ 1/3  0          ]
  //     [-1/3  1   0      ]
  //     [  1  -1   1   0  ]
  // b = [ 1/8 3/8 3/8 1/8 ]'
  // This method is 4th order
  typedef ScalarTraits<Scalar> ST;
  int numStages = 4;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  Scalar one = ST::one();
  Scalar zero = ST::zero();
  Scalar one_third    = Scalar(ST::one()/(3*ST::one()));
  Scalar two_third    = Scalar(2*ST::one()/(3*ST::one()));
  Scalar one_eighth   = Scalar(ST::one()/(8*ST::one()));
  Scalar three_eighth = Scalar(3*ST::one()/(8*ST::one()));

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;
  A(0,3) = zero;

  A(1,0) = one_third;
  A(1,1) = zero;
  A(1,2) = zero;
  A(1,3) = zero;

  A(2,0) = Scalar(-one_third);
  A(2,1) = one;
  A(2,2) = zero;
  A(2,3) = zero;

  A(3,0) = one;
  A(3,1) = Scalar(-one);
  A(3,2) = one;
  A(3,3) = zero;

  // Fill b:
  b(0) = one_eighth;
  b(1) = three_eighth;
  b(2) = three_eighth;
  b(3) = one_eighth;
  
  // Fill c:
  c(0) = zero;
  c(1) = one_third;
  c(2) = two_third;
  c(3) = one;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit4Stage3rdOrderRungeRKBT()
{
  // 4 Stage 3rd Order Explicit RK method, Runge
  // "Solving Ordinary Differential Equations I:  Nonstiff Problems", 2nd Edition
  // E. Hairer, S.P. Norsett, G. Wanner
  // Table 1.1, pg 135
  // c = [  0  1/2  1   1  ]'
  // A = [  0              ]
  //     [ 1/2  0          ]
  //     [  0   1   0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 2/3  0  1/6 ]'
  // This method is 3rd order
  typedef ScalarTraits<Scalar> ST;
  int numStages = 4;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  Scalar one = ST::one();
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar onesixth = ST::one()/(6*ST::one());
  Scalar twothirds = 2*ST::one()/(3*ST::one());
  Scalar zero = ST::zero();

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;
  A(0,3) = zero;

  A(1,0) = onehalf;
  A(1,1) = zero;
  A(1,2) = zero;
  A(1,3) = zero;

  A(2,0) = zero;
  A(2,1) = one;
  A(2,2) = zero;
  A(2,3) = zero;

  A(3,0) = zero;
  A(3,1) = zero;
  A(3,2) = one;
  A(3,3) = zero;

  // Fill b:
  b(0) = onesixth;
  b(1) = twothirds;
  b(2) = zero;
  b(3) = onesixth;
  
  // Fill c:
  c(0) = zero;
  c(1) = onehalf;
  c(2) = one;
  c(3) = one;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit3Stage3rdOrderRKBT()
{
  // 3 Stage 3rd Order Explicit RK Method
  // c = [  0  1/2  1  ]'
  // A = [  0          ] 
  //     [ 1/2  0      ]
  //     [ -1   2   0  ]
  // b = [ 1/6 4/6 1/6 ]'
  // This method is 3rd order
  typedef ScalarTraits<Scalar> ST;
  Scalar one = ST::one();
  Scalar two = Scalar(2*ST::one());
  Scalar zero = ST::zero();
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar onesixth = ST::one()/(6*ST::one());
  Scalar foursixth = 4*ST::one()/(6*ST::one());

  int numStages = 3;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;

  A(1,0) = onehalf;
  A(1,1) = zero;
  A(1,2) = zero;

  A(2,0) = -one;
  A(2,1) = two;
  A(2,2) = zero;

  // Fill b:
  b(0) = onesixth;
  b(1) = foursixth;
  b(2) = onesixth;
  
  // fill b_c_
  c(0) = zero;
  c(1) = onehalf;
  c(2) = one;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit3Stage3rdOrderHeunRKBT()
{
  // 3 Stage 3rd Order Explicit RK Method, Heun
  // "Solving Ordinary Differential Equations I:  Nonstiff Problems", 2nd Edition
  // E. Hairer, S.P. Norsett, G. Wanner
  // Table 1.1, pg 135
  // c = [  0  1/3 2/3 ]'
  // A = [  0          ] 
  //     [ 1/3  0      ]
  //     [  0  2/3  0  ]
  // b = [ 1/4  0  3/4 ]'
  // This method is 3rd order
  typedef ScalarTraits<Scalar> ST;
  Scalar one = ST::one();
  Scalar zero = ST::zero();
  Scalar onethird = one/(3*one);
  Scalar twothirds = 2*one/(3*one);
  Scalar onefourth = one/(4*one);
  Scalar threefourths = 3*one/(4*one);

  int numStages = 3;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;

  A(1,0) = onethird;
  A(1,1) = zero;
  A(1,2) = zero;

  A(2,0) = zero;
  A(2,1) = twothirds;
  A(2,2) = zero;

  // Fill b:
  b(0) = onefourth;
  b(1) = zero;
  b(2) = threefourths;
  
  // fill b_c_
  c(0) = zero;
  c(1) = onethird;
  c(2) = twothirds;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit2Stage2ndOrderRungeRKBT()
{
  // 2 Stage 2nd Order Explicit RK Method, Runge
  // "Solving Ordinary Differential Equations I:  Nonstiff Problems", 2nd Edition
  // E. Hairer, S.P. Norsett, G. Wanner
  // Table 1.1, pg 135
  // c = [  0  1/2 ]'
  // A = [  0      ] 
  //     [ 1/2  0  ]
  // b = [  0   1  ]'
  // This method is 2nd order
  typedef ScalarTraits<Scalar> ST;
  Scalar one = ST::one();
  Scalar zero = ST::zero();
  Scalar onehalf = ST::one()/(2*ST::one());

  int numStages = 2;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;

  A(1,0) = onehalf;
  A(1,1) = zero;

  // Fill b:
  b(0) = zero;
  b(1) = one;
  
  // fill b_c_
  c(0) = zero;
  c(1) = onehalf;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createExplicit4Stage3rdOrderRKBT()
{
  // 4 Stage 3rd Order Explicit RK Method
  // c = [  0  1/2 1/2  1  ]'
  // A = [  0              ]
  //     [ 1/2  0          ]
  //     [  0  1/2  0      ]
  //     [  0   0   1   0  ]
  // b = [ 1/6 2/6 2/6 1/6 ]'
  // This method is 3rd order
  typedef ScalarTraits<Scalar> ST;
  int numStages = 4;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);

  Scalar one = ST::one();
  Scalar zero = ST::zero();
  Scalar one_half = Scalar(ST::one()/(2*ST::one()));
  Scalar two_sixth    = Scalar(2*ST::one()/(6*ST::one()));
  Scalar one_sixth   = Scalar(ST::one()/(6*ST::one()));

  // Fill A:
  A(0,0) = zero;
  A(0,1) = zero;
  A(0,2) = zero;
  A(0,3) = zero;

  A(1,0) = one_half;
  A(1,1) = zero;
  A(1,2) = zero;
  A(1,3) = zero;

  A(2,0) = zero;
  A(2,1) = one_half;
  A(2,2) = zero;
  A(2,3) = zero;

  A(3,0) = zero;
  A(3,1) = zero;
  A(3,2) = one;
  A(3,3) = zero;

  // Fill b:
  b(0) = one_sixth;
  b(1) = two_sixth;
  b(2) = two_sixth;
  b(3) = one_sixth;
  
  // Fill c:
  c(0) = zero;
  c(1) = one_half;
  c(2) = one_half;
  c(3) = one;

  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createImplicit1Stage2ndOrderGaussRKBT()
{
  // 1 Stage 2nd order Gauss Implicit RK Method
  // "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems", 2nd Revised Edition
  // E. Hairer and G. Wanner
  // Table 5.2, pg 72
  // c = [ 1/2 ]'
  // A = [ 1/2 ]
  // b = [  1  ]'
  typedef ScalarTraits<Scalar> ST;
  int numStages = 1;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);
  Scalar onehalf = ST::one()/(2*ST::one());
  Scalar one = ST::one();
  A(0,0) = onehalf;
  b(0) = one;
  c(0) = onehalf;
  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createImplicit2Stage4thOrderGaussRKBT()
{
  // "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems", 2nd Revised Edition
  // E. Hairer and G. Wanner
  // Table 5.2, pg 72
  // c = [ 1/2-sqrt(3)/6  1/2+sqrt(3)/6 ]'
  // A = [ 1/4            1/4-sqrt(3)/6 ]
  //     [ 1/4+sqrt(3)/6  1/4           ]
  // b = [ 1/2 1/2 ]'
  typedef ScalarTraits<Scalar> ST;
  int numStages = 2;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);
  Scalar one = ST::one();
  Scalar onehalf = Scalar(one/(2*one));
  Scalar three = Scalar(3*one);
  Scalar six = Scalar(6*one);
  Scalar onefourth = Scalar(one/(4*one));
  Scalar alpha = sqrt(three)/six;

  A(0,0) = onefourth;
  A(0,1) = onefourth-alpha;
  A(1,0) = onefourth+alpha;
  A(1,1) = onefourth;
  b(0) = onehalf;
  b(1) = onehalf;
  c(0) = onehalf-alpha;
  c(1) = onehalf+alpha;
  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> createImplicit3Stage6thOrderGaussRKBT()
{
  // "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems", 2nd Revised Edition
  // E. Hairer and G. Wanner
  // Table 5.2, pg 72
  // c = [ 1/2-sqrt(15)/10   1/2              1/2+sqrt(15)/10  ]'
  // A = [ 5/36              2/9-sqrt(15)/15  5/36-sqrt(15)/30 ]
  //     [ 5/36+sqrt(15)/24  2/9              5/36-sqrt(15)/24 ]
  //     [ 5/36+sqrt(15)/30  2/9+sqrt(15)/15  5/36             ]
  // b = [ 5/18              4/9              5/18             ]'
  typedef ScalarTraits<Scalar> ST;
  int numStages = 3;
  Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
  Teuchos::SerialDenseVector<int,Scalar> b(numStages);
  Teuchos::SerialDenseVector<int,Scalar> c(numStages);
  Scalar one = ST::one();
  Scalar ten = Scalar(10*one);
  Scalar fifteen = Scalar(15*one);
  Scalar twentyfour = Scalar(24*one);
  Scalar thirty = Scalar(30*one);
  Scalar sqrt15over10 = Scalar(sqrt(fifteen)/ten);
  Scalar sqrt15over15 = Scalar(sqrt(fifteen)/fifteen);
  Scalar sqrt15over24 = Scalar(sqrt(fifteen)/twentyfour);
  Scalar sqrt15over30 = Scalar(sqrt(fifteen)/thirty);

  A(0,0) = Scalar(5*one/(36*one));
  A(0,1) = Scalar(2*one/(9*one))-sqrt15over15;
  A(0,2) = Scalar(5*one/(36*one))-sqrt15over30;
  A(1,0) = Scalar(5*one/(36*one))+sqrt15over24;
  A(1,1) = Scalar(2*one/(9*one));
  A(1,2) = Scalar(5*one/(36*one))-sqrt15over24;
  A(2,0) = Scalar(5*one/(36*one))+sqrt15over30;
  A(2,1) = Scalar(2*one/(9*one))+sqrt15over15;
  A(2,2) = Scalar(5*one/(36*one));
  b(0) = Scalar(5*one/(18*one));
  b(1) = Scalar(4*one/(9*one));
  b(2) = Scalar(5*one/(18*one));
  c(0) = Scalar(one/(2*one))-sqrt15over10;
  c(1) = Scalar(one/(2*one));
  c(2) = Scalar(one/(2*one))+sqrt15over10;
  return RKButcherTableau<Scalar>(A,b,c);
}

template<class Scalar>
RKButcherTableau<Scalar> DefaultRKButcherTableauFactory<Scalar>::create(Teuchos::ParameterList& paramList) const
{
  paramList.validateParameters(*this->getValidParameters());

  RKButcherTableau<Scalar> rkbt_out;

  E_RKButcherTableauSelectionTypes selectionTypeEnum = rkButcherTableauSelectionTypeValidator->getIntegralValue(
      paramList, SelectionType_name, SelectionType_default
      );
  if (selectionTypeEnum == RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_BY_NAME) {
    E_RKButcherTableauSelectionMethodNames selectionMethodEnum = rkButcherTableauMethodNameValidator->getIntegralValue(
        paramList, SelectionTypeByName_name, SelectionTypeByName_default
        );
    if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_FORWARD_EULER) {
      rkbt_out = createForwardEulerRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_BACKWARD_EULER) {
      rkbt_out = createBackwardEulerRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE) {
      rkbt_out = createExplicit4StageRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_8_RULE) {
      rkbt_out = createExplicit3_8RuleRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_2_STAGE_2_ORDER_RUNGE) {
      rkbt_out = createExplicit2Stage2ndOrderRungeRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER_HEUN) {
      rkbt_out = createExplicit3Stage3rdOrderHeunRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_3_STAGE_3_ORDER) {
      rkbt_out = createExplicit3Stage3rdOrderRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER_RUNGE) {
      rkbt_out = createExplicit4Stage3rdOrderRungeRKBT<Scalar>();
    } else if (selectionMethodEnum == RYTHMOS_RKBUTCHERTABLEAU_METHOD_NAME_EXPLICIT_4_STAGE_3_ORDER) {
      rkbt_out = createExplicit4Stage3rdOrderRKBT<Scalar>();
    } else {
      // Should never get here.
      TEST_FOR_EXCEPT(true);
    }
  } else if (selectionTypeEnum == RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_EXPLICIT_BY_ORDER) {
    int order = paramList.get<int>(SelectionTypeExplicitByOrder_name);
    if (order == 1) {
      rkbt_out = createForwardEulerRKBT<Scalar>();
    } else if (order == 2) {
      rkbt_out = createExplicit2Stage2ndOrderRungeRKBT<Scalar>();
    } else if (order == 3) {
      rkbt_out = createExplicit3Stage3rdOrderHeunRKBT<Scalar>();
    } else if (order == 4) {
      rkbt_out = createExplicit4StageRKBT<Scalar>();
    } else {
      TEST_FOR_EXCEPTION( order > 4, std::logic_error,
          "Error!  Explicit method by order is only implemented for orders 1 - 4\n"
          );
    }
  } else if (selectionTypeEnum == RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_IMPLICIT_BY_ORDER) {
    int order = paramList.get<int>(SelectionTypeImplicitByOrder_name);
    if (order == 1) {
      rkbt_out = createBackwardEulerRKBT<Scalar>();
    } else if (order == 2) {
      rkbt_out = createImplicit1Stage2ndOrderGaussRKBT<Scalar>();
    } else if (order == 4) {
      rkbt_out = createImplicit2Stage4thOrderGaussRKBT<Scalar>();
    } else if (order == 6) {
      rkbt_out = createImplicit3Stage6thOrderGaussRKBT<Scalar>();
    } else {
      TEST_FOR_EXCEPTION( order>1, std::logic_error,
          "Error!  Implicit method by order is only implemented for orders = 1,2,4,6 \n"
          );
    }
  } else if (selectionTypeEnum == RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_DIRK_BY_ORDER) {
    int order = paramList.get<int>(SelectionTypeDIRKByOrder_name);
    if (order == 1) {
      rkbt_out = createBackwardEulerRKBT<Scalar>();
    } else {
      TEST_FOR_EXCEPTION( order>1, std::logic_error,
          "Error!  DIRK method by order is not implemented for order > 1\n"
          );
    }
  } else if (selectionTypeEnum == RYTHMOS_RKBUTCHERTABLEAU_SELECTION_TYPE_SDIRK_BY_ORDER) {
    int order = paramList.get<int>(SelectionTypeSDIRKByOrder_name);
    if (order == 1) {
      rkbt_out = createBackwardEulerRKBT<Scalar>();
    } else {
      TEST_FOR_EXCEPTION( order>1, std::logic_error,
          "Error!  SDIRK method by order is not implemented for order > 1\n"
          );
    }
  } else {
    // Should never get here.
    TEST_FOR_EXCEPT(true);
  }
  
  return rkbt_out;
}

} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_HPP
