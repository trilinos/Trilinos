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

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

namespace Rythmos {

using Teuchos::SerialDenseMatrix;
using Teuchos::SerialDenseVector;
using Thyra::VectorBase;
using Thyra::ProductVectorBase;

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validInitialize ) {
  int numStages = 1;
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  RKButcherTableau<double> rkButcherTableau(A,b,c,1);

  TEST_EQUALITY( rkButcherTableau.numStages(), numStages );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidInitialize ) {

  RKButcherTableau<double> rkButcherTableau;

  int numStages = 1;
  {
    SerialDenseMatrix<int,double> A(numStages+1, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c,0), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages+1);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages+1);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages+1);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidAssembleIRKState ) {
  int N = 10;
  int numStages = 1;

  {
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages,numStages);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, outArg(*x)), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages+1,numStages);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, outArg(*x)), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages,numStages+1);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, outArg(*x)), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages,numStages);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages+1,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, outArg(*x)), 
        std::logic_error 
        );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleIRKState ) {

  {
    int N = 1;
    int numStages = 1;
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 5.0;
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,1.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,2.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,3.0);
    assembleIRKState(stageIndex, A, dt, *x_base, *x_stage_bar, outArg(*x));
    // What should x be?
    // x = x_base == 1.0      
    // x += dt*A(0,0)*x_stage_bar(0), so x = 1 + 0.1*5.0*2.0 = 2.0
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 2.0 );
  }
  {
    int N = 1;
    int numStages = 2;
    int stageIndex = 0;
    SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 1.0;
    A(0,1) = 2.0;
    A(1,0) = 3.0;
    A(1,1) = 4.0;
    double dt = 10.0;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,5.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,6.0);
    V_S(&*(x_stage_bar->getNonconstVectorBlock(1)),7.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,8.0);
    assembleIRKState(stageIndex, A, dt, *x_base, *x_stage_bar, outArg(*x));
    // What should x be?
    // x = x_base == 5.0               so x = 5
    // x += dt*A(0,0)*x_stage_bar(0)   so x = 5 + 10.0*1.0*6.0 = 65
    // x += dt*A(0,1)*x_stage_bar(1)   so x = 65 + 10.0*2.0*7.0 = 205
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 205 );

    assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, outArg(*x));
    // What should x be?
    // x = x_base == 5.0               so x = 5
    // x += dt*A(1,0)*x_stage_bar(0)   so x = 5 + 10.0*3.0*6.0 = 185
    // x += dt*A(1,1)*x_stage_bar(1)   so x = 185 + 10.0*4.0*7.0 = 465
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 465 );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidAssembleIRKSolution ) {
  int N = 10;
  int numStages = 1;

  {
    SerialDenseVector<int,double> b(numStages+1);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKSolution(b, dt, *x_base, *x_stage_bar, outArg(*x)),
        std::logic_error 
        );
  }
  {
    SerialDenseVector<int,double> b(numStages);
    double dt = 0.1;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,0.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages+1,N,0.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKSolution(b, dt, *x_base, *x_stage_bar, outArg(*x)),
        std::logic_error 
        );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleIRKSolution ) {

  {
    int N = 1;
    int numStages = 1;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 2.0;
    double dt = 10.0;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,3.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,4.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,5.0);
    assembleIRKSolution(b, dt, *x_base, *x_stage_bar, outArg(*x));
    // What should x be?
    // x = x_base == 3.0             so x = 3
    // x += dt*b(0)*x_stage_bar(0)   so x = 3 + 10.0*2.0*4.0 = 83
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 83 );
  }
  {
    int N = 1;
    int numStages = 2;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 2.0;
    b(1) = 3.0;
    double dt = 10.0;
    RCP<VectorBase<double> > x_base = createDefaultVector(N,4.0);
    RCP<ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,5.0);
    V_S(&*(x_stage_bar->getNonconstVectorBlock(1)),6.0);
    RCP<VectorBase<double> > x = createDefaultVector(N,7.0);
    assembleIRKSolution(b, dt, *x_base, *x_stage_bar, outArg(*x));
    // What should x be?
    // x = x_base == 4.0             so x = 4
    // x += dt*b(0)*x_stage_bar(0)   so x = 4 + 10.0*2.0*5.0 = 104
    // x += dt*b(1)*x_stage_bar(1)   so x = 104 + 10.0*3.0*6.0 = 284
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 284 );
  }
}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleERKState ) { 
//  // TODO:  Fill in tests
//}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleERKSolution ) { 
//  // TODO:  Fill in tests
//}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assertNonEmptyRKButcherTableau ) {
  {
    // Check that an empty tableau is thrown
    int numStages = 0;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // Check that a zero tableau is thrown
    int numStages = 1;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // Check that a non-zero A, non-zero c, but zero b throws.
    int numStages = 2;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(0,0) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    c(0) = 0.5;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // Check that a valid tableau, valid b, valid c is NOT thrown
    int numStages = 1;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(0,0) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 1.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    assertNonEmptyRKButcherTableau(rkButcherTableau);
  }

}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validateERKButcherTableau ) {
  {
    // Check that the Backward Euler tableau is invalid
    int numStages = 1;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(0,0) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 1.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // Check that the Forward Euler tableau is valid
    int numStages = 1;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(0,0) = 0.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 0.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    validateERKButcherTableau(rkButcherTableau);
  }
  {
    // Check that c(0) != 0 throws
    int numStages = 3;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(1,0) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 1.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // verify throw for entries on diagonal
    int numStages = 3;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(2,2) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 0.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(rkButcherTableau), std::logic_error );
  }
  {
    // verify throw for entries in upper triangle
    int numStages = 3;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(1,2) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    b(0) = 1.0;
    SerialDenseVector<int,double> c(numStages);
    c(0) = 0.0;
    RKButcherTableau<double> rkButcherTableau = RKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(rkButcherTableau), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createForwardEulerRKBT ) {
  RKButcherTableau<double> rkbt = createForwardEulerRKBT<double>();
  TEST_EQUALITY_CONST( rkbt.numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt.A()(0,0), 0.0 );
  TEST_EQUALITY_CONST( rkbt.b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt.c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createBackwardEulerRKBT ) {
  RKButcherTableau<double> rkbt = createBackwardEulerRKBT<double>();
  TEST_EQUALITY_CONST( rkbt.numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt.A()(0,0), 1.0 );
  TEST_EQUALITY_CONST( rkbt.b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt.c()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit4StageRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit4StageRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate4thOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 4 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0), 0.5 );
  TEST_EQUALITY_CONST( A(2,0), 0.0 );
  TEST_EQUALITY_CONST( A(2,1), 0.5 );
  TEST_EQUALITY_CONST( A(3,0), 0.0 );
  TEST_EQUALITY_CONST( A(3,1), 0.0 );
  TEST_EQUALITY_CONST( A(3,2), 1.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/6.0 );
  TEST_EQUALITY_CONST( b(1), 1.0/3.0 );
  TEST_EQUALITY_CONST( b(2), 1.0/3.0 );
  TEST_EQUALITY_CONST( b(3), 1.0/6.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( c(2), 0.5 );
  TEST_EQUALITY_CONST( c(3), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3_8RuleRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit3_8RuleRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate4thOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 4 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  1.0/3.0 );
  TEST_EQUALITY_CONST( A(2,0), -1.0/3.0 );
  TEST_EQUALITY_CONST( A(2,1),  1.0 );
  TEST_EQUALITY_CONST( A(3,0),  1.0 );
  TEST_EQUALITY_CONST( A(3,1), -1.0 );
  TEST_EQUALITY_CONST( A(3,2),  1.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/8.0 );
  TEST_EQUALITY_CONST( b(1), 3.0/8.0 );
  TEST_EQUALITY_CONST( b(2), 3.0/8.0 );
  TEST_EQUALITY_CONST( b(3), 1.0/8.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 1.0/3.0 );
  TEST_EQUALITY_CONST( c(2), 2.0/3.0 );
  TEST_EQUALITY_CONST( c(3), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit2Stage2ndOrderRungeRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit2Stage2ndOrderRungeRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate2ndOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  0.5 );
  TEST_EQUALITY_CONST( b(0), 0.0 );
  TEST_EQUALITY_CONST( b(1), 1.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( rkbt.order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit4Stage3rdOrderRungeRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit4Stage3rdOrderRungeRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 4 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  0.5 );
  TEST_EQUALITY_CONST( A(2,0),  0.0 );
  TEST_EQUALITY_CONST( A(2,1),  1.0 );
  TEST_EQUALITY_CONST( A(3,0),  0.0 );
  TEST_EQUALITY_CONST( A(3,1),  0.0 );
  TEST_EQUALITY_CONST( A(3,2),  1.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/6.0 );
  TEST_EQUALITY_CONST( b(1), 2.0/3.0 );
  TEST_EQUALITY_CONST( b(2), 0.0     );
  TEST_EQUALITY_CONST( b(3), 1.0/6.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( c(2), 1.0 );
  TEST_EQUALITY_CONST( c(3), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3Stage3rdOrderHeunRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit3Stage3rdOrderHeunRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  1.0/3.0 );
  TEST_EQUALITY_CONST( A(2,0),  0.0     );
  TEST_EQUALITY_CONST( A(2,1),  2.0/3.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/4.0 );
  TEST_EQUALITY_CONST( b(1), 0.0     );
  TEST_EQUALITY_CONST( b(2), 3.0/4.0 );
  TEST_EQUALITY_CONST( c(0), 0.0     );
  TEST_EQUALITY_CONST( c(1), 1.0/3.0 );
  TEST_EQUALITY_CONST( c(2), 2.0/3.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3Stage3rdOrderRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit3Stage3rdOrderRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  0.5 );
  TEST_EQUALITY_CONST( A(2,0), -1.0 );
  TEST_EQUALITY_CONST( A(2,1),  2.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/6.0 );
  TEST_EQUALITY_CONST( b(1), 4.0/6.0 );
  TEST_EQUALITY_CONST( b(2), 1.0/6.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( c(2), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit4Stage3rdOrderRKBT ) {
  RKButcherTableau<double> rkbt = createExplicit4Stage3rdOrderRKBT<double>();
  validateERKButcherTableau(rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 4 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(1,0),  0.5 );
  TEST_EQUALITY_CONST( A(2,0),  0.0 );
  TEST_EQUALITY_CONST( A(2,1),  0.5 );
  TEST_EQUALITY_CONST( A(3,0),  0.0 );
  TEST_EQUALITY_CONST( A(3,1),  0.0 );
  TEST_EQUALITY_CONST( A(3,2),  1.0 );
  TEST_EQUALITY_CONST( b(0), 1.0/6.0 );
  TEST_EQUALITY_CONST( b(1), 2.0/6.0 );
  TEST_EQUALITY_CONST( b(2), 2.0/6.0 );
  TEST_EQUALITY_CONST( b(3), 1.0/6.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( c(2), 0.5 );
  TEST_EQUALITY_CONST( c(3), 1.0 );
  TEST_EQUALITY_CONST( rkbt.order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit1Stage2ndOrderGaussRKBT ) {
  RKButcherTableau<double> rkbt = createImplicit1Stage2ndOrderGaussRKBT<double>();
  validateIRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 1 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_EQUALITY_CONST( A(0,0),  0.5 );
  TEST_EQUALITY_CONST( b(0), 1.0 );
  TEST_EQUALITY_CONST( c(0), 0.5 );
  TEST_EQUALITY_CONST( rkbt.order(), 2 );
}
TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage4thOrderGaussRKBT ) {
  double tol = 1.0e-10;
  RKButcherTableau<double> rkbt = createImplicit2Stage4thOrderGaussRKBT<double>();
  validateIRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_FLOATING_EQUALITY( A(0,0),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.25-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  0.25+sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  0.25, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.5, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.5-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5+sqrt(3.0)/6.0, tol );
  TEST_EQUALITY_CONST( rkbt.order(), 4 );
}
TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage6thOrderGaussRKBT ) {
  double tol = 1.0e-10;
  RKButcherTableau<double> rkbt = createImplicit3Stage6thOrderGaussRKBT<double>();
  validateIRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_FLOATING_EQUALITY( A(0,0),  5.0/36.0, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  2.0/9.0-sqrt(15.0)/15.0, tol );
  TEST_FLOATING_EQUALITY( A(0,2),  5.0/36.0-sqrt(15.0)/30.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  5.0/36.0+sqrt(15.0)/24.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  2.0/9.0, tol );
  TEST_FLOATING_EQUALITY( A(1,2),  5.0/36.0-sqrt(15.0)/24.0, tol );
  TEST_FLOATING_EQUALITY( A(2,0),  5.0/36.0+sqrt(15.0)/30.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  2.0/9.0+sqrt(15.0)/15.0, tol );
  TEST_FLOATING_EQUALITY( A(2,2),  5.0/36.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 5.0/18.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 4.0/9.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 5.0/18.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.5-sqrt(15.0)/10.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(2), 0.5+sqrt(15.0)/10.0, tol );
  TEST_EQUALITY_CONST( rkbt.order(), 6 );
}

/*
TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createSDIRK5Stage5thOrderRKBT ) {
  double tol = 1.0e-10;
  RKButcherTableau<double> rkbt = createSDIRK5Stage5thOrderRKBT<double>();
  validateSDIRKButcherTableau(rkbt);
  TEST_EQUALITY_CONST( rkbt.numStages(), 5 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt.A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt.b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt.c();
  TEST_FLOATING_EQUALITY( A(0,0),  (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,2),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(1,0),  (-6.0+5.0*sqrt(6.0))/14.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( A(1,2),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(2,0),  (888.0+607.0*sqrt(6.0))/2850.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  (126.0-161.0*sqrt(6.0))/1425.0, tol );
  TEST_FLOATING_EQUALITY( A(2,2),  (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( A(2,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(2,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(3,0),  (3153.0-3082.0*sqrt(6.0))/14250.0, tol );
  TEST_FLOATING_EQUALITY( A(3,1),  (3213.0+1148.0*sqrt(6.0))/28500.0, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  (-267.0+88.0*sqrt(6.0))/500.0, tol );
  TEST_FLOATING_EQUALITY( A(3,3),  (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( A(3,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(4,0),  (-32583.0+14638.0*sqrt(6.0))/71250.0, tol );
  TEST_FLOATING_EQUALITY( A(4,1),  (-17199.0+364.0*sqrt(6.0))/142500.0, tol );
  TEST_FLOATING_EQUALITY( A(4,2),  (1329.0-544.0*sqrt(6.0))/2500.0, tol );
  TEST_FLOATING_EQUALITY( A(4,3),  (-96.0+131.0*sqrt(6.0))/625.0, tol );
  TEST_FLOATING_EQUALITY( A(4,4),  (6.0-sqrt(6.0))/10.0, tol );

  TEST_FLOATING_EQUALITY( b(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( b(3), (16.0-sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( b(4), (16.0+sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( c(0), (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( c(1), (6.0+9.0*sqrt(6.0))/35.0, tol );
  TEST_FLOATING_EQUALITY( c(2), 1.0, tol );
  TEST_FLOATING_EQUALITY( c(3), (4.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( c(4), (4.0+sqrt(6.0))/10.0, tol );
  TEST_EQUALITY_CONST( rkbt.order(), 5 );
}
*/

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validateDIRKButcherTableau ) {
  {
    // Entries above the diagonal should throw
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,2) = 1.0;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    RKButcherTableau<double> rkbt(A,b,c,1);
    TEST_THROW(validateDIRKButcherTableau(rkbt),std::logic_error);
  }
  {
    // Valid DIRK tableau
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 1.0;
    A(1,0) = 0.5;
    A(1,1) = 2.0;
    A(2,0) = 3.0;
    A(2,1) = 4.0;
    A(2,2) = 5.0;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(0) = 0.5;
    b(1) = 0.75;
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    c(0) = 0.2;
    c(1) = 0.5;
    c(2) = 0.75;
    RKButcherTableau<double> rkbt(A,b,c,1);
    validateDIRKButcherTableau(rkbt);
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validateSDIRKButcherTableau ) {
  {
    // Entries above the diagonal should throw
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,2) = 1.0;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    RKButcherTableau<double> rkbt(A,b,c,1);
    TEST_THROW(validateSDIRKButcherTableau(rkbt),std::logic_error);
  }
  {
    // The diagonal values must all be the same.
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 1.0;
    A(1,0) = 0.5;
    A(1,1) = 2.0;
    A(2,0) = 3.0;
    A(2,1) = 4.0;
    A(2,2) = 5.0;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(0) = 0.5;
    b(1) = 0.75;
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    c(0) = 0.2;
    c(1) = 0.5;
    c(2) = 0.75;
    RKButcherTableau<double> rkbt(A,b,c,1);
    TEST_THROW(validateSDIRKButcherTableau(rkbt),std::logic_error);
  }
  {
    // valid SDIRK tableau
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 0.9;
    A(1,0) = 0.5;
    A(1,1) = 0.9;
    A(2,0) = 3.0;
    A(2,1) = 4.0;
    A(2,2) = 0.9;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(0) = 0.5;
    b(1) = 0.75;
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    c(0) = 0.2;
    c(1) = 0.5;
    c(2) = 0.75;
    RKButcherTableau<double> rkbt(A,b,c,1);
    validateSDIRKButcherTableau(rkbt);
  }
}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableauFactory, validParameters ) {
//  TEST_EQUALITY_CONST( true, false );
//}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableauFactory, create ) {
//  TEST_EQUALITY_CONST( true, false );
//}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, operatorEqualEqual ) {
  {
    RKButcherTableau<double> rkbt_1 = createExplicit4StageRKBT<double>();
    {
      RKButcherTableau<double> rkbt_2 = createExplicit4StageRKBT<double>();
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, true ); // same tableau
    }
    {
      RKButcherTableau<double> rkbt_2 = createExplicit3_8RuleRKBT<double>();
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // different values in A,b,c
    }
    {
      RKButcherTableau<double> rkbt_2 = createBackwardEulerRKBT<double>();
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // different number of stages
    }
    {
      RKButcherTableau<double> rkbt_2 = createImplicit3Stage6thOrderGaussRKBT<double>();
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // different number of stages
    }
  }

  {
    // differing values in A then b then c:
    int numStages = 2;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    Teuchos::SerialDenseVector<int,double> b(numStages);
    Teuchos::SerialDenseVector<int,double> c(numStages);
    A(0,0) = 2.0/3.0;
    A(0,1) = 0.5;
    A(1,0) = 0.5;
    A(1,1) = 1.0;
    b(0) = 0.0;
    b(1) = 1.0;
    c(0) = 0.0;
    c(1) = 1.0;
    int order = 2;
    RKButcherTableau<double> rkbt_1(A,b,c,order);
    {
      Teuchos::SerialDenseMatrix<int,double> Aprime(numStages,numStages);
      Aprime(0,0) = 1.0/1.5;
      Aprime(0,1) = 0.5;
      Aprime(1,0) = 0.5;
      Aprime(1,1) = 1.0;
      RKButcherTableau<double> rkbt_2(Aprime,b,c,order);
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, true ); // Same values in A
    }
    {
      Teuchos::SerialDenseMatrix<int,double> Aprime(numStages,numStages);
      Aprime(0,0) = 2.0/3.0;
      Aprime(0,1) = 0.5;
      Aprime(1,0) = 0.5;
      Aprime(1,1) = 1.1;
      RKButcherTableau<double> rkbt_2(Aprime,b,c,order);
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // differing values in A
    }
    {
      Teuchos::SerialDenseVector<int,double> bprime(numStages);
      bprime(0) = 0.0;
      bprime(1) = 0.9;
      RKButcherTableau<double> rkbt_2(A,bprime,c,order);
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // differing values in b
    }
    {
      Teuchos::SerialDenseVector<int,double> cprime(numStages);
      cprime(0) = 0.0;
      cprime(1) = 0.9;
      RKButcherTableau<double> rkbt_2(A,b,cprime,order);
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // differing values in c
    }
    {
      int orderprime = 3;
      RKButcherTableau<double> rkbt_2(A,b,c,orderprime);
      TEST_EQUALITY_CONST( rkbt_1 == rkbt_2, false ); // differing values in order
    }
  }
}

} // namespace Rythmos

