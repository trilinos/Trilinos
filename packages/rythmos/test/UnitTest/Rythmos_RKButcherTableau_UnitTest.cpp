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
#include "Rythmos_RKButcherTableauHelpers.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

namespace Rythmos {

using Teuchos::SerialDenseMatrix;
using Teuchos::SerialDenseVector;
using Teuchos::outArg;
using Thyra::VectorBase;
using Thyra::ProductVectorBase;

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validInitialize ) {
  int numStages = 1;
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);

  TEST_EQUALITY( rkButcherTableau->numStages(), numStages );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidInitialize ) {

  RCP<RKButcherTableauBase<double> > rkButcherTableau;

  int numStages = 1;
  {
    SerialDenseMatrix<int,double> A(numStages+1, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = rKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = rKButcherTableau<double>(A,b,c,0), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages+1);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = rKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages+1);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = rKButcherTableau<double>(A,b,c,1), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages+1);
    TEST_THROW( rkButcherTableau = rKButcherTableau<double>(A,b,c,1), std::logic_error );
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
    V_S(x_stage_bar->getNonconstVectorBlock(1).ptr(),7.0);
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
    V_S(x_stage_bar->getNonconstVectorBlock(1).ptr(),6.0);
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(*rkButcherTableau), std::logic_error );
  }
  {
    // Check that a zero tableau is thrown
    int numStages = 1;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(*rkButcherTableau), std::logic_error );
  }
  {
    // Check that a non-zero A, non-zero c, but zero b throws.
    int numStages = 2;
    SerialDenseMatrix<int,double> A(numStages, numStages);
    A(0,0) = 1.0;
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    c(0) = 0.5;
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( assertNonEmptyRKButcherTableau(*rkButcherTableau), std::logic_error );
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    assertNonEmptyRKButcherTableau(*rkButcherTableau);
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(*rkButcherTableau), std::logic_error );
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    validateERKButcherTableau(*rkButcherTableau);
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(*rkButcherTableau), std::logic_error );
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(*rkButcherTableau), std::logic_error );
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
    RCP<RKButcherTableauBase<double> > rkButcherTableau = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW( validateERKButcherTableau(*rkButcherTableau), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createForwardEuler_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new ForwardEuler_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createBackwardEuler_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new BackwardEuler_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit4Stage4thOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit4Stage4thOrder_RKBT<double>());
  double tol = 1.0e-10;
  validateERKButcherTableau(*rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( A(1,0),  0.5, tol );
  TEST_FLOATING_EQUALITY( A(2,0),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  0.5, tol );
  TEST_FLOATING_EQUALITY( A(3,0),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(3,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  1.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( b(3), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(2), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3_8Rule_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit3_8Rule_RKBT<double>());
  double tol = 1.0e-10;
  validateERKButcherTableau(*rkbt);
  //validate4thOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( A(1,0),  1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( A(2,0), -1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  1.0, tol );
  TEST_FLOATING_EQUALITY( A(3,0),  1.0, tol );
  TEST_FLOATING_EQUALITY( A(3,1), -1.0, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  1.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 1.0/8.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 3.0/8.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 3.0/8.0, tol );
  TEST_FLOATING_EQUALITY( b(3), 1.0/8.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( c(2), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( c(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit2Stage2ndOrderRunge_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit2Stage2ndOrderRunge_RKBT<double>());
  validateERKButcherTableau(*rkbt);
  //validate2ndOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_EQUALITY_CONST( A(1,0),  0.5 );
  TEST_EQUALITY_CONST( b(0), 0.0 );
  TEST_EQUALITY_CONST( b(1), 1.0 );
  TEST_EQUALITY_CONST( c(0), 0.0 );
  TEST_EQUALITY_CONST( c(1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit4Stage3rdOrderRunge_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit4Stage3rdOrderRunge_RKBT<double>());
  double tol = 1.0e-10;
  validateERKButcherTableau(*rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( A(1,0),  0.5, tol );
  TEST_FLOATING_EQUALITY( A(2,0),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  1.0, tol );
  TEST_FLOATING_EQUALITY( A(3,0),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(3,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  1.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 0.0    , tol );
  TEST_FLOATING_EQUALITY( b(3), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(2), 1.0, tol );
  TEST_FLOATING_EQUALITY( c(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3Stage3rdOrderHeun_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit3Stage3rdOrderHeun_RKBT<double>());
  double tol = 1.0e-10;
  validateERKButcherTableau(*rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( A(1,0),  1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( A(2,0),  0.0    , tol );
  TEST_FLOATING_EQUALITY( A(2,1),  2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 1.0/4.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.0    , tol );
  TEST_FLOATING_EQUALITY( b(2), 3.0/4.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0    , tol );
  TEST_FLOATING_EQUALITY( c(1), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( c(2), 2.0/3.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createExplicit3Stage3rdOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Explicit3Stage3rdOrder_RKBT<double>());
  double tol = 1.0e-10;
  validateERKButcherTableau(*rkbt);
  //validate3rdOrderRKButcherTableau(rkbt);
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( A(1,0),  0.5, tol );
  TEST_FLOATING_EQUALITY( A(2,0), -1.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  2.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( b(1), 4.0/6.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(2), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createSDIRK2Stage3rdOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new SDIRK2Stage3rdOrder_RKBT<double>());
  validateSDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double gamma = (3.0+sqrt(3.0))/6.0; // could also be (3-sqrt(3))/6
  double tol=1.0e-10;
  TEST_FLOATING_EQUALITY( A(0,0),  gamma, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  1.0-2.0*gamma, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  gamma, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.5, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(0), gamma, tol );
  TEST_FLOATING_EQUALITY( c(1), 1.0-gamma, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, SDIRK2Stage3rdOrder_RKBT_pl ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new SDIRK2Stage3rdOrder_RKBT<double>());
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set<int>("gamma coefficient",-1);
  rkbt->setParameterList(pl);

  validateSDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double gamma = (3.0-sqrt(3.0))/6.0; 
  double tol=1.0e-10;
  TEST_FLOATING_EQUALITY( A(0,0),  gamma, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  1.0-2.0*gamma, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  gamma, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.5, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(0), gamma, tol );
  TEST_FLOATING_EQUALITY( c(1), 1.0-gamma, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createDIRK2Stage3rdOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new DIRK2Stage3rdOrder_RKBT<double>());
  validateDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double tol=1.0e-10;
  TEST_FLOATING_EQUALITY( A(0,0),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.25, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.75, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 2.0/3.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}


TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage4thOrderHammerHollingsworth_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage4thOrderHammerHollingsworth_RKBT<double>());
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double tol=1.0e-10;
  TEST_FLOATING_EQUALITY( A(0,0),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.25-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  0.25+sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  0.25, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.5, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.5-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5+sqrt(3.0)/6.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage6thOrderKuntzmannButcher_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage6thOrderKuntzmannButcher_RKBT<double>());
  double tol = 1.0e-10;
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
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
  TEST_EQUALITY_CONST( rkbt->order(), 6 );
}


TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit4Stage8thOrderKuntzmannButcher_RKBT ) {
  // This RKBT doesn't pass convergence testing, so its excluded from the factory.
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit4Stage8thOrderKuntzmannButcher_RKBT<double>());
  double tol = 1.0e-10;
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double w1 = 1.0/8.0-sqrt(30.0)/144.0;
  double w2 = 0.5*sqrt((15.0+2.0*sqrt(30.0))/35.0);
  double w3 = w2*(1.0/6.0+sqrt(30.0)/24.0);
  double w4 = w2*(1.0/21.0+5.0*sqrt(30.0)/168.0);
  double w5 = w2-2*w3;
  double w1p = 1.0/8.0+sqrt(30.0)/144.0;
  double w2p = 0.5*sqrt((15.0-2.0*sqrt(30.0))/35.0);
  double w3p = w2p*(1.0/6.0-sqrt(30.0)/24.0);
  double w4p = w2p*(1.0/21.0-5.0*sqrt(30.0)/168.0);
  double w5p = w2p-2*w3p;
  TEST_FLOATING_EQUALITY( A(0,0),  w1, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  w1p-w3+w4p, tol );
  TEST_FLOATING_EQUALITY( A(0,2),  w1p-w3-w4p, tol );
  TEST_FLOATING_EQUALITY( A(0,3),  w1-w5, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  w1-w3p+w4, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  w1p, tol );
  TEST_FLOATING_EQUALITY( A(1,2),  w1p-w5p, tol );
  TEST_FLOATING_EQUALITY( A(1,3),  w1-w3p-w4, tol );
  TEST_FLOATING_EQUALITY( A(2,0),  w1+w3p+w4, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  w1p+w5p, tol );
  TEST_FLOATING_EQUALITY( A(2,2),  w1p, tol );
  TEST_FLOATING_EQUALITY( A(2,3),  w1+w3p-w4, tol );
  TEST_FLOATING_EQUALITY( A(3,0),  w1+w5, tol );
  TEST_FLOATING_EQUALITY( A(3,1),  w1p+w3+w4p, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  w1p+w3-w4p, tol );
  TEST_FLOATING_EQUALITY( A(3,3),  w1, tol );
  TEST_FLOATING_EQUALITY( b(0), 2.0*w1, tol );
  TEST_FLOATING_EQUALITY( b(1), 2.0*w1p, tol );
  TEST_FLOATING_EQUALITY( b(2), 2.0*w1p, tol );
  TEST_FLOATING_EQUALITY( b(3), 2.0*w1, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.5-w2, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5-w2p, tol );
  TEST_FLOATING_EQUALITY( c(2), 0.5+w2p, tol );
  TEST_FLOATING_EQUALITY( c(3), 0.5+w2, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 8 );
}


TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit1Stage2ndOrderGauss_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit1Stage2ndOrderGauss_RKBT<double>());
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 1 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_EQUALITY_CONST( A(0,0),  0.5 );
  TEST_EQUALITY_CONST( b(0), 1.0 );
  TEST_EQUALITY_CONST( c(0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage4thOrderGauss_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage4thOrderGauss_RKBT<double>());
  double tol = 1.0e-10;
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_FLOATING_EQUALITY( A(0,0),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.25-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,0),  0.25+sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  0.25, tol );
  TEST_FLOATING_EQUALITY( b(0), 0.5, tol );
  TEST_FLOATING_EQUALITY( b(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(0), 0.5-sqrt(3.0)/6.0, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5+sqrt(3.0)/6.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage6thOrderGauss_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage6thOrderGauss_RKBT<double>());
  double tol = 1.0e-10;
  validateIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
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
  TEST_EQUALITY_CONST( rkbt->order(), 6 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit1Stage1stOrderRadauA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit1Stage1stOrderRadauA_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage3rdOrderRadauA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage3rdOrderRadauA_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), -0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 0.75, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), 2.0/3.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage5thOrderRadauA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage5thOrderRadauA_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), (-1.0-sqrt(6.0))/18.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), (-1.0+sqrt(6.0))/18.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), (88.0+7.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), (88.0-43.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), (88.0+43.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), (88.0-7.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), (16.0+sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), (16.0-sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), (6.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), (6.0+sqrt(6.0))/10.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 5 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit1Stage1stOrderRadauB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit1Stage1stOrderRadauB_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 1 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 1.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 1 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage3rdOrderRadauB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage3rdOrderRadauB_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), -1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 0.75, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 0.75, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 3 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage5thOrderRadauB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage5thOrderRadauB_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), (88.0-7.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), (296.0-169.0*sqrt(6.0))/1800.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), (-2.0+3.0*sqrt(6.0))/225.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), (296.0+169.0*sqrt(6.0))/1800.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), (88.0+7.0*sqrt(6.0))/360.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), (-2.0-3.0*sqrt(6.0))/225.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), (16.0-sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), (16.0+sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), (16.0-sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), (16.0+sqrt(6.0))/36.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 1.0/9.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), (4.0-sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), (4.0+sqrt(6.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 5 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage2ndOrderLobattoA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage2ndOrderLobattoA_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->A()(0,1), 0.0 );
  TEST_EQUALITY_CONST( rkbt->A()(1,0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(1,1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->b()(1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->c()(1), 1.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage4thOrderLobattoA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage4thOrderLobattoA_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 5.0/24.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), -1.0/24.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit4Stage6thOrderLobattoA_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit4Stage6thOrderLobattoA_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,3), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), (11.0+sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), (25.0-sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), (25.0-13.0*sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,3), (-1.0+sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), (11.0-sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), (25.0+13.0*sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), (25.0+sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,3), (-1.0-sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,2), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,3), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(3), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), (5.0-sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), (5.0+sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 6 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage2ndOrderLobattoB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage2ndOrderLobattoB_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(0,1), 0.0 );
  TEST_EQUALITY_CONST( rkbt->A()(1,0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(1,1), 0.0 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->b()(1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->c()(1), 1.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage4thOrderLobattoB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage4thOrderLobattoB_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), -1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), 5.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit4Stage6thOrderLobattoB_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit4Stage6thOrderLobattoB_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), (-1.0-sqrt(5.0))/24.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), (-1.0+sqrt(5.0))/24.0 , tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,3), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), (25.0+sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), (25.0-13.0*sqrt(5.0))/120.0 , tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,3), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), (25.0+13.0*sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), (25.0-sqrt(5.0))/120.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,3), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,1), (11.0-sqrt(5.0))/24.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,2), (11.0+sqrt(5.0))/24.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,3), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(3), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), (5.0-sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), (5.0+sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 6 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit2Stage2ndOrderLobattoC_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit2Stage2ndOrderLobattoC_RKBT<double>());
  TEST_EQUALITY_CONST( rkbt->numStages(), 2 );
  TEST_EQUALITY_CONST( rkbt->A()(0,0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(0,1), -0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(1,0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->A()(1,1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->b()(0), 0.5 );
  TEST_EQUALITY_CONST( rkbt->b()(1), 0.5 );
  TEST_EQUALITY_CONST( rkbt->c()(0), 0.0 );
  TEST_EQUALITY_CONST( rkbt->c()(1), 1.0 );
  TEST_EQUALITY_CONST( rkbt->order(), 2 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit3Stage4thOrderLobattoC_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit3Stage4thOrderLobattoC_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), -1.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), -1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 2.0/3.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 1.0/6.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createImplicit4Stage6thOrderLobattoC_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new Implicit4Stage6thOrderLobattoC_RKBT<double>());
  double tol = 1.0e-10;
  TEST_EQUALITY_CONST( rkbt->numStages(), 4 );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,1), -sqrt(5.0)/12, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,2), +sqrt(5.0)/12, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(0,3), -1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,1), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,2), (10.0-7.0*sqrt(5.0))/60.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(1,3), sqrt(5.0)/60.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,1), (10.0+7.0*sqrt(5.0))/60.0 , tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,2), 0.25, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(2,3), -sqrt(5.0)/60.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,2), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->A()(3,3), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(0), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(1), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(2), 5.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->b()(3), 1.0/12.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(0), 0.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(1), (5.0-sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(2), (5.0+sqrt(5.0))/10.0, tol );
  TEST_FLOATING_EQUALITY( rkbt->c()(3), 1.0, tol );
  TEST_EQUALITY_CONST( rkbt->order(), 6 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createSDIRK5Stage5thOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new SDIRK5Stage5thOrder_RKBT<double>());
  //TEST_THROW(validateSDIRKButcherTableau<double>(rkbt), std::logic_error );
  double tol = 1.0e-10;
  validateSDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 5 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
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

  TEST_EQUALITY_CONST( rkbt->order(), 5 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createSDIRK5Stage4thOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new SDIRK5Stage4thOrder_RKBT<double>());
  double tol = 1.0e-10;
  validateSDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 5 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  TEST_FLOATING_EQUALITY( A(0,0),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,2),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(1,0),  0.5, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(1,2),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(1,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(2,0),  17.0/50.0, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  -1.0/25.0, tol );
  TEST_FLOATING_EQUALITY( A(2,2),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(2,3),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(2,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(3,0),  371.0/1360.0, tol );
  TEST_FLOATING_EQUALITY( A(3,1),  -137.0/2720.0, tol );
  TEST_FLOATING_EQUALITY( A(3,2),  15.0/544.0, tol );
  TEST_FLOATING_EQUALITY( A(3,3),  0.25, tol );
  TEST_FLOATING_EQUALITY( A(3,4),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(4,0),  25.0/24.0, tol );
  TEST_FLOATING_EQUALITY( A(4,1),  -49.0/48.0, tol );
  TEST_FLOATING_EQUALITY( A(4,2),  125.0/16.0, tol );
  TEST_FLOATING_EQUALITY( A(4,3),  -85.0/12.0, tol );
  TEST_FLOATING_EQUALITY( A(4,4),  0.25, tol );

  TEST_FLOATING_EQUALITY( b(0), 25.0/24.0, tol );
  TEST_FLOATING_EQUALITY( b(1), -49.0/48.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 125.0/16.0, tol );
  TEST_FLOATING_EQUALITY( b(3), -85.0/12.0, tol );
  TEST_FLOATING_EQUALITY( b(4), 0.25, tol );

  /* 
  //Alternate version
  TEST_FLOATING_EQUALITY( b(0), 59.0/48.0, tol );
  TEST_FLOATING_EQUALITY( b(1), -17.0/96.0, tol );
  TEST_FLOATING_EQUALITY( b(2), 225.0/32.0, tol );
  TEST_FLOATING_EQUALITY( b(3), -85.0/12.0, tol );
  TEST_FLOATING_EQUALITY( b(4), 0.0, tol );
  */

  TEST_FLOATING_EQUALITY( c(0), 0.25, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.75, tol );
  TEST_FLOATING_EQUALITY( c(2), 11.0/20.0, tol );
  TEST_FLOATING_EQUALITY( c(3), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(4), 1.0, tol );

  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createSDIRK3Stage4thOrder_RKBT ) {
  RCP<RKButcherTableauBase<double> > rkbt = rcp(new SDIRK3Stage4thOrder_RKBT<double>());
  //TEST_THROW(validateSDIRKButcherTableau<double>(rkbt), std::logic_error );
  double tol = 1.0e-10;
  validateSDIRKButcherTableau(*rkbt);
  TEST_EQUALITY_CONST( rkbt->numStages(), 3 );
  const Teuchos::SerialDenseMatrix<int,double> A = rkbt->A();
  const Teuchos::SerialDenseVector<int,double> b = rkbt->b();
  const Teuchos::SerialDenseVector<int,double> c = rkbt->c();
  double pi = 4.0*atan(1.0);
  double gamma = (1.0/sqrt(3.0))*cos(pi/18.0)+0.5;
  double delta = 1.0/(6.0*pow(2.0*gamma-1.0,2.0));
  TEST_FLOATING_EQUALITY( A(0,0),  gamma, tol );
  TEST_FLOATING_EQUALITY( A(0,1),  0.0, tol );
  TEST_FLOATING_EQUALITY( A(0,2),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(1,0),  0.5-gamma, tol );
  TEST_FLOATING_EQUALITY( A(1,1),  gamma, tol );
  TEST_FLOATING_EQUALITY( A(1,2),  0.0, tol );

  TEST_FLOATING_EQUALITY( A(2,0),  2.0*gamma, tol );
  TEST_FLOATING_EQUALITY( A(2,1),  1.0-4.0*gamma, tol );
  TEST_FLOATING_EQUALITY( A(2,2),  gamma, tol );

  TEST_FLOATING_EQUALITY( b(0), delta, tol );
  TEST_FLOATING_EQUALITY( b(1), 1.0-2.0*delta, tol );
  TEST_FLOATING_EQUALITY( b(2), delta, tol );

  TEST_FLOATING_EQUALITY( c(0), gamma, tol );
  TEST_FLOATING_EQUALITY( c(1), 0.5, tol );
  TEST_FLOATING_EQUALITY( c(2), 1.0-gamma, tol );

  TEST_EQUALITY_CONST( rkbt->order(), 4 );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validateDIRKButcherTableau ) {
  {
    // Entries above the diagonal should throw
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,2) = 1.0;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(2) = 1.0;
    Teuchos::SerialDenseVector<int,double> c(numStages);
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW(validateDIRKButcherTableau(*rkbt),std::logic_error);
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
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,1);
    validateDIRKButcherTableau(*rkbt);
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
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW(validateSDIRKButcherTableau(*rkbt),std::logic_error);
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
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,1);
    TEST_THROW(validateSDIRKButcherTableau(*rkbt),std::logic_error);
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
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,1);
    validateSDIRKButcherTableau(*rkbt);
  }
}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableauFactory, validParameters ) {
//  TEST_EQUALITY_CONST( true, false );
//}

//TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableauFactory, create ) {
//  TEST_EQUALITY_CONST( true, false );
//}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, operatorEqualEqual ) {
  RCP<RKButcherTableauBuilder<double> > rkbtB = rKButcherTableauBuilder<double>();
  {
    RCP<RKButcherTableauBase<double> > rkbt_1 = rkbtB->create("Explicit 4 Stage");
    {
      RCP<RKButcherTableauBase<double> > rkbt_2 = rkbtB->create("Explicit 4 Stage");
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, true ); // same tableau
    }
    {
      RCP<RKButcherTableauBase<double> > rkbt_2 = rkbtB->create("Explicit 3/8 Rule");
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // different values in A,b,c
    }
    {
      RCP<RKButcherTableauBase<double> > rkbt_2 = rkbtB->create("Backward Euler");
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // different number of stages
    }
    {
      RCP<RKButcherTableauBase<double> > rkbt_2 = rkbtB->create("Implicit 3 Stage 6th order Gauss");
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // different number of stages
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
    RCP<RKButcherTableauBase<double> > rkbt_1 = rKButcherTableau<double>(A,b,c,order);
    {
      Teuchos::SerialDenseMatrix<int,double> Aprime(numStages,numStages);
      Aprime(0,0) = 1.0/1.5;
      Aprime(0,1) = 0.5;
      Aprime(1,0) = 0.5;
      Aprime(1,1) = 1.0;
      RCP<RKButcherTableauBase<double> > rkbt_2 = rKButcherTableau<double>(Aprime,b,c,order);
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, true ); // Same values in A
    }
    {
      Teuchos::SerialDenseMatrix<int,double> Aprime(numStages,numStages);
      Aprime(0,0) = 2.0/3.0;
      Aprime(0,1) = 0.5;
      Aprime(1,0) = 0.5;
      Aprime(1,1) = 1.1;
      RCP<RKButcherTableauBase<double> > rkbt_2 = rKButcherTableau<double>(Aprime,b,c,order);
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // differing values in A
    }
    {
      Teuchos::SerialDenseVector<int,double> bprime(numStages);
      bprime(0) = 0.0;
      bprime(1) = 0.9;
      RCP<RKButcherTableauBase<double> > rkbt_2 = rKButcherTableau<double>(A,bprime,c,order);
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // differing values in b
    }
    {
      Teuchos::SerialDenseVector<int,double> cprime(numStages);
      cprime(0) = 0.0;
      cprime(1) = 0.9;
      RCP<RKButcherTableauBase<double> > rkbt_2 = rKButcherTableau<double>(A,b,cprime,order);
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // differing values in c
    }
    {
      int orderprime = 3;
      RCP<RKButcherTableauBase<double> > rkbt_2 = rKButcherTableau<double>(A,b,c,orderprime);
      TEST_EQUALITY_CONST( *rkbt_1 == *rkbt_2, false ); // differing values in order
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, determineType ) {
  RCP<RKButcherTableauBuilder<double> > rkbtB = rKButcherTableauBuilder<double>();
  { 
    //RCP<RKButcherTableauBase<double> > rkbt = rcp(new RKButcherTableau<double>());
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>();
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_INVALID );
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt  = rkbtB->create("Explicit 4 Stage");
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_ERK );
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt  = rkbtB->create("Forward Euler");
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_ERK );
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt  = rkbtB->create("Backward Euler");
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK );
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt  = rkbtB->create("Implicit 2 Stage 4th order Gauss");
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_IRK );
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt  = rkbtB->create("Singly Diagonal IRK 5 Stage 5th order");
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_SDIRK );
  }
  {
    int numStages = 3;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    Teuchos::SerialDenseVector<int,double> b(numStages);
    Teuchos::SerialDenseVector<int,double> c(numStages);
    A(0,0) = 1.0;
    A(1,1) = 2.0;
    A(2,2) = 3.0;
    b(0) = 0.0;
    b(1) = 0.0;
    b(2) = 1.0;
    c(0) = 0.0;
    c(1) = 0.5;
    c(2) = 1.0;
    RCP<RKButcherTableauBase<double> > rkbt = rKButcherTableau<double>(A,b,c,2);
    E_RKButcherTableauTypes type = determineRKBTType(*rkbt);
    TEST_EQUALITY_CONST( type, RYTHMOS_RK_BUTCHER_TABLEAU_TYPE_DIRK );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, createRKBT ) {
  RCP<RKButcherTableauBase<double> > rkbtA = createRKBT<double>("Forward Euler");
  RCP<ForwardEuler_RKBT<double> > rkbtB = Teuchos::rcp_dynamic_cast<ForwardEuler_RKBT<double> >(rkbtA,false);
  TEST_ASSERT( !is_null(rkbtB) );
}

/*
TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validateERKOrder3 ) {
  // This is only implemented for order=3 at the moment.
  {
    RCP<RKButcherTableauBase<double> > rkbt = createExplicit4Stage3rdOrderRunge_RKBT<double>();
    TEST_THROW(validateERKOrder<double>(*rkbt,1), std::logic_error);
    TEST_THROW(validateERKOrder<double>(*rkbt,2), std::logic_error);
               validateERKOrder<double>(*rkbt,3);
    TEST_THROW(validateERKOrder<double>(*rkbt,4), std::logic_error);
  }
  { 
    RCP<RKButcherTableauBase<double> > rkbt = createExplicit4Stage4thOrder_RKBT<double>();
    TEST_THROW(validateERKOrder<double>(*rkbt,1), std::logic_error);
    TEST_THROW(validateERKOrder<double>(*rkbt,2), std::logic_error);
    TEST_THROW(validateERKOrder<double>(*rkbt,3), std::logic_error);
    TEST_THROW(validateERKOrder<double>(*rkbt,4), std::logic_error);
  }
  {
    RCP<RKButcherTableauBase<double> > rkbt;
    TEST_THROW(validateERKOrder<double>(*rkbt,1), std::logic_error);
  }
}
*/

// This test is just for diagnostic purposes to see the valid parameter list
/*
#ifdef RYTHMOS_DEBUG
TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, FactoryGetValidPL ) {
  RCP<RKButcherTableauBuilder<double> > rkbtF = rKButcherTableauBuilder<double>();
  RCP<const ParameterList> pl = rkbtF->getValidParameters();
  std::cout << "Valid Parameter List for RKButcherTableauBuilder:" << std::endl;
  pl->print(std::cout,Teuchos::ParameterList::PrintOptions().showDoc(true).indent(4));
  TEST_ASSERT( true );
}
#endif // RYTHMOS_DEBUG
*/

} // namespace Rythmos

