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

#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validInitialize ) {
  int numStages = 1;
  Teuchos::SerialDenseMatrix<int,double> A(numStages, numStages);
  Teuchos::SerialDenseVector<int,double> b(numStages);
  Teuchos::SerialDenseVector<int,double> c(numStages);
  A(0,0) = 1.0;
  b(0) = 1.0;
  c(0) = 1.0;
  RKButcherTableau<double> rkButcherTableau(A,b,c);

  TEST_EQUALITY( rkButcherTableau.numStages(), numStages );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidInitialize ) {

  RKButcherTableau<double> rkButcherTableau;

  int numStages = 1;
  {
    Teuchos::SerialDenseMatrix<int,double> A(numStages+1, numStages);
    Teuchos::SerialDenseVector<int,double> b(numStages);
    Teuchos::SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    Teuchos::SerialDenseMatrix<int,double> A(numStages, numStages+1);
    Teuchos::SerialDenseVector<int,double> b(numStages);
    Teuchos::SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    Teuchos::SerialDenseMatrix<int,double> A(numStages, numStages);
    Teuchos::SerialDenseVector<int,double> b(numStages+1);
    Teuchos::SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    Teuchos::SerialDenseMatrix<int,double> A(numStages, numStages);
    Teuchos::SerialDenseVector<int,double> b(numStages);
    Teuchos::SerialDenseVector<int,double> c(numStages+1);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidAssembleIRKState ) {
  int N = 10;
  int numStages = 1;

  {
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, &*x), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages+1,numStages);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, &*x), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages+1);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, &*x), 
        std::logic_error 
        );
  }
  {
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages+1,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, &*x), 
        std::logic_error 
        );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleIRKState ) {

  {
    int N = 1;
    int numStages = 1;
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 5.0;
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,1.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,2.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,3.0);
    assembleIRKState(stageIndex, A, dt, *x_base, *x_stage_bar, &*x);
    // What should x be?
    // x = x_base == 1.0      
    // x += dt*A(0,0)*x_stage_bar(0), so x = 1 + 0.1*5.0*2.0 = 2.0
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 2.0 );
  }
  {
    int N = 1;
    int numStages = 2;
    int stageIndex = 0;
    Teuchos::SerialDenseMatrix<int,double> A(numStages,numStages);
    A(0,0) = 1.0;
    A(0,1) = 2.0;
    A(1,0) = 3.0;
    A(1,1) = 4.0;
    double dt = 10.0;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,5.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,6.0);
    V_S(&*(x_stage_bar->getNonconstVectorBlock(1)),7.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,8.0);
    assembleIRKState(stageIndex, A, dt, *x_base, *x_stage_bar, &*x);
    // What should x be?
    // x = x_base == 5.0               so x = 5
    // x += dt*A(0,0)*x_stage_bar(0)   so x = 5 + 10.0*1.0*6.0 = 65
    // x += dt*A(0,1)*x_stage_bar(1)   so x = 65 + 10.0*2.0*7.0 = 205
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 205 );

    assembleIRKState(stageIndex+1, A, dt, *x_base, *x_stage_bar, &*x);
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
    Teuchos::SerialDenseVector<int,double> b(numStages+1);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKSolution(b, dt, *x_base, *x_stage_bar, &*x),
        std::logic_error 
        );
  }
  {
    Teuchos::SerialDenseVector<int,double> b(numStages);
    double dt = 0.1;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,0.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages+1,N,0.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,0.0);
    
    TEST_THROW( 
        assembleIRKSolution(b, dt, *x_base, *x_stage_bar, &*x),
        std::logic_error 
        );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, assembleIRKSolution ) {

  {
    int N = 1;
    int numStages = 1;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(0) = 2.0;
    double dt = 10.0;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,3.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,4.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,5.0);
    assembleIRKSolution(b, dt, *x_base, *x_stage_bar, &*x);
    // What should x be?
    // x = x_base == 3.0             so x = 3
    // x += dt*b(0)*x_stage_bar(0)   so x = 3 + 10.0*2.0*4.0 = 83
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 83 );
  }
  {
    int N = 1;
    int numStages = 2;
    Teuchos::SerialDenseVector<int,double> b(numStages);
    b(0) = 2.0;
    b(1) = 3.0;
    double dt = 10.0;
    Teuchos::RCP<Thyra::VectorBase<double> > x_base = createDefaultVector(N,4.0);
    Teuchos::RCP<Thyra::ProductVectorBase<double> > x_stage_bar = createDefaultProductVector(numStages,N,5.0);
    V_S(&*(x_stage_bar->getNonconstVectorBlock(1)),6.0);
    Teuchos::RCP<Thyra::VectorBase<double> > x = createDefaultVector(N,7.0);
    assembleIRKSolution(b, dt, *x_base, *x_stage_bar, &*x);
    // What should x be?
    // x = x_base == 4.0             so x = 4
    // x += dt*b(0)*x_stage_bar(0)   so x = 4 + 10.0*2.0*5.0 = 104
    // x += dt*b(1)*x_stage_bar(1)   so x = 104 + 10.0*3.0*6.0 = 284
    TEST_EQUALITY_CONST( Thyra::get_ele(*x,0), 284 );
  }
}

} // namespace Rythmos

