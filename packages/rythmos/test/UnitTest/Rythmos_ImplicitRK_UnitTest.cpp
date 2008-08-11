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
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_ImplicitRKStepper.hpp"
#include "Rythmos_UnitTestModels.hpp"

namespace Rythmos {

using Teuchos::SerialDenseMatrix;
using Teuchos::SerialDenseVector;
using Thyra::VectorBase;
using Thyra::ProductVectorBase;
using Thyra::VectorSpaceBase;
using Thyra::ProductVectorSpaceBase;

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, validInitialize ) {
  int numStages = 1;
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  RKButcherTableau<double> rkButcherTableau(A,b,c);

  TEST_EQUALITY( rkButcherTableau.numStages(), numStages );
}

TEUCHOS_UNIT_TEST( Rythmos_RKButcherTableau, invalidInitialize ) {

  RKButcherTableau<double> rkButcherTableau;

  int numStages = 1;
  {
    SerialDenseMatrix<int,double> A(numStages+1, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages+1);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages+1);
    SerialDenseVector<int,double> c(numStages);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
  }
  {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages+1);
    TEST_THROW( rkButcherTableau = RKButcherTableau<double>(A,b,c), std::logic_error );
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

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, emptyCreate ) {
  ImplicitRKModelEvaluator<double> irkME;
  TEST_EQUALITY_CONST( irkME.get_x_space(), Teuchos::null );
  TEST_EQUALITY_CONST( irkME.get_f_space(), Teuchos::null );
  TEST_THROW( irkME.create_W_op(), std::logic_error );
  TEST_EQUALITY_CONST( irkME.get_W_factory(), Teuchos::null );

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<double> nominalValues = irkME.getNominalValues();
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_t), false );

  MEB::InArgs<double> inArgs = irkME.createInArgs();
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), false );

  RCP<VectorBase<double> > x_old;
  double t_old = 0.1;
  double delta_t = -0.1;
  TEST_THROW(
      irkME.setTimeStepPoint( x_old, t_old, delta_t ),
      std::logic_error
      );
  RCP<VectorBase<double> > valid_x_old = createDefaultVector<double>(10, 0.1);
  TEST_THROW(
      irkME.setTimeStepPoint( valid_x_old, t_old, delta_t ),
      std::logic_error
      );
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, emptyInitialize ) {
  RCP<Thyra::ModelEvaluator<double> > daeModel;
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory;
  RKButcherTableau<double> irkButcherTableau;

  ImplicitRKModelEvaluator<double> irkME;
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = irkME.createInArgs();
  // The empty RCPs will result in an exception
  TEST_THROW(
    irkME.initializeIRKModel(
        daeModel, basePoint, irk_W_factory, irkButcherTableau 
        ),
    std::logic_error
    );
  RCP<ImplicitRKModelEvaluator<double> > irkMEptr;
  TEST_THROW(
    irkMEptr = implicitRKModelEvaluator<double>( 
        daeModel, basePoint, irk_W_factory, irkButcherTableau 
        ),
    std::logic_error
    );
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, validInitialize ) {
  // create valid model with W
  RCP<ParameterList> paramList = Teuchos::parameterList();
  sublist(paramList,Stratimikos_name);
  sublist(paramList,DiagonalTransientModel_name);
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  
  // create an empty but valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > x = Thyra::createMember(model->get_x_space());
  basePoint.set_x(x);
  basePoint.set_t(0.5);

  // create valid W_factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  
  // create valid RKButcherTableau 
  for (int numStages = 1; numStages <= 5 ; ++numStages) {
    SerialDenseMatrix<int,double> A(numStages, numStages);
    SerialDenseVector<int,double> b(numStages);
    SerialDenseVector<int,double> c(numStages);
    RKButcherTableau<double> irkButcherTableau(A,b,c);

    // initialize irk model evaluator
    RCP<ImplicitRKModelEvaluator<double> > irkME = 
      implicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );
    
    // Check that x_space is a product vector space of the correct dimension.
    RCP<const VectorSpaceBase<double> > x_space = irkME->get_x_space();
    RCP<const ProductVectorSpaceBase<double> > x_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(x_space);
    TEST_EQUALITY( x_pspace->numBlocks(), numStages );

    RCP<const VectorSpaceBase<double> > f_space = irkME->get_f_space();
    RCP<const ProductVectorSpaceBase<double> > f_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(f_space);
    TEST_EQUALITY( f_pspace->numBlocks(), numStages );
    
    // Check get_W_factory returns the same factory we passed in
    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > out_irk_W_factory = irkME->get_W_factory();
    TEST_EQUALITY( out_irk_W_factory, irk_W_factory );
    
    // Check create_W_op creates a block matrix of model->create_W_op objects with dimension (stages,stages)
    RCP<Thyra::LinearOpBase<double> > irk_W_op = irkME->create_W_op();
    RCP<Thyra::BlockedLinearOpBase<double> > block_irk_W_op = 
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(irk_W_op);
    RCP<const ProductVectorSpaceBase<double> > block_op_range =
      block_irk_W_op->productRange();
    TEST_EQUALITY( block_op_range->numBlocks(), numStages );
    RCP<const ProductVectorSpaceBase<double> > block_op_domain =
      block_irk_W_op->productDomain();
    TEST_EQUALITY( block_op_domain->numBlocks(), numStages );
    for (int i=0 ; i<numStages ; ++i) {
      for (int j=0 ; j<numStages ; ++j) {
      RCP<const Thyra::LinearOpBase<double> > subBlock =
        block_irk_W_op->getBlock(i,j);
      RCP<const VectorSpaceBase<double> > subBlockRange =
        subBlock->range();
      TEST_EQUALITY_CONST( subBlockRange->isCompatible(*(model->get_f_space())), true );
      RCP<const VectorSpaceBase<double> > subBlockDomain =
        subBlock->domain();
      TEST_EQUALITY_CONST( subBlockDomain->isCompatible(*(model->get_x_space())), true );
      }
    }

    // Call setTimeStepPoint with product vector and check that it throws an exception
    // setTimeStepPoint is supposed to get a vector from the underlying DAE model
    // vector space, NOT a product vector.
    RCP<const ProductVectorSpaceBase<double> > large_pvs = 
      Thyra::productVectorSpace( x_pspace->getBlock(0), numStages+1 );
    RCP<VectorBase<double> > x_old = createMember( *large_pvs );
    double t_old = numStages + 0.123;
    double delta_t = numStages + 0.321;
    TEST_THROW( irkME->setTimeStepPoint( x_old, t_old, delta_t ), std::logic_error );
    if (numStages > 2) { // numStages = 2 means one less is the correct size vector space.
      RCP<const ProductVectorSpaceBase<double> > large_pvs = 
        Thyra::productVectorSpace( x_pspace->getBlock(0), numStages-1 );
      RCP<VectorBase<double> > x_old = createMember( *large_pvs );
      double t_old = numStages + 0.123;
      double delta_t = numStages + 0.321;
      TEST_THROW( irkME->setTimeStepPoint( x_old, t_old, delta_t ), std::logic_error );
    }
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, evalModel ) {
  
  // Create a Model Evaluator
  RCP<ParameterList> paramList = Teuchos::parameterList();
  sublist(paramList,Stratimikos_name);
  sublist(paramList,DiagonalTransientModel_name);
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  TEST_EQUALITY_CONST( is_null(model), false );
  
  // create an empty but valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x), false );
  V_S(&*base_x,2.0);
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x_dot), false );
  V_S(&*base_x_dot,3.0);
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  // create valid W_factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  TEST_EQUALITY_CONST( is_null(irk_W_factory), false );
  
  int numStages = 1;
  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  // Backward Euler
  A(0,0) = 1.0;
  b(0) = 1.0;
  c(0) = 1.0;
  RKButcherTableau<double> irkButcherTableau(A,b,c);

  // Create an IRKModelEvaluator
  RCP<ImplicitRKModelEvaluator<double> > irkME = 
    implicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );

  // Create IRKModelEvaluator InArgs and OutArgs
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = irkME->createInArgs();
  RCP<VectorBase<double> > x_in = createMember(*(irkME->get_x_space()));
  TEST_EQUALITY_CONST( is_null(x_in), false );
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = irkME->createOutArgs();
  RCP<VectorBase<double> > f_out = createMember(*(irkME->get_f_space()));
  TEST_EQUALITY_CONST( is_null(f_out), false );
  RCP<Thyra::LinearOpBase<double> > W_op_out = irkME->create_W_op();
  TEST_EQUALITY_CONST( is_null(W_op_out), false );
  outArgs.set_f(f_out);
  outArgs.set_W_op(W_op_out);

  // Fill x_in with valid data.
  // Note:  These are product vectors in the numStages>1 case
  V_S(outArg(*x_in),5.0);
  inArgs.set_x(x_in);
  
  // We need to call setTimeStepPoint before evalModel
  TEST_THROW( irkME->evalModel(inArgs,outArgs), std::logic_error );
      
  double delta_t = 6.0;
  irkME->setTimeStepPoint( base_x, base_t, delta_t );
  irkME->evalModel(inArgs,outArgs);

  // Verify contents of f_out and W_op_out.
  // TODO

}

} // namespace Rythmos

