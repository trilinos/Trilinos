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
#include "Rythmos_RKButcherTableau.hpp"
#include "Rythmos_RKButcherTableauBuilder.hpp"
#include "Rythmos_SingleResidualModelEvaluator.hpp"
#include "Rythmos_ImplicitRKModelEvaluator.hpp"
#include "Rythmos_DiagonalImplicitRKModelEvaluator.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

namespace Rythmos {

using Teuchos::SerialDenseMatrix;
using Teuchos::SerialDenseVector;
using Thyra::VectorBase;
using Thyra::ProductVectorBase;
using Thyra::VectorSpaceBase;
using Thyra::ProductVectorSpaceBase;

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, emptyCreate ) {
  ImplicitRKModelEvaluator<double> irkME;
  TEST_THROW( irkME.get_x_space(), std::logic_error );
  TEST_THROW( irkME.get_f_space(), std::logic_error );
  TEST_THROW( irkME.create_W_op(), std::logic_error );
  TEST_THROW( irkME.get_W_factory(), std::logic_error );

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
  RCP<RKButcherTableauBase<double> > irkButcherTableau;

  ImplicitRKModelEvaluator<double> irkME;
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint; 
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
    RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,numStages);

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
    {
      RCP<const ProductVectorSpaceBase<double> > large_pvs = 
        Thyra::productVectorSpace( x_pspace->getBlock(0), numStages+1 );
      RCP<VectorBase<double> > x_old = createMember( *large_pvs );
      double t_old = numStages + 0.123;
      double delta_t = numStages + 0.321;
      TEST_THROW( irkME->setTimeStepPoint( x_old, t_old, delta_t ), std::logic_error );
    }
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

RCP<ImplicitRKModelEvaluator<double> > getImplicitRKModelEvaluator(int numStages) {
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
  
  // create a valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  V_S(base_x.ptr(),2.0);
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  V_S(base_x_dot.ptr(),3.0);
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,1);

  // Create an IRKModelEvaluator
  RCP<ImplicitRKModelEvaluator<double> > irkME = 
    implicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );

  return irkME;
}


TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, evalModelOneStage ) {
  double tol = 1.0e-10;
  
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
//  RCP<ParameterList> paramList = Teuchos::parameterList();
//  sublist(paramList,Stratimikos_name);
//  sublist(paramList,DiagonalTransientModel_name);
//  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  TEST_EQUALITY_CONST( is_null(model), false );
  
  // create an valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x), false );
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 2.0;
    base_x_view[1] = 2.5;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x_dot), false );
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 3.0;
    base_x_dot_view[1] = 3.5;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  // create valid W_factory
//  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  TEST_EQUALITY_CONST( is_null(irk_W_factory), false );
  
  int numStages = 1;
  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  // Backward Euler
  A(0,0) = 7.0;
  b(0) = 8.0;
  c(0) = 9.0;
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,1);

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
  RCP<Thyra::ProductVectorBase<double> > x_in_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(x_in,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(x_in_b), false);
  RCP<VectorBase<double> > x_in_b0 = x_in_b->getNonconstVectorBlock(0);
  {
    Thyra::DetachedVectorView<double> x_in_b0_view( *x_in_b0 );
    x_in_b0_view[0] = 5.0;
    x_in_b0_view[1] = 5.5;
  }
  inArgs.set_x(x_in);
  
  // We need to call setTimeStepPoint before evalModel
  TEST_THROW( irkME->evalModel(inArgs,outArgs), std::logic_error );
      
  double delta_t = 6.0;
  irkME->setTimeStepPoint( base_x, base_t, delta_t );
  irkME->evalModel(inArgs,outArgs);

  // Verify contents of f_out and W_op_out.
  // 
  // What should f_out be?
  // F(t,x,xdot)=0, this is from Petzold/Brennan/Campbell.
  // h = delta_t
  // M = numStages
  // F_i = F( t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^' )
  // For numSages = 1:
  // F_1 = F( t_{n-1}+c*h, base_x + h*a_{1,1}*Y', Y' )
  // for A = 7, c = 9:
  // F_1 = F( t_{n-1}+9*h, base_x + 7*h*Y', Y' )
  // The SinCos Model's implicit formulation is:
  // f(t,x,xdot) (0) = xdot(0) - x(1)
  // f(t,x,xdot) (1) = xdot(1) + x(0)
  // F_1_exact (0) = Y'(0) - (base_x(1) + 7*h*Y'(1))
  // F_1_exact (1) = Y'(1) + (base_x(0) + 7*h*Y'(0))
  // base_x = (2.0,2.5)
  // Y' = (5.0,5.5)
  // h = 6.0
  // F_1_exact (0) = 5   - (2.5 + 7*6*5.5) = -228.5
  // F_1_exact (1) = 5.5 + (2   + 7*6*5  ) =  217.5
  RCP<Thyra::ProductVectorBase<double> > f_out_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(f_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(f_out_b), false);
  RCP<const VectorBase<double> > f_out_b0 = f_out_b->getVectorBlock(0);
  {
    Thyra::ConstDetachedVectorView<double> f_out_b0_view( *f_out_b0 );
    TEST_FLOATING_EQUALITY( f_out_b0_view[0], -228.5, tol );
    TEST_FLOATING_EQUALITY( f_out_b0_view[1],  217.5, tol );
  }

  // What should W_op_out be?
  // W_op_out_{i,k} = \frac{\partial F_i}{\partial Y_k^'}
  //                =                \frac{\partial F}{\partial x   }(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*h*a_{i,k}
  //                  + \delta_{i==k}\frac{\partial F}{\partial xdot}(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*1
  // For numStages = 1:
  // W_op_out_{1,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+c*h, base_x + h*a*Y', Y')*h*a 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c*h, base_x + h*a*Y', Y')*1
  // For A=a=7, c=9:
  // W_op_out_{1,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+9*h, base_x + 7*h*Y', Y')*h*7
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+9*h, base_x + 7*h*Y', Y')
  // The SinCos Model's implicit formulation's W_op is \alpha*\frac{\partial F}{\partial xdot} + \beta*\frac{\partial F}{\partial x}:
  // W_op (0,0) = alpha
  // W_op (0,1) = -beta
  // W_op (1,0) =  beta
  // W_op (1,1) = alpha
  // In the IRK method:  alpha = 1.0, and beta = a*h = 7*6.0
  // W_op_out_exact (0,0) = alpha   =   1.0
  // W_op_out_exact (0,1) = -beta*h = -42.0
  // W_op_out_exact (1,0) =  beta*h =  42.0
  // W_op_out_exact (1,1) = alpha   =   1.0
  RCP<Thyra::BlockedLinearOpBase<double> > W_op_out_b = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b), false);
  RCP<const Thyra::LinearOpBase<double> > W_op_out_b0 = W_op_out_b->getBlock(0,0);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_b0_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out_b0,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b0_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_b0_mv_view( *W_op_out_b0_mv );
    TEST_FLOATING_EQUALITY( W_op_out_b0_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b0_mv_view(0,1), -42.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b0_mv_view(1,0),  42.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b0_mv_view(1,1),   1.0, tol );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, evalModelTwoStage ) {
  double tol=1.0e-10; 
  
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
//  RCP<ParameterList> paramList = Teuchos::parameterList();
//  sublist(paramList,Stratimikos_name);
//  sublist(paramList,DiagonalTransientModel_name);
//  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  TEST_EQUALITY_CONST( is_null(model), false );
  
  // create an valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x), false );
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 2.0;
    base_x_view[1] = 2.5;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x_dot), false );
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 3.0;
    base_x_dot_view[1] = 3.5;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  // create valid W_factory
//  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  TEST_EQUALITY_CONST( is_null(irk_W_factory), false );
  
  int numStages = 2;
  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  // Backward Euler
  A(0,0) = 7.1;
  A(0,1) = 7.2;
  A(1,0) = 7.3;
  A(1,1) = 7.4;
  b(0) = 8.0;
  b(1) = 8.1;
  c(0) = 9.0;
  c(1) = 9.1;
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,2);

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
  RCP<Thyra::ProductVectorBase<double> > x_in_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(x_in,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(x_in_b), false);
  RCP<VectorBase<double> > x_in_b0 = x_in_b->getNonconstVectorBlock(0);
  {
    Thyra::DetachedVectorView<double> x_in_b0_view( *x_in_b0 );
    x_in_b0_view[0] = 5.1;
    x_in_b0_view[1] = 5.2;
  }
  RCP<VectorBase<double> > x_in_b1 = x_in_b->getNonconstVectorBlock(1);
  {
    Thyra::DetachedVectorView<double> x_in_b1_view( *x_in_b1 );
    x_in_b1_view[0] = 5.6;
    x_in_b1_view[1] = 5.7;
  }

  inArgs.set_x(x_in);
  
  // We need to call setTimeStepPoint before evalModel
  TEST_THROW( irkME->evalModel(inArgs,outArgs), std::logic_error );
      
  double delta_t = 6.0;
  irkME->setTimeStepPoint( base_x, base_t, delta_t );
  irkME->evalModel(inArgs,outArgs);

  // Verify contents of f_out and W_op_out.
  // 
  // What should f_out be?
  // F(t,x,xdot)=0, this is from Petzold/Brennan/Campbell.
  // h = delta_t
  // M = numStages
  // F_i = F( t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^' )
  // For numSages = 2:
  // F_1 = F( t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'+a_{1,2}*Y_2^'), Y_1^' )
  // F_2 = F( t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^' )
  // for our specific A and c:
  // F_1 = F( t_{n-1}+9.0*h, base_x + h*(7.1*Y_1^'+7.2*Y_2^'), Y_1^' )
  // F_2 = F( t_{n-1}+9.1*h, base_x + h*(7.3*Y_1^'+7.4*Y_2^'), Y_2^' )
  // The SinCos Model's implicit formulation is:
  // f(t,x,xdot) (0) = xdot(0) - x(1)
  // f(t,x,xdot) (1) = xdot(1) + x(0)
  // F_1_exact (0) = Y_1^'(0) - (base_x(1) + h*(7.1*Y_1^'(1)+7.2*Y_2^'(1)))
  // F_1_exact (1) = Y_1^'(1) + (base_x(0) + h*(7.1*Y_1^'(0)+7.2*Y_2^'(0)))
  // F_2_exact (0) = Y_2^'(0) - (base_x(1) + h*(7.3*Y_1^'(1)+7.4*Y_2^'(1)))
  // F_2_exact (1) = Y_2^'(1) + (base_x(0) + h*(7.3*Y_1^'(0)+7.4*Y_2^'(0)))
  // base_x = (2.0,2.5)
  // Y' = [ (5.1,5.2), (5.6,5.7) ]
  // h = 6.0
  // F_1_exact (0) = 5.1 - ( 2.5 + 6*(7.1*5.2 + 7.2*5.7 ) ) = -465.16
  // F_1_exact (1) = 5.2 + ( 2.0 + 6*(7.1*5.1 + 7.2*5.6 ) ) =  466.38
  // F_2_exact (0) = 5.6 - ( 2.5 + 6*(7.3*5.2 + 7.4*5.7 ) ) = -477.74
  // F_2_exact (1) = 5.7 + ( 2.0 + 6*(7.3*5.1 + 7.4*5.6 ) ) =  479.72
  //
  RCP<Thyra::ProductVectorBase<double> > f_out_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(f_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(f_out_b), false);
  RCP<const VectorBase<double> > f_out_b0 = f_out_b->getVectorBlock(0);
  {
    Thyra::ConstDetachedVectorView<double> f_out_b0_view( *f_out_b0 );
    TEST_FLOATING_EQUALITY( f_out_b0_view[0], -465.16, tol );
    TEST_FLOATING_EQUALITY( f_out_b0_view[1],  466.38, tol );
  }
  RCP<const VectorBase<double> > f_out_b1 = f_out_b->getVectorBlock(1);
  {
    Thyra::ConstDetachedVectorView<double> f_out_b1_view( *f_out_b1 );
    TEST_FLOATING_EQUALITY( f_out_b1_view[0], -477.74, tol );
    TEST_FLOATING_EQUALITY( f_out_b1_view[1],  479.72, tol );
  }

  // What should W_op_out be?
  // W_op_out_{i,k} = \frac{\partial F_i}{\partial Y_k^'}
  //                =                \frac{\partial F}{\partial x   }(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*h*a_{i,k}
  //                  + \delta_{i==k}\frac{\partial F}{\partial xdot}(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*1
  // For numStages = 2:
  // W_op_out_{1,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'+a_{1,2}*Y_2^'), Y_1^')*h*a_{1,1} 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'+a_{1,2}*Y_2^'), Y_1^')*1
  // W_op_out_{1,2} =   \frac{\partial F}{\partial x   }(t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'+a_{1,2}*Y_2^'), Y_1^')*h*a_{1,2} 
  // W_op_out_{2,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^')*h*a_{2,1} 
  // W_op_out_{2,2} =   \frac{\partial F}{\partial x   }(t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^')*h*a_{2,2} 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^')*1
  //
  // The SinCos Model's implicit formulation's W_op is \alpha*\frac{\partial F}{\partial xdot} + \beta*\frac{\partial F}{\partial x}:
  // W_op (0,0) = alpha
  // W_op (0,1) = -beta
  // W_op (1,0) =  beta
  // W_op (1,1) = alpha
  // In the IRK method:  alpha and beta vary by the block
  //   block (1,1):  alpha = 1  beta = h*a_{1,1} = 6*7.1 = 42.6
  //   block (1,2):  alpha = 0  beta = h*a_{1,2} = 6*7.2 = 43.2
  //   block (2,1):  alpha = 0  beta = h*a_{2,1} = 6*7.3 = 43.8
  //   block (2,2):  alpha = 1  beta = h*a_{2,2} = 6*7.4 = 44.4
  //
  // W_op_out_exact_{1,1} (0,0) = alpha =   1.0
  // W_op_out_exact_{1,1} (0,1) = -beta = -42.6
  // W_op_out_exact_{1,1} (1,0) =  beta =  42.6
  // W_op_out_exact_{1,1} (1,1) = alpha =   1.0
  //
  // W_op_out_exact_{1,2} (0,0) = alpha =   0.0
  // W_op_out_exact_{1,2} (0,1) = -beta = -43.2
  // W_op_out_exact_{1,2} (1,0) =  beta =  43.2
  // W_op_out_exact_{1,2} (1,1) = alpha =   0.0
  //
  // W_op_out_exact_{2,1} (0,0) = alpha =   0.0
  // W_op_out_exact_{2,1} (0,1) = -beta = -43.8
  // W_op_out_exact_{2,1} (1,0) =  beta =  43.8
  // W_op_out_exact_{2,1} (1,1) = alpha =   0.0
  //
  // W_op_out_exact_{2,2} (0,0) = alpha =   1.0
  // W_op_out_exact_{2,2} (0,1) = -beta = -44.4
  // W_op_out_exact_{2,2} (1,0) =  beta =  44.4
  // W_op_out_exact_{2,2} (1,1) = alpha =   1.0
  RCP<Thyra::BlockedLinearOpBase<double> > W_op_out_b = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b), false);
  RCP<const Thyra::LinearOpBase<double> > W_op_out_b00 = W_op_out_b->getBlock(0,0);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_b00_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out_b00,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b00_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_b00_mv_view( *W_op_out_b00_mv );
    TEST_FLOATING_EQUALITY( W_op_out_b00_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b00_mv_view(0,1), -42.6, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b00_mv_view(1,0),  42.6, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b00_mv_view(1,1),   1.0, tol );
  }
  RCP<const Thyra::LinearOpBase<double> > W_op_out_b01 = W_op_out_b->getBlock(0,1);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_b01_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out_b01,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b01_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_b01_mv_view( *W_op_out_b01_mv );
    TEST_FLOATING_EQUALITY( W_op_out_b01_mv_view(0,0),   0.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b01_mv_view(0,1), -43.2, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b01_mv_view(1,0),  43.2, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b01_mv_view(1,1),   0.0, tol );
  }
  RCP<const Thyra::LinearOpBase<double> > W_op_out_b10 = W_op_out_b->getBlock(1,0);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_b10_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out_b10,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b10_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_b10_mv_view( *W_op_out_b10_mv );
    TEST_FLOATING_EQUALITY( W_op_out_b10_mv_view(0,0),   0.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b10_mv_view(0,1), -43.8, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b10_mv_view(1,0),  43.8, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b10_mv_view(1,1),   0.0, tol );
  }
  RCP<const Thyra::LinearOpBase<double> > W_op_out_b11 = W_op_out_b->getBlock(1,1);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_b11_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out_b11,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b11_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_b11_mv_view( *W_op_out_b11_mv );
    TEST_FLOATING_EQUALITY( W_op_out_b11_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b11_mv_view(0,1), -44.4, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b11_mv_view(1,0),  44.4, tol );
    TEST_FLOATING_EQUALITY( W_op_out_b11_mv_view(1,1),   1.0, tol );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, createInArgs ) {
  RCP<ImplicitRKModelEvaluator<double> > irkME = getImplicitRKModelEvaluator(1);
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<double> inArgs = irkME->createInArgs();
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), false );
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, createOutArgs ) {
  RCP<ImplicitRKModelEvaluator<double> > irkME = getImplicitRKModelEvaluator(1);
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<double> outArgs = irkME->createOutArgs();
  TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
  TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, spaces ) {
  for (int numStages=1 ; numStages <= 5 ; ++numStages) {
    RCP<ImplicitRKModelEvaluator<double> > irkME = getImplicitRKModelEvaluator(numStages);
    RCP<const VectorSpaceBase<double> > x_space = irkME->get_x_space();
    RCP<const Thyra::ProductVectorSpaceBase<double> > x_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( Teuchos::is_null(x_pspace), false );
    TEST_EQUALITY( x_pspace->numBlocks(), numStages );

    RCP<const VectorSpaceBase<double> > f_space = irkME->get_f_space();
    RCP<const Thyra::ProductVectorSpaceBase<double> > f_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(f_space,false);
    TEST_EQUALITY_CONST( Teuchos::is_null(f_pspace), false );
    TEST_EQUALITY( f_pspace->numBlocks(), numStages );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, nominalValues ) {
  int numStages = 2;
  RCP<ImplicitRKModelEvaluator<double> > irkME = getImplicitRKModelEvaluator(numStages);
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<double> nominalValues = irkME->getNominalValues();
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_t), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x), true );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_alpha), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_beta), false );
//  RCP<const Thyra::VectorBase<double> > x_nom = nominalValues->get_x();
//  RCP<const Thyra::VectorSpaceBase<double> > x_nom_space = x_nom->space();
//  RCP<const Thyra::ProductVectorSpaceBase<double> > x_nom_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(x_nom_space,false);
//  TEST_EQUALITY_CONST( Teuchos::is_null(x_nom_pspace), false);
//  TEST_EQUALITY( x_nom_pspace->numBlocks(), numStages );

//  RCP<const Thyra::ProductVectorBase<double> > x_nom_b = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double> >(x_nom,false);
//  TEST_EQUALITY_CONST( Teuchos::is_null(x_nom_b), false );

  // What is a valid nominal value for this ME?

}

TEUCHOS_UNIT_TEST( Rythmos_IRKModelEvaluator, create_W_op ) {
  int numStages = 3;
  RCP<ImplicitRKModelEvaluator<double> > irkME = getImplicitRKModelEvaluator(numStages);
  RCP<Thyra::LinearOpBase<double> > W_op = irkME->create_W_op();
  RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > W_op_b = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(W_op,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_b), false );
  TEST_EQUALITY( W_op_b->productRange()->numBlocks(), numStages );
  TEST_EQUALITY( W_op_b->productDomain()->numBlocks(), numStages );
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, create ) {
  RCP<ImplicitRKStepper<double> > irkStepper = rcp(new ImplicitRKStepper<double> );
  TEST_EQUALITY_CONST( Teuchos::is_null(irkStepper), false );
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, invalidRKBT ) {
  RCP<ImplicitRKStepper<double> > stepper = implicitRKStepper<double>();
  RCP<RKButcherTableauBase<double> > rkbt;
  TEST_THROW( stepper->setRKButcherTableau(rkbt), std::logic_error ); // empty RKBT
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, assertValidModel ) {
  {
    // implicit model, OK
    RCP<SinCosModel> model = sinCosModel(true);
    RCP<ImplicitRKStepper<double> > stepper = implicitRKStepper<double>();
    TEST_NOTHROW( stepper->setModel(model) );
  }
  {
    // explicit model, throw
    RCP<SinCosModel> model = sinCosModel(false);
    RCP<ImplicitRKStepper<double> > stepper = implicitRKStepper<double>();
    TEST_THROW( stepper->setModel(model), std::logic_error );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, clone ) {
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  RCP<Teuchos::ParameterList> stratPl = sublist(pl,Stratimikos_name);
  RCP<Teuchos::ParameterList> modelPl = sublist(pl,DiagonalTransientModel_name);
  stratPl->set("Linear Solver Type","AztecOO");
  stratPl->set("Preconditioner Type","None");
  modelPl->set("NumElements",2);
  modelPl->set("Gamma_min",-2.5);
  modelPl->set("Gamma_max",-0.5);
  //RCP<Teuchos::ParameterList> model_voPl = sublist(modelPl,"VerboseObject");
  //model_voPl->set("Verbosity Level","none");
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(pl);
  // Create the nonlinear solver
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  // Create the IRK W factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory =
    getWFactory<double>(pl);
  // Create the RKBT
  RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Implicit 3 Stage 6th order Gauss");
  // Create the IRK Stepper
  RCP<ImplicitRKStepper<double> > stepper = 
    implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );

  TEST_ASSERT( !is_null(stepper) );
  TEST_ASSERT( stepper->supportsCloning() );
  RCP<StepperBase<double> > otherStepper = stepper->cloneStepperAlgorithm();
  TEST_ASSERT( !is_null(otherStepper) );
  TEST_ASSERT( otherStepper.ptr() != stepper.ptr() );
  {
    RCP<ImplicitRKStepper<double> > irkStepper = Teuchos::rcp_dynamic_cast<ImplicitRKStepper<double> >(otherStepper,false);
    TEST_ASSERT( !is_null(irkStepper) );
    RCP<const RKButcherTableauBase<double> > rkbt_out = irkStepper->getRKButcherTableau();
    TEST_ASSERT( !is_null(rkbt_out) );
    TEST_ASSERT( rkbt.ptr() == rkbt_out.ptr() );
    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory_out = irkStepper->get_W_factory();
    TEST_ASSERT( !is_null(irk_W_factory_out) );
    TEST_ASSERT( irk_W_factory.ptr() == irk_W_factory_out.ptr() );
    RCP<const Thyra::NonlinearSolverBase<double> > solver_out = irkStepper->getSolver();
    TEST_ASSERT( !is_null(solver_out) );
    TEST_ASSERT( solver_out.ptr() != nonlinearSolver.ptr() );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, takeStep ) {
  // Create the model
  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  RCP<Teuchos::ParameterList> stratPl = sublist(pl,Stratimikos_name);
  RCP<Teuchos::ParameterList> modelPl = sublist(pl,DiagonalTransientModel_name);
  stratPl->set("Linear Solver Type","AztecOO");
  stratPl->set("Preconditioner Type","None");
  modelPl->set("NumElements",2);
  modelPl->set("Gamma_min",-2.5);
  modelPl->set("Gamma_max",-0.5);
  //RCP<Teuchos::ParameterList> model_voPl = sublist(modelPl,"VerboseObject");
  //model_voPl->set("Verbosity Level","none");
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(pl);
  // create a valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 0.2;
    base_x_view[1] = 0.3;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 0.6;
    base_x_dot_view[1] = 0.7;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 0.5;
  basePoint.set_t(base_t);
  // Create the nonlinear solver
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  // Create the IRK W factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory =
    getWFactory<double>(pl);

  // Create the IRK Butcher Tableau
  RCP<RKButcherTableauBase<double> > irkbt = createRKBT<double>("Backward Euler");
  // Create the IRK Stepper
  RCP<ImplicitRKStepper<double> > irkStepper = 
    implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, irkbt );
  irkStepper->setDirk(false);
  RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
  RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
  stepperVOPL->set("Verbosity Level","none");
  irkStepper->setParameterList(stepperPL);
  irkStepper->setInitialCondition(basePoint);
  double h = 1.2;
  double stepTaken = irkStepper->takeStep(h, STEP_TYPE_FIXED);
  TEST_EQUALITY_CONST( stepTaken, h );
  // What should the answer be?
  // The DE for the DiagonalTransient problem:
  // \dot{x} = lambda*x
  // By setting the parameters above, we get:
  // \lambda_0 = -2.5
  // \lambda_1 = -0.5
  // BE applied to this problem should result in the step:
  // x_n = (1/(1-h*\lambda))*x_{n-1}
  // Therefore, 
  // x_n (0) = x_{n-1}(0)/(1+2.5*h) = 0.2/(1+2.5*1.2) = 0.05
  // x_n (1) = x_{n-1}(1)/(1+0.5*h) = 0.3/(1+0.5*1.2) = 0.1875
  double tol = 1.0e-10;
  RCP<const VectorBase<double> > sol = irkStepper->getStepStatus().solution;
  {
    Thyra::ConstDetachedVectorView<double> sol_view( *sol );
    TEST_FLOATING_EQUALITY( sol_view[0], 0.05, tol );
    TEST_FLOATING_EQUALITY( sol_view[1], 0.1875, tol );
  }
  
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, setDirk ) {
  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>();
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory =
    getWFactory<double>();
  {
    // create stepper with DIRK tableau
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Singly Diagonal IRK 5 Stage 4th order");
    RCP<ImplicitRKStepper<double> > irkStepper = 
      implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );
    irkStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
    RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
    stepperVOPL->set("Verbosity Level","none");
    irkStepper->setParameterList(stepperPL);
    irkStepper->setDirk(true);
    // take step
    irkStepper->takeStep(1.0, STEP_TYPE_FIXED);
    // get_x_space, verify it is NOT a product space
    RCP<const VectorSpaceBase<double> > x_space = irkStepper->getSolver()->getModel()->get_x_space();
    RCP<const ProductVectorSpaceBase<double> > x_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( is_null(x_pspace), true );
    TEST_EQUALITY_CONST( x_space->isCompatible( *model->get_x_space() ), true );
  }
  {
    // create same stepper with same DIRK tableau
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Singly Diagonal IRK 5 Stage 4th order");
    RCP<ImplicitRKStepper<double> > irkStepper = 
      implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );
    irkStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
    RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
    stepperVOPL->set("Verbosity Level","none");
    irkStepper->setParameterList(stepperPL);
    irkStepper->setDirk(false);
    // take step
    irkStepper->takeStep(1.0, STEP_TYPE_FIXED);
    // get_x_space, verfiy it IS a product space of the correct dimension
    RCP<const VectorSpaceBase<double> > x_space = irkStepper->getSolver()->getModel()->get_x_space();
    RCP<const ProductVectorSpaceBase<double> > x_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( is_null(x_pspace), false );
    TEST_EQUALITY( x_pspace->numBlocks(), rkbt->numStages() );
  }
  {
    // create same stepper with same BE tableau
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Backward Euler");
    RCP<ImplicitRKStepper<double> > irkStepper = 
      implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );
    irkStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
    RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
    stepperVOPL->set("Verbosity Level","none");
    irkStepper->setParameterList(stepperPL);
    irkStepper->setDirk(false);
    // take step
    irkStepper->takeStep(1.0, STEP_TYPE_FIXED);
    // get_x_space, verfiy it IS a product space of the correct dimension
    RCP<const VectorSpaceBase<double> > x_space = irkStepper->getSolver()->getModel()->get_x_space();
    RCP<const ProductVectorSpaceBase<double> > x_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( is_null(x_pspace), false );
    TEST_EQUALITY( x_pspace->numBlocks(), rkbt->numStages() );
  }
  {
    // create stepper with fully dense RK Tableau
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Implicit 2 Stage 4th order Gauss");
    RCP<ImplicitRKStepper<double> > irkStepper = 
      implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );
    irkStepper->setInitialCondition(ic);
    TEST_THROW(irkStepper->setDirk(true), std::logic_error);
  }
  {
    // create sepper with any RKBT
    RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Implicit 2 Stage 4th order Gauss");
    RCP<ImplicitRKStepper<double> > irkStepper = 
      implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, rkbt );
    irkStepper->setInitialCondition(ic);
    RCP<Teuchos::ParameterList> stepperPL = Teuchos::parameterList();
    RCP<Teuchos::ParameterList> stepperVOPL = Teuchos::sublist(stepperPL,"VerboseObject");
    stepperVOPL->set("Verbosity Level","none");
    irkStepper->setParameterList(stepperPL);
    // take step
    irkStepper->takeStep(1.0, STEP_TYPE_FIXED);
    TEST_THROW(irkStepper->setDirk(true), std::logic_error)
    TEST_THROW(irkStepper->setDirk(false), std::logic_error)
  }
}


TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, emptyCreate ) {
  DiagonalImplicitRKModelEvaluator<double> dirkME;
  TEST_THROW( dirkME.get_x_space(), std::logic_error );
  TEST_THROW( dirkME.get_f_space(), std::logic_error );
  TEST_THROW( dirkME.create_W_op(), std::logic_error );
  TEST_THROW( dirkME.get_W_factory(), std::logic_error );

  RCP<VectorBase<double> > x_old;
  double t_old = 0.1;
  double delta_t = -0.1;
  TEST_THROW(
      dirkME.setTimeStepPoint( x_old, t_old, delta_t ),
      std::logic_error
      );
  RCP<VectorBase<double> > valid_x_old = createDefaultVector<double>(10, 0.1);
  TEST_THROW(
      dirkME.setTimeStepPoint( valid_x_old, t_old, delta_t ),
      std::logic_error
      );
}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, emptyInitialize ) {
  RCP<Thyra::ModelEvaluator<double> > daeModel;
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory;
  RCP<RKButcherTableauBase<double> > irkButcherTableau;

  DiagonalImplicitRKModelEvaluator<double> dirkME;
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint; 
  // The empty RCPs will result in an exception
  TEST_THROW(
    dirkME.initializeDIRKModel(
        daeModel, basePoint, irk_W_factory, irkButcherTableau 
        ),
    std::logic_error
    );
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkMEptr;
  TEST_THROW(
    dirkMEptr = diagonalImplicitRKModelEvaluator<double>( 
        daeModel, basePoint, irk_W_factory, irkButcherTableau 
        ),
    std::logic_error
    );
}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, validInitialize ) {
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
    RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,numStages);

    // initialize irk model evaluator
    RCP<DiagonalImplicitRKModelEvaluator<double> > irkME = 
      diagonalImplicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );
    
    // In this case, check that x_space is the same space as the underlying model's x_space.
    RCP<const VectorSpaceBase<double> > x_space = irkME->get_x_space();
    RCP<const ProductVectorSpaceBase<double> > x_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( is_null(x_pspace), true );
    TEST_EQUALITY_CONST( x_space->isCompatible( *model->get_x_space() ), true );

    RCP<const VectorSpaceBase<double> > f_space = irkME->get_f_space();
    RCP<const ProductVectorSpaceBase<double> > f_pspace = 
      Teuchos::rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(f_space,false);
    TEST_EQUALITY_CONST( is_null(f_pspace), true );
    TEST_EQUALITY_CONST( f_space->isCompatible( *model->get_f_space() ), true );

    // In this case, check that create_W_op is exactly the same matrix as model->create_W_op
    RCP<Thyra::LinearOpBase<double> > irk_W_op = irkME->create_W_op();
    RCP<Thyra::BlockedLinearOpBase<double> > block_irk_W_op = 
      Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(irk_W_op,false);
    TEST_EQUALITY_CONST( is_null(block_irk_W_op), true );
    TEST_EQUALITY_CONST( irk_W_op->domain()->isCompatible( *model->get_x_space() ), true );
    TEST_EQUALITY_CONST( irk_W_op->range()->isCompatible( *model->get_f_space() ), true );
  }
}

RCP<DiagonalImplicitRKModelEvaluator<double> > getDiagonalImplicitRKModelEvaluator(int numStages) {
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
  
  // create a valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  V_S(base_x.ptr(),2.0);
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  V_S(base_x_dot.ptr(),3.0);
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,1);

  // Create an IRKModelEvaluator
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = 
    diagonalImplicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );

  return dirkME;
}


TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, evalModelOneStage ) {
  double tol = 1.0e-10;
  
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
//  RCP<ParameterList> paramList = Teuchos::parameterList();
//  sublist(paramList,Stratimikos_name);
//  sublist(paramList,DiagonalTransientModel_name);
//  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  TEST_EQUALITY_CONST( is_null(model), false );
  
  // create an valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x), false );
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 2.0;
    base_x_view[1] = 2.5;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x_dot), false );
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 3.0;
    base_x_dot_view[1] = 3.5;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  // create valid W_factory
//  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  TEST_EQUALITY_CONST( is_null(irk_W_factory), false );
  
  int numStages = 1;
  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  // Backward Euler
  A(0,0) = 7.0;
  b(0) = 8.0;
  c(0) = 9.0;
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,1);

  // Create an IRKModelEvaluator
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = 
    diagonalImplicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );

  // Create IRKModelEvaluator InArgs and OutArgs
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = dirkME->createInArgs();
  RCP<VectorBase<double> > x_in = createMember(*(dirkME->get_x_space()));
  TEST_EQUALITY_CONST( is_null(x_in), false );
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = dirkME->createOutArgs();
  RCP<VectorBase<double> > f_out = createMember(*(dirkME->get_f_space()));
  TEST_EQUALITY_CONST( is_null(f_out), false );
  RCP<Thyra::LinearOpBase<double> > W_op_out = dirkME->create_W_op();
  TEST_EQUALITY_CONST( is_null(W_op_out), false );
  outArgs.set_f(f_out);
  outArgs.set_W_op(W_op_out);

  // Fill x_in with valid data.
  // Note:  These are NOT product vectors in the Diagonal Implicit RK case
  RCP<Thyra::ProductVectorBase<double> > x_in_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(x_in,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(x_in_b), true );
  {
    Thyra::DetachedVectorView<double> x_in_view( *x_in );
    x_in_view[0] = 5.0;
    x_in_view[1] = 5.5;
  }
  inArgs.set_x(x_in);
  
  // We need to call setTimeStepPoint before evalModel
  TEST_THROW( dirkME->evalModel(inArgs,outArgs), std::logic_error );
      
  double delta_t = 6.0;
  dirkME->setTimeStepPoint( base_x, base_t, delta_t );
  // We need to call setCurrentStage first!
  TEST_THROW( dirkME->evalModel(inArgs,outArgs), std::logic_error) ;
  dirkME->setCurrentStage(0);
  dirkME->evalModel(inArgs,outArgs);

  // Verify contents of f_out and W_op_out.
  // 
  // What should f_out be?
  // F(t,x,xdot)=0, this is from Petzold/Brennan/Campbell.
  // h = delta_t
  // M = numStages
  // F_i = F( t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^' )
  // For numSages = 1:
  // F_1 = F( t_{n-1}+c*h, base_x + h*a_{1,1}*Y', Y' )
  // for A = 7, c = 9:
  // F_1 = F( t_{n-1}+9*h, base_x + 7*h*Y', Y' )
  // The SinCos Model's implicit formulation is:
  // f(t,x,xdot) (0) = xdot(0) - x(1)
  // f(t,x,xdot) (1) = xdot(1) + x(0)
  // F_1_exact (0) = Y'(0) - (base_x(1) + 7*h*Y'(1))
  // F_1_exact (1) = Y'(1) + (base_x(0) + 7*h*Y'(0))
  // base_x = (2.0,2.5)
  // Y' = (5.0,5.5)
  // h = 6.0
  // F_1_exact (0) = 5   - (2.5 + 7*6*5.5) = -228.5
  // F_1_exact (1) = 5.5 + (2   + 7*6*5  ) =  217.5
  RCP<Thyra::ProductVectorBase<double> > f_out_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(f_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(f_out_b), true );
  {
    Thyra::ConstDetachedVectorView<double> f_out_view( *f_out );
    TEST_FLOATING_EQUALITY( f_out_view[0], -228.5, tol );
    TEST_FLOATING_EQUALITY( f_out_view[1],  217.5, tol );
  }

  // What should W_op_out be?
  // W_op_out_{i,k} = \frac{\partial F_i}{\partial Y_k^'}
  //                =                \frac{\partial F}{\partial x   }(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*h*a_{i,k}
  //                  + \delta_{i==k}\frac{\partial F}{\partial xdot}(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*1
  // For numStages = 1:
  // W_op_out_{1,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+c*h, base_x + h*a*Y', Y')*h*a 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c*h, base_x + h*a*Y', Y')*1
  // For A=a=7, c=9:
  // W_op_out_{1,1} =   \frac{\partial F}{\partial x   }(t_{n-1}+9*h, base_x + 7*h*Y', Y')*h*7
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+9*h, base_x + 7*h*Y', Y')
  // The SinCos Model's implicit formulation's W_op is \alpha*\frac{\partial F}{\partial xdot} + \beta*\frac{\partial F}{\partial x}:
  // W_op (0,0) = alpha
  // W_op (0,1) = -beta
  // W_op (1,0) =  beta
  // W_op (1,1) = alpha
  // In the IRK method:  alpha = 1.0, and beta = a*h = 7*6.0
  // W_op_out_exact (0,0) = alpha   =   1.0
  // W_op_out_exact (0,1) = -beta*h = -42.0
  // W_op_out_exact (1,0) =  beta*h =  42.0
  // W_op_out_exact (1,1) = alpha   =   1.0
  RCP<Thyra::BlockedLinearOpBase<double> > W_op_out_b = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b), true);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_mv_view( *W_op_out_mv );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,1), -42.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,0),  42.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,1),   1.0, tol );
  }

}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, setCurrentStage ) {
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = 
    getDiagonalImplicitRKModelEvaluator(4);
  TEST_THROW(dirkME->setCurrentStage(-1), std::logic_error);
  for (int stage=0 ; stage<4 ; ++stage) {
    dirkME->setCurrentStage(stage);
  }
  TEST_THROW(dirkME->setCurrentStage(4), std::logic_error);
}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, evalModelTwoStage ) {
  double tol=1.0e-10; 
  
  // Create a Model Evaluator
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
//  RCP<ParameterList> paramList = Teuchos::parameterList();
//  sublist(paramList,Stratimikos_name);
//  sublist(paramList,DiagonalTransientModel_name);
//  RCP<Thyra::ModelEvaluator<double> > model = getDiagonalModel<double>(paramList);
  TEST_EQUALITY_CONST( is_null(model), false );
  
  // create an valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x), false );
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 2.0;
    base_x_view[1] = 2.5;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  TEST_EQUALITY_CONST( is_null(base_x_dot), false );
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 3.0;
    base_x_dot_view[1] = 3.5;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 4.0;
  basePoint.set_t(base_t);

  // create valid W_factory
//  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(paramList);
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  TEST_EQUALITY_CONST( is_null(irk_W_factory), false );
  
  int numStages = 2;
  // Create an IRK Butcher Tableau
  SerialDenseMatrix<int,double> A(numStages, numStages);
  SerialDenseVector<int,double> b(numStages);
  SerialDenseVector<int,double> c(numStages);
  // Backward Euler
  A(0,0) = 7.1;
  A(0,1) = 0.0; // DIRK tableau 
  A(1,0) = 7.3;
  A(1,1) = 7.4;
  b(0) = 8.0;
  b(1) = 8.1;
  c(0) = 9.0;
  c(1) = 9.1;
  RCP<RKButcherTableauBase<double> > irkButcherTableau = rKButcherTableau<double>(A,b,c,2);

  // Create an IRKModelEvaluator
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = 
    diagonalImplicitRKModelEvaluator<double>( model, basePoint, irk_W_factory, irkButcherTableau );

  // Create IRKModelEvaluator InArgs and OutArgs
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = dirkME->createInArgs();
  RCP<VectorBase<double> > x_in = createMember(*(dirkME->get_x_space()));
  TEST_EQUALITY_CONST( is_null(x_in), false );
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = dirkME->createOutArgs();
  RCP<VectorBase<double> > f_out = createMember(*(dirkME->get_f_space()));
  TEST_EQUALITY_CONST( is_null(f_out), false );
  RCP<Thyra::LinearOpBase<double> > W_op_out = dirkME->create_W_op();
  TEST_EQUALITY_CONST( is_null(W_op_out), false );
  outArgs.set_f(f_out);
  outArgs.set_W_op(W_op_out);

  // Fill x_in with valid data.
  // Note:  These are product vectors in the numStages>1 case
  RCP<Thyra::ProductVectorBase<double> > x_in_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(x_in,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(x_in_b), true);
  // stage = 0
  {
    Thyra::DetachedVectorView<double> x_in_view( *x_in );
    x_in_view[0] = 5.1;
    x_in_view[1] = 5.2;
  }

  inArgs.set_x(x_in);
  
  // We need to call setTimeStepPoint before evalModel
  TEST_THROW( dirkME->evalModel(inArgs,outArgs), std::logic_error );
      
  double delta_t = 6.0;
  dirkME->setTimeStepPoint( base_x, base_t, delta_t );
  // We need to call setCurrentStage before evalModel
  TEST_THROW(dirkME->evalModel(inArgs,outArgs), std::logic_error);
  dirkME->setCurrentStage(0);
  dirkME->evalModel(inArgs,outArgs);

  // Verify contents of f_out and W_op_out.
  // 
  // What should f_out be?
  // F(t,x,xdot)=0, this is from Petzold/Brennan/Campbell.
  // h = delta_t
  // M = numStages
  // F_i = F( t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^' )
  // For numSages = 2:
  // F_1 = F( t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'+a_{1,2}*Y_2^'), Y_1^' )
  // F_2 = F( t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^' )
  // for our specific A and c:
  // F_1 = F( t_{n-1}+9.0*h, base_x + h*(7.1*Y_1^'          ), Y_1^' )
  // F_2 = F( t_{n-1}+9.1*h, base_x + h*(7.3*Y_1^'+7.4*Y_2^'), Y_2^' )
  // The SinCos Model's implicit formulation is:
  // f(t,x,xdot) (0) = xdot(0) - x(1)
  // f(t,x,xdot) (1) = xdot(1) + x(0)
  // F_1_exact (0) = Y_1^'(0) - (base_x(1) + h*(7.1*Y_1^'(1)             ))
  // F_1_exact (1) = Y_1^'(1) + (base_x(0) + h*(7.1*Y_1^'(0)             ))
  // F_2_exact (0) = Y_2^'(0) - (base_x(1) + h*(7.3*Y_1^'(1)+7.4*Y_2^'(1)))
  // F_2_exact (1) = Y_2^'(1) + (base_x(0) + h*(7.3*Y_1^'(0)+7.4*Y_2^'(0)))
  // base_x = (2.0,2.5)
  // Y' = [ (5.1,5.2), (5.6,5.7) ]
  // h = 6.0
  // F_1_exact (0) = 5.1 - ( 2.5 + 6*(7.1*5.2           ) ) = -218.92
  // F_1_exact (1) = 5.2 + ( 2.0 + 6*(7.1*5.1           ) ) =  224.46
  // Here we'll setCurrentStage(1) and setStageSolution(0,[5.3,5.4])
  // F_2_exact (0) = 5.6 - ( 2.5 + 6*(7.3*5.4 + 7.4*5.7 ) ) = -486.50
  // F_2_exact (1) = 5.7 + ( 2.0 + 6*(7.3*5.3 + 7.4*5.6 ) ) =  488.48
  //
  RCP<Thyra::ProductVectorBase<double> > f_out_b = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(f_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(f_out_b), true);
  // f_out for stage = 0
  {
    Thyra::ConstDetachedVectorView<double> f_out_view( *f_out );
    TEST_FLOATING_EQUALITY( f_out_view[0], -218.92, tol );
    TEST_FLOATING_EQUALITY( f_out_view[1],  224.46, tol );
  }

  // What should W_op_out be?
  // W_op_out_{i,k} = \frac{\partial F_i}{\partial Y_k^'}
  //                =                \frac{\partial F}{\partial x   }(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*h*a_{i,k}
  //                  + \delta_{i==k}\frac{\partial F}{\partial xdot}(t_{n-1}+c_i*h, base_x + h*sum_{j=1..M}(a_{i,j}*Y_i^'), Y_i^')*1
  // For numStages = 2:
  // DIRK stage = 0:
  // W_op_out       =   \frac{\partial F}{\partial x   }(t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'              ), Y_1^')*h*a_{1,1} 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c_1*h, base_x + h*(a_{1,1}*Y_1^'              ), Y_1^')*1
  // DIRK stage = 1:
  // W_op_out       =   \frac{\partial F}{\partial x   }(t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^')*h*a_{2,2} 
  //                  + \frac{\partial F}{\partial xdot}(t_{n-1}+c_2*h, base_x + h*(a_{2,1}*Y_1^'+a_{2,2}*Y_2^'), Y_2^')*1
  //
  // The SinCos Model's implicit formulation's W_op is \alpha*\frac{\partial F}{\partial xdot} + \beta*\frac{\partial F}{\partial x}:
  // W_op (0,0) = alpha
  // W_op (0,1) = -beta
  // W_op (1,0) =  beta
  // W_op (1,1) = alpha
  // In the DIRK method:  alpha and beta vary by the stage:
  // stage = 0:
  //   alpha = 1  beta = h*a_{1,1} = 6*7.1 = 42.6
  // stage = 1:
  //   alpha = 1  beta = h*a_{2,2} = 6*7.4 = 44.4
  //
  // Stage = 0:
  // W_op_out_exact (0,0) = alpha =   1.0
  // W_op_out_exact (0,1) = -beta = -42.6
  // W_op_out_exact (1,0) =  beta =  42.6
  // W_op_out_exact (1,1) = alpha =   1.0
  //
  // Stage = 1:
  // W_op_out_exact (0,0) = alpha =   1.0
  // W_op_out_exact (0,1) = -beta = -44.4
  // W_op_out_exact (1,0) =  beta =  44.4
  // W_op_out_exact (1,1) = alpha =   1.0
  RCP<Thyra::BlockedLinearOpBase<double> > W_op_out_b = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_b), true);
  RCP<const Thyra::MultiVectorBase<double> > W_op_out_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_mv), false);
  // W_op_out for stage = 0
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_mv_view( *W_op_out_mv );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,1), -42.6, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,0),  42.6, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,1),   1.0, tol );
  }

  // Now we move to the next stage by setting a converged stage solution,
  // updating the currentStage, and calling evalModel again.

  // x_in for "converged" stage = 0
  {
    Thyra::DetachedVectorView<double> x_in_view( *x_in );
    x_in_view[0] = 5.3;
    x_in_view[1] = 5.4;
  }
  dirkME->setStageSolution(0,*x_in);

  // x_in for stage = 1
  {
    Thyra::DetachedVectorView<double> x_in_view( *x_in );
    x_in_view[0] = 5.6;
    x_in_view[1] = 5.7;
  }
  dirkME->setCurrentStage(1);
  dirkME->evalModel(inArgs,outArgs);

  // f_out for stage = 1
  {
    Thyra::ConstDetachedVectorView<double> f_out_view( *f_out );
    TEST_FLOATING_EQUALITY( f_out_view[0], -486.50, tol );
    TEST_FLOATING_EQUALITY( f_out_view[1], 488.48, tol );
  }

  // W_op_out for stage = 1
  //RCP<const Thyra::MultiVectorBase<double> > W_op_out_mv = Teuchos::rcp_dynamic_cast<const Thyra::MultiVectorBase<double> >(W_op_out,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_out_mv), false);
  {
    Thyra::ConstDetachedMultiVectorView<double> W_op_out_mv_view( *W_op_out_mv );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,0),   1.0, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(0,1), -44.4, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,0),  44.4, tol );
    TEST_FLOATING_EQUALITY( W_op_out_mv_view(1,1),   1.0, tol );
  }


}


TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, createInArgs ) {
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = getDiagonalImplicitRKModelEvaluator(1);
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<double> inArgs = dirkME->createInArgs();
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_t), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x), true );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_alpha), false );
  TEST_EQUALITY_CONST( inArgs.supports(MEB::IN_ARG_beta), false );
}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, createOutArgs ) {
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = getDiagonalImplicitRKModelEvaluator(1);
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<double> outArgs = dirkME->createOutArgs();
  TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_f), true );
  TEST_EQUALITY_CONST( outArgs.supports(MEB::OUT_ARG_W_op), true );
}


TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, spaces ) {
  for (int numStages=1 ; numStages <= 5 ; ++numStages) {
    RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = getDiagonalImplicitRKModelEvaluator(numStages);
    RCP<const VectorSpaceBase<double> > x_space = dirkME->get_x_space();
    RCP<const Thyra::ProductVectorSpaceBase<double> > x_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(x_space,false);
    TEST_EQUALITY_CONST( Teuchos::is_null(x_pspace), true );

    RCP<const VectorSpaceBase<double> > f_space = dirkME->get_f_space();
    RCP<const Thyra::ProductVectorSpaceBase<double> > f_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(f_space,false);
    TEST_EQUALITY_CONST( Teuchos::is_null(f_pspace), true );
  }
}

TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, nominalValues ) {
  int numStages = 2;
  RCP<DiagonalImplicitRKModelEvaluator<double> > dirkME = getDiagonalImplicitRKModelEvaluator(numStages);
  typedef ModelEvaluatorBase MEB;
  MEB::InArgs<double> nominalValues = dirkME->getNominalValues();
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_t), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x), true );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_x_dot), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_alpha), false );
  TEST_EQUALITY_CONST( nominalValues.supports(MEB::IN_ARG_beta), false );
//  RCP<const Thyra::VectorBase<double> > x_nom = nominalValues->get_x();
//  RCP<const Thyra::VectorSpaceBase<double> > x_nom_space = x_nom->space();
//  RCP<const Thyra::ProductVectorSpaceBase<double> > x_nom_pspace = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(x_nom_space,false);
//  TEST_EQUALITY_CONST( Teuchos::is_null(x_nom_pspace), false);
//  TEST_EQUALITY( x_nom_pspace->numBlocks(), numStages );

//  RCP<const Thyra::ProductVectorBase<double> > x_nom_b = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorBase<double> >(x_nom,false);
//  TEST_EQUALITY_CONST( Teuchos::is_null(x_nom_b), false );

  // What is a valid nominal value for this ME?

}


TEUCHOS_UNIT_TEST( Rythmos_DIRKModelEvaluator, create_W_op ) {
  int numStages = 3;
  RCP<DiagonalImplicitRKModelEvaluator<double> > irkME = getDiagonalImplicitRKModelEvaluator(numStages);
  RCP<Thyra::LinearOpBase<double> > W_op = irkME->create_W_op();
  RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > W_op_b = Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<double> >(W_op,false);
  TEST_EQUALITY_CONST( Teuchos::is_null(W_op_b), true );
}

TEUCHOS_UNIT_TEST( Rythmos_ImplicitRKStepper, DIRKtakeStep ) {
  // Create the model
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
  // create a valid base point
  Thyra::ModelEvaluatorBase::InArgs<double> basePoint = model->createInArgs();
  RCP<VectorBase<double> > base_x = Thyra::createMember(model->get_x_space());
  {
    Thyra::DetachedVectorView<double> base_x_view( *base_x );
    base_x_view[0] = 0.2;
    base_x_view[1] = 0.3;
  }
  basePoint.set_x(base_x);
  RCP<VectorBase<double> > base_x_dot = Thyra::createMember(model->get_x_space());
  {
    Thyra::DetachedVectorView<double> base_x_dot_view( *base_x_dot );
    base_x_dot_view[0] = 0.6;
    base_x_dot_view[1] = 0.7;
  }
  basePoint.set_x_dot(base_x_dot);
  double base_t = 0.5;
  basePoint.set_t(base_t);
  // Create the nonlinear solver
  RCP<Rythmos::TimeStepNonlinearSolver<double> >
    nonlinearSolver = Rythmos::timeStepNonlinearSolver<double>();
  // Create the IRK W factory
  RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = 
    Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();

  // Create the IRK Butcher Tableau
  RCP<RKButcherTableauBase<double> > irkbt = createRKBT<double>("Backward Euler");
  // Create the IRK Stepper
  RCP<ImplicitRKStepper<double> > irkStepper = 
    implicitRKStepper<double>( model, nonlinearSolver, irk_W_factory, irkbt );
  irkStepper->setInitialCondition(basePoint);
  double dt = 0.75;
  double stepTaken = irkStepper->takeStep(dt, STEP_TYPE_FIXED);
  TEST_EQUALITY( stepTaken, dt );
  // What should the answer be?
  // By writing down Backward Euler applied to the SinCosModel and solving the
  // linear system for the stage-derivative and then substituting this into the
  // assembled IRK solution, we should get:
  // X_bar = [X,Y]^T  this is the stage derivative we're solving for in the nonlinear solver
  // Y = (-1/(1+h^2))*(x_{n-1}+h*y_{n-1})  
  // X = y_{n-1}+h*Y
  // x_{n} = x_{n-1} + h*X
  // y_{n} = y_{n-1} + h*Y
  //
  // x_{n} = x_{n-1} + h*(y_{n-1}+h*(-1/(1+h^2))*(x_{n-1}+h*y_{n-1}))
  // y_{n} = y_{n-1} + h*((-1/(1+h^2))*(x_{n-1}+h*y_{n-1}))
  //
  // x_{n} = 0.2 + 0.75*(0.3+0.75*(-1/(1+0.75^2))*(0.2+0.75*0.3)) = 0.272
  // y_{n} = 0.3 + 0.75*((-1/(1+0.75^2))*(0.2+0.75*0.3)) = 0.096
  //
  //
  // Alternative approach.  Try solving the SinCosModel with Backward Euler
  // directly, without looking at the RK equations.
  // In this approach, I get:
  // x_{n} = x_{n-1}+h*(1/(1+h^2))*(y_{n-1}-h*x_{n-1})
  // y_{n} = (1/(1+h^2))*(y_{n-1}-h*x_{n-1})
  //
  // x_{n} = 0.2+0.75*(1/(1+0.75^2))*(0.3-0.75*0.2) = 0.272
  // y_{n} = (1/(1+0.75^2))*(0.3-0.75*0.2) = 0.096
  //
  //
  double tol = 1.0e-10;
  RCP<const VectorBase<double> > sol = irkStepper->getStepStatus().solution;
  {
    Thyra::ConstDetachedVectorView<double> sol_view( *sol );
    TEST_FLOATING_EQUALITY( sol_view[0], 0.272, tol );
    TEST_FLOATING_EQUALITY( sol_view[1], 0.096, tol );
  }
}


} // namespace Rythmos

