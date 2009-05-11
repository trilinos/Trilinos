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

//#include "Rythmos_GlobalErrorEstimator.hpp"
#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "Rythmos_PointwiseInterpolationBufferAppender.hpp"
#include "Rythmos_CubicSplineInterpolator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_TrailingInterpolationBufferAcceptingIntegratorBase.hpp"
#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearNonlinearSolver.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

//TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, create ) {
//  RCP<GlobalErrorEstimator<double> > gee = globalErrorEstimator<double>();
//  TEST_ASSERT( !is_null(gee) );
//}

TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, SinCos ) {
  typedef Teuchos::ScalarTraits<double> ST;
  // Forward Solve, storing data in linear interpolation buffer
  int storageLimit = 100;
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  {
    RCP<ParameterList> ibPL = Teuchos::parameterList();
    ibPL->sublist("Integrator Settings").sublist("Integrator Selection").set("Integrator Type","Default Integrator");
    ibPL->sublist("Integrator Settings").set("Final Time",1.0);
    ibPL->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",0.1);

    ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit RK");
    ibPL->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Backward Euler");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").set("Interpolation Buffer Type","Interpolation Buffer");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").sublist("Interpolation Buffer").set("StorageLimit",storageLimit);
    ibPL->sublist("Interpolation Buffer Settings").sublist("Interpolator Selection").set("Interpolator Type","Linear Interpolator");
    ib->setParameterList(ibPL);
  }
  RCP<SinCosModel> fwdModel = sinCosModel(true); // implicit formulation
  Thyra::ModelEvaluatorBase::InArgs<double> fwdIC = fwdModel->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > fwdNLSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > fwdIntegrator = ib->create(fwdModel,fwdIC,fwdNLSolver);
  {
    Array<double> time_vec;
    time_vec.push_back(1.0);
    fwdIntegrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  }
  // Copy InterpolationBuffer data into Cubic Spline interpolation buffer for use in Adjoint Solve
  TimeRange<double> fwdTimeRange; 
  RCP<InterpolationBufferBase<double> > fwdCubicSplineInterpBuffer;
  {
    RCP<PointwiseInterpolationBufferAppender<double> > piba = pointwiseInterpolationBufferAppender<double>();
    RCP<InterpolationBuffer<double> > sinkInterpBuffer = interpolationBuffer<double>();
    sinkInterpBuffer->setStorage(storageLimit);
    RCP<CubicSplineInterpolator<double> > csi = cubicSplineInterpolator<double>();
    sinkInterpBuffer->setInterpolator(csi);
    RCP<const InterpolationBufferBase<double> > sourceInterpBuffer;
    {
      RCP<TrailingInterpolationBufferAcceptingIntegratorBase<double> > tibaib = 
        Teuchos::rcp_dynamic_cast<TrailingInterpolationBufferAcceptingIntegratorBase<double> >(fwdIntegrator,true);
      sourceInterpBuffer = tibaib->getTrailingInterpolationBuffer();
    }
    fwdTimeRange = sourceInterpBuffer->getTimeRange();
    piba->append(*sourceInterpBuffer, fwdTimeRange, Teuchos::outArg(*sinkInterpBuffer));
    fwdCubicSplineInterpBuffer = sinkInterpBuffer;

    TimeRange<double> sourceRange = sourceInterpBuffer->getTimeRange();
    TimeRange<double> sinkRange = sinkInterpBuffer->getTimeRange();
    TEST_EQUALITY( sourceRange.lower(), sinkRange.lower() );
    TEST_EQUALITY( sourceRange.upper(), sinkRange.upper() );
  }
  // Adjoint Solve, reading forward solve data from Cubic Spline interpolation buffer
  {
    RCP<ParameterList> ibPL = Teuchos::parameterList();
    ibPL->sublist("Integrator Settings").sublist("Integrator Selection").set("Integrator Type","Default Integrator");
    ibPL->sublist("Integrator Settings").set("Final Time",1.0);
    ibPL->sublist("Integration Control Strategy Selection").set("Integration Control Strategy Type","Simple Integration Control Strategy");
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Take Variable Steps",false);
    ibPL->sublist("Integration Control Strategy Selection").sublist("Simple Integration Control Strategy").set("Fixed dt",0.1);

    ibPL->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit RK");
    ibPL->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Implicit 1 Stage 2nd order Gauss");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").set("Interpolation Buffer Type","Interpolation Buffer");
    ibPL->sublist("Interpolation Buffer Settings").sublist("Trailing Interpolation Buffer Selection").sublist("Interpolation Buffer").set("StorageLimit",storageLimit);
    ibPL->sublist("Interpolation Buffer Settings").sublist("Interpolator Selection").set("Interpolator Type","Linear Interpolator");
    ib->setParameterList(ibPL);
  }
  RCP<Thyra::ModelEvaluator<double> > adjModel;
  {
    RCP<Rythmos::AdjointModelEvaluator<double> > model = 
      Rythmos::adjointModelEvaluator<double>(
          fwdModel, fwdTimeRange
          );
    model->setFwdStateSolutionBuffer(fwdCubicSplineInterpBuffer);
    adjModel = model;
  }
  Thyra::ModelEvaluatorBase::InArgs<double> adjIC = adjModel->getNominalValues();
  {
    // Initial conditions for adjoint:
    const RCP<const Thyra::VectorSpaceBase<double> >
      f_space = fwdModel->get_f_space();
    const RCP<Thyra::VectorBase<double> > x_ic = createMember(f_space);
    V_S( Teuchos::outArg(*x_ic), ST::one() );
    const RCP<Thyra::VectorBase<double> > xdot_ic = createMember(f_space);
    V_S( Teuchos::outArg(*xdot_ic), ST::zero() );
    adjIC.set_x(x_ic);
    adjIC.set_x_dot(xdot_ic);
  }
  RCP<Thyra::LinearNonlinearSolver<double> > adjNLSolver = Thyra::linearNonlinearSolver<double>();
  RCP<IntegratorBase<double> > adjIntegrator = ib->create(adjModel,adjIC,adjNLSolver);
  RCP<const VectorBase<double> > phi_final;
  /*
  {
    Array<double> time_vec;
    time_vec.push_back(1.0);
    Array<RCP<const Thyra::VectorBase<double> > > phi_final_array;
    adjIntegrator->getFwdPoints(time_vec,&phi_final_array,NULL,NULL);
    phi_final = phi_final_array[0];
  }
  {
    Thyra::ConstDetachedVectorView<double> phi_view ( *phi_final );
    TEST_EQUALITY_CONST( phi_view[0], 0.0 );
    TEST_EQUALITY_CONST( phi_view[1], 0.0 );
  }
  */
  // Compute error estimate
  //TEST_ASSERT( false );
}

} // namespace Rythmos 


