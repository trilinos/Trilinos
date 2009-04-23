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

namespace Rythmos {

//TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, create ) {
//  RCP<GlobalErrorEstimator<double> > gee = globalErrorEstimator<double>();
//  TEST_ASSERT( !is_null(gee) );
//}

TEUCHOS_UNIT_TEST( Rythmos_GlobalErrorEstimator, SinCos ) {
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
  RCP<SinCosModel> model = sinCosModel(true); // implicit formulation
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  RCP<Thyra::NonlinearSolverBase<double> > nlSolver = timeStepNonlinearSolver<double>();
  RCP<IntegratorBase<double> > fwdIntegrator = ib->create(model,ic,nlSolver);
  {
    Array<double> time_vec;
    time_vec.push_back(1.0);
    fwdIntegrator->getFwdPoints(time_vec,NULL,NULL,NULL);
  }
  // Copy InterpolationBuffer data into Cubic Spline interpolation buffer for use in Adjoint Solve
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
    TimeRange<double> tr = sourceInterpBuffer->getTimeRange();
    piba->append(*sourceInterpBuffer, tr, Teuchos::outArg(*sinkInterpBuffer));
    fwdCubicSplineInterpBuffer = sinkInterpBuffer;

    TimeRange<double> sourceRange = sourceInterpBuffer->getTimeRange();
    TimeRange<double> sinkRange = sinkInterpBuffer->getTimeRange();
    TEST_EQUALITY( sourceRange.lower(), sinkRange.lower() );
    TEST_EQUALITY( sourceRange.upper(), sinkRange.upper() );
  }
  // Adjoint Solve, reading forward solve data from Cubic Spline interpolation buffer
    // Initial conditions for adjoint:
  // Compute error estimate
}

} // namespace Rythmos 


