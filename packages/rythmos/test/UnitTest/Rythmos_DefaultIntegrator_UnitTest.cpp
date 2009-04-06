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

#include "Rythmos_ExplicitRKStepper.hpp"

#include "../SinCos/SinCosModel.hpp"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Thyra::VectorBase;

// Test the ERK stepper through the integrator
TEUCHOS_UNIT_TEST( Rythmos_DefaultIntegrator, ExplicitRKStepper ) {
  // Integrator
  RCP<DefaultIntegrator<double> > integrator = defaultIntegrator<double>();

  // Stepper
  double finalTime = 1.0;
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<RKButcherTableauBase<double> > rkbt = createRKBT<double>("Forward Euler");
  TEST_ASSERT( !is_null(rkbt) );
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  // Set initial condition on stepper
  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();
  stepper->setInitialCondition(ic);
  // Set stepper on integrator:
  integrator->setStepper(stepper, finalTime);

  // Trailing Interpolation Buffer
  // Do not use with ExplicitRK yet as the InterpolationBuffer does not work without x_dot
  
  // IntegrationControlStrategy to specify fixed steps for this stepper
  double dt = 0.1;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Take Variable Steps",false);
  pl->set("Fixed dt", dt);
  RCP<SimpleIntegrationControlStrategy<double> > intCont = simpleIntegrationControlStrategy<double>(pl);
  TimeRange<double> tr(0.0,finalTime);
  intCont->resetIntegrationControlStrategy(tr); // ??? Why do I need to do this?
  integrator->setIntegrationControlStrategy(intCont);

  // Ask integrator for points
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  double N = 10;
  for (int i=0 ; i<=N ; ++i) {
    double t = 0.0 + i*finalTime/N;
    time_vec.push_back(t);
  }
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

  // Verify that these points are accurate
  // Since we're using Forward Euler, we can write down the exact numerical solution
  double tol = 1.0e-10;
  double exact_x0 = 0.0; // nominal values on SinCosModel
  double exact_x1 = 1.0;
  for (int i=0 ; i<=N ; ++i) {
    {
      Thyra::ConstDetachedVectorView<double> x_vec_view( *(x_vec[i]) );
      TEST_FLOATING_EQUALITY( exact_x0, x_vec_view[0], tol );
      TEST_FLOATING_EQUALITY( exact_x1, x_vec_view[1], tol );
    }
    double x0 = exact_x0;
    double x1 = exact_x1;
    exact_x0 += dt*x1;
    exact_x1 -= dt*x0;
  }
}

} // namespace Rythmos

