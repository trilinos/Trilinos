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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "../SinCos/SinCosModel.hpp"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

#include "Thyra_DetachedVectorView.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_ImplicitBDFStepperStepControl.hpp"
#include "Rythmos_MockStepControlStrategyDecorator.hpp"

namespace Rythmos {


using Teuchos::getParametersFromXmlString;
using Thyra::VectorBase;

TEUCHOS_UNIT_TEST( Rythmos_DefaultIntegrator,
                   ImplicitBDFStepperRampingStepControl )
{
  const RCP<SinCosModel> model = sinCosModel(true);

  const RCP<TimeStepNonlinearSolver<double> > nonlinearSolver =
    timeStepNonlinearSolver<double>();

  const RCP<Teuchos::ParameterList> params = 
    Teuchos::parameterList();
  const RCP<ImplicitBDFStepper<double> > stepper =
    implicitBDFStepper<double>(model, nonlinearSolver, params);

//   const RCP<ImplicitBDFStepperStepControl<double> > scs = 
//     Teuchos::rcp(new ImplicitBDFStepperStepControl<double>);

  const RCP<ImplicitBDFStepperRampingStepControl<double> > rscs = 
    Teuchos::rcp(new ImplicitBDFStepperRampingStepControl<double>);

  rscs->setVerbLevel(Teuchos::VERB_HIGH);

  const RCP<MockStepControlStrategyDecorator<double> > mscsd = 
    mockStepControlStrategyDecorator<double>();
  mscsd->initialize(rscs);
  mscsd->addNonlinearSolverFailureOnStep(4,3);

  stepper->setStepControlStrategy(mscsd);

  stepper->setInitialCondition(model->getNominalValues());
  const RCP<DefaultIntegrator<double> > integrator =
    defaultIntegrator<double>();

  const double finalTime = 1.0;
  integrator->setStepper(stepper, finalTime);
  integrator->setVerbLevel(Teuchos::VERB_MEDIUM);
  integrator->setOStream(Teuchos::rcpFromRef(out));
  const RCP<const Thyra::VectorBase<double> > x_final =
    get_fwd_x<double>(*integrator, finalTime);

  // not the best way to test, but will do for now
  TEST_EQUALITY(rscs->numberOfSteps(), 47);
  TEST_EQUALITY(rscs->numberOfFailedSteps(), 3);
  TEST_FLOATING_EQUALITY(rscs->currentStepSize(),
                         Teuchos::as<double>(0.07864277), 1.0e-7);
  TEST_EQUALITY(rscs->currentOrder(), 5);

}


} // namespace Rythmos

