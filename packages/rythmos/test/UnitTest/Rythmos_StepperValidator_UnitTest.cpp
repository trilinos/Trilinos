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
#include "Teuchos_ObjectBuilder.hpp"

#include "Rythmos_StepperValidator.hpp"
#include "Rythmos_IntegratorBuilder.hpp"
#include "Rythmos_UnitTestModels.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ExplicitRK ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Explicit RK");
  pl->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Forward Euler");
  ib->setParameterList(pl);
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, BackwardEuler ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Backward Euler");
  ib->setParameterList(pl);
  sv->setIntegratorBuilder(ib);
  
  sv->validateStepper();
  TEST_ASSERT( true );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ImplicitRK ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  {
    // Stepper Settings
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit RK");
    pl->sublist("Stepper Settings").sublist("Runge Kutta Butcher Tableau Selection").set("Runge Kutta Butcher Tableau Type","Backward Euler");
    ib->setParameterList(pl);
  }
  {
    // W Factory
    RCP<ParameterList> wPL = Teuchos::parameterList();
    sublist(wPL,"Stratimikos");
    RCP<Thyra::LinearOpWithSolveFactoryBase<double> > irk_W_factory = getWFactory<double>(wPL);
    ib->setWFactoryObject(irk_W_factory);
  }
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}

TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ImplicitBDF ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  {
    // Stepper Settings
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Implicit BDF");
    ib->setParameterList(pl);
  }
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}

/*
TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ForwardSensitivityStepper ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<StepperBuilder<double> > sb = ib->getStepperBuilder();
  // Add the Forward Sensitivity Stepper to the Stepper Builder.
  sb->setStepperFactory(
      Teuchos::abstractFactoryStd< StepperBase<double>, ForwardSensitivityStepper<double> >(),
      "Forward Sensitivity"
      );
  {
    // Stepper Settings
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Forward Sensitivity");
    ib->setParameterList(pl);
  }
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}
*/


#ifdef HAVE_THYRA_ME_POLYNOMIAL


TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ExplicitTaylorPolynomial ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Explicit Taylor Polynomial");
  ib->setParameterList(pl);
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}


#endif // HAVE_THYRA_ME_POLYNOMIAL


TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, ForwardEuler ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Forward Euler");
  ib->setParameterList(pl);
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}

#ifdef HAVE_RYTHMOS_EXPERIMENTAL
TEUCHOS_UNIT_TEST( Rythmos_StepperValidator, Theta ) {
  RCP<StepperValidator<double> > sv = stepperValidator<double>();
  RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->sublist("Stepper Settings").sublist("Stepper Selection").set("Stepper Type","Theta");
  ib->setParameterList(pl);
  sv->setIntegratorBuilder(ib);

  sv->validateStepper();
  TEST_ASSERT( true );
}
#endif // HAVE_RYTHMOS_EXPERIMENTAL

} // namespace Rythmos 


