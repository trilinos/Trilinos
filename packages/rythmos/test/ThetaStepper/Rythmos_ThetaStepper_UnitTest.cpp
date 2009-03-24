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
#include "Rythmos_ThetaStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"
#include "../SinCos/SinCosModel.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_ThetaStepper, create) {
  // Model
  RCP<SinCosModel> model = sinCosModel(true);

  Thyra::ModelEvaluatorBase::InArgs<double> ic = model->getNominalValues();

  // Solver
  RCP<TimeStepNonlinearSolver<double> > nlSolver = timeStepNonlinearSolver<double>();

  // Stepper
  RCP<ThetaStepper<double> > stepper = thetaStepper<double>();
  TEST_ASSERT( !is_null(stepper) );
  stepper->setModel(model);
  stepper->setSolver(nlSolver);
  stepper->setInitialCondition(ic);

  double dt = 0.1;
  double dt_taken = 0.0;
  TEST_NOTHROW( dt_taken = stepper->takeStep(dt,STEP_TYPE_FIXED) );
  TEST_EQUALITY_CONST( dt_taken, dt );

}

} // namespace Rythmos

