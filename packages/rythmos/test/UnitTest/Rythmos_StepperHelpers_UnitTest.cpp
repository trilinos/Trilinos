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

#include "Rythmos_StepperHelpers.hpp"
#include "Rythmos_StepperBuilder.hpp"
#include "../SinCos/SinCosModel.hpp"

namespace Rythmos {

// Functions to test:
// assertValidModel
TEUCHOS_UNIT_TEST( Rythmos_StepperHelpers, assertValidModel ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<StepperBuilder<double> > builder = stepperBuilder<double>();
  RCP<StepperBase<double> > stepper = builder->create("Implicit BDF");
  TEST_EQUALITY_CONST( is_null(stepper), false );
  // implicit Stepper and explicit model, throws
  TEST_THROW(
      assertValidModel( *stepper, *model ),
      std::logic_error
      ); 
  stepper = builder->create("Explicit RK");
  // explicit stepper and explicit model, OK
  TEST_NOTHROW(
      assertValidModel( *stepper, *model )
      );
  model = sinCosModel(true);
  // explicit stepper and implicit model, throws
//  TEST_THROW( 
//      assertValidModel( *stepper, *model ),
//      std::logic_error
//      )
  stepper = builder->create("Implicit RK");
  // implicit stepper and implicit model, OK
  TEST_NOTHROW(
      assertValidModel( *stepper, *model )
      );

}

// setDefaultInitialConditionFromNominalValues
// restart

} // namespace Rythmos 


