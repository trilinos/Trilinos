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

#include "Rythmos_AdjointModelEvaluator.hpp"
#include "Rythmos_UnitTestHelpers.hpp"
#include "Rythmos_Types.hpp"

#include "Thyra_DetachedVectorView.hpp"

#include "../SinCos/SinCosModel.hpp"

namespace Rythmos  {

TEUCHOS_UNIT_TEST( Rythmos_AdjointModelEvaluator, create ) {
  RCP<AdjointModelEvaluator<double> > ame = rcp(new AdjointModelEvaluator<double>());
  TEST_ASSERT(!is_null(ame));
}

TEUCHOS_UNIT_TEST( Rythmos_AdjointModelEvaluator, evalLinearModel ) {
  RCP<SinCosModel> fwdModel = sinCosModel(true);
  TimeRange<double> fwdTimeRange(0.0,1.0);
  RCP<AdjointModelEvaluator<double> > ame = adjointModelEvaluator<double>(fwdModel,fwdTimeRange);

  // in args:
  Thyra::ModelEvaluatorBase::InArgs<double> inArgs = ame->createInArgs();
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_x_space();
    RCP<VectorBase<double> > x = Thyra::createMember(space);
    RCP<VectorBase<double> > xdot = Thyra::createMember(space);
    {
      Thyra::DetachedVectorView<double> x_view( *x );
      x_view[0] = 3.0;
      x_view[1] = 4.0;
      Thyra::DetachedVectorView<double> xdot_view( *xdot );
      xdot_view[0] = 5.0;
      xdot_view[1] = 6.0;
    }
    inArgs.set_x(x);
    inArgs.set_x_dot(xdot);
  }

  // out args:
  Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = ame->createOutArgs();
  RCP<VectorBase<double> > f_out;
  {
    RCP<const Thyra::VectorSpaceBase<double> > space = ame->get_f_space();
    f_out = Thyra::createMember(space);
    V_S(Teuchos::outArg(*f_out),0.0);
    outArgs.set_f(f_out);
  }

  /*
  // eval model:
  ame->evalModel(inArgs,outArgs);

  // verify output:
  {
    Thyra::ConstDetachedVectorView<double> f_view( *f_out );
    TEST_EQUALITY_CONST( f_view[0], 1.0 );
    TEST_EQUALITY_CONST( f_view[1], 1.0 );
  }
  */
}

} // namespace Rythmos



