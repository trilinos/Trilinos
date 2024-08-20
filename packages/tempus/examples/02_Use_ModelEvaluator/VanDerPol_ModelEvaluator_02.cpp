//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"

#include <iostream>

#include "VanDerPol_ModelEvaluator_02.hpp"


template<class Scalar>
VanDerPol_ModelEvaluator_02<Scalar>::
VanDerPol_ModelEvaluator_02()
 : dim_(2),
   t0_ic_  (Scalar(0.0)),
   epsilon_(Scalar(1.0e-01)),
   x0_ic_  (Scalar(2.0)),
   x1_ic_  (Scalar(0.0))
{
  using Teuchos::RCP;
  typedef ::Thyra::ModelEvaluatorBase MEB;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = x_space_;

  {
    // Set up prototypical InArgs
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( MEB::IN_ARG_t );
    inArgs.setSupports( MEB::IN_ARG_x );
    inArgs.setSupports( MEB::IN_ARG_x_dot );
    prototypicalInArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( MEB::OUT_ARG_f );
    prototypicalOutArgs_ = outArgs;
  }

  // Set the initial conditions
  nominalValues_ = prototypicalInArgs_;
  nominalValues_.set_t(t0_ic_);
  const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
    x_ic_view[0] = x0_ic_;
    x_ic_view[1] = x1_ic_;
  }
  nominalValues_.set_x(x_ic);

  const RCP<Thyra::VectorBase<Scalar> > xDot_ic = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> xDot_ic_view( *xDot_ic );
    xDot_ic_view[0] = x1_ic_;
    xDot_ic_view[1] = ((1.0-x0_ic_*x0_ic_)*x1_ic_-x0_ic_)/epsilon_;
  }
  nominalValues_.set_x_dot(xDot_ic);

}


template<class Scalar>
void
VanDerPol_ModelEvaluator_02<Scalar>::
evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;

  const RCP<const Thyra::VectorBase<Scalar> > x_in =
    inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view( *x_in );

  const RCP<Thyra::VectorBase<Scalar> > f_out =
    outArgs.get_f().assert_not_null();

  if (inArgs.get_x_dot().is_null()) {
    // Evaluate the Explicit ODE f(x,t) [= xdot]
    Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
    f_out_view[0] = x_in_view[1];
    f_out_view[1] =
      ((1.0-x_in_view[0]*x_in_view[0])*x_in_view[1]-x_in_view[0])/epsilon_;

  } else {
    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const Thyra::VectorBase<Scalar> > x_dot_in;
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
    Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view( *x_dot_in );
    f_out_view[0] = x_dot_in_view[0] - x_in_view[1];
    f_out_view[1] = x_dot_in_view[1]
      - ((1.0-x_in_view[0]*x_in_view[0])*x_in_view[1]-x_in_view[0])/epsilon_;
  }
}

template class VanDerPol_ModelEvaluator_02<double>;
