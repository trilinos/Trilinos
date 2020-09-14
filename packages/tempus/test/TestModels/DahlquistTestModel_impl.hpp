// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP
#define TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
//#include "Thyra_DetachedMultiVectorView.hpp"
//#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
//#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
//#include "Thyra_DefaultLinearOpSource.hpp"
//#include "Thyra_VectorStdOps.hpp"
//#include "Thyra_MultiVectorStdOps.hpp"
//#include "Thyra_DefaultMultiVectorProductVector.hpp"

//#include <iostream>


namespace Tempus_Test {

template<class Scalar>
DahlquistTestModel<Scalar>::DahlquistTestModel()
{ constructDahlquistTestModel(-1.0, false); }


template<class Scalar>
DahlquistTestModel<Scalar>::
DahlquistTestModel(Scalar lambda, bool includeXDot)
{ constructDahlquistTestModel(lambda, includeXDot); }


template<class Scalar>
void
DahlquistTestModel<Scalar>::
constructDahlquistTestModel(Scalar lambda, bool includeXDot)
{
  lambda_        = lambda;
  includeXDot_   = includeXDot;
  isInitialized_ = false;
  xIC_           = Scalar(   1.0);
  xDotIC_        = Scalar(lambda);
  int dim = 1;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim);

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;
  {
    // Set up prototypical InArgs
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports( MEB::IN_ARG_t );
    inArgs.setSupports( MEB::IN_ARG_x );
    inArgs.setSupports( MEB::IN_ARG_x_dot );
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports( MEB::OUT_ARG_f );
    outArgs_ = outArgs;
  }

  // Set up nominal values (initial conditions)
  nominalValues_ = inArgs_;
  {
    nominalValues_.set_t(Scalar(0.0));
    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_ic_view( *x_ic );
      x_ic_view[0] = xIC_;
    }
    nominalValues_.set_x(x_ic);
  }

  if (includeXDot_) {
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view( *x_dot_ic );
      x_dot_ic_view[0] = xDotIC_;
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::
getExactSolution(double t) const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  double exact_t = t;
  inArgs.set_t(exact_t);

  // Set the exact solution, x.
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x = createMember(x_space_);
  { // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    exact_x_view[0] = exp(lambda_*exact_t);
  }
  inArgs.set_x(exact_x);

  // Set the exact solution time derivative, xDot.
  if (includeXDot_) {
    Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x_dot = createMember(x_space_);
    { // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
      exact_x_dot_view[0] = lambda_ * exp(lambda_*exact_t);
    }
    inArgs.set_x_dot(exact_x_dot);
  }

  return inArgs;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::
get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::
get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::
getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::
createInArgs() const
{
  return inArgs_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DahlquistTestModel<Scalar>::
createOutArgsImpl() const
{
  return outArgs_;
}


template<class Scalar>
void
DahlquistTestModel<Scalar>::
evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Thyra::VectorBase;

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view( *x_in );
  const RCP<VectorBase<Scalar> > f_out = outArgs.get_f();

  if (inArgs.get_x_dot().is_null()) {

    // Evaluate the Explicit ODE f(x,t) [= 0]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view( *f_out );
      f_out_view[0] = lambda_*x_in_view[0];
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
        "Error -- Dahlquist Test Model requires f_out!\n");
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      "Error -- Dahlquist Test Model only setup for explicit ODE!\n");
  }

}


} // namespace Tempus_Test
#endif // TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP
