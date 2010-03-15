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

#ifndef RYTHMOS_STEPPER_HELPERS_DECL_HPP
#define RYTHMOS_STEPPER_HELPERS_DECL_HPP


#include "Rythmos_Types.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Rythmos_InterpolatorBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Rythmos {


/** \brief Assert valid transient model for StepperBase.
 *
 */
template<class Scalar>
void assertValidModel(
  const StepperBase<Scalar>& stepper,
  const Thyra::ModelEvaluator<Scalar>& model
  );


/** \brief Set an initial condition on a stepper from a model if the stepper
 * does not already have an initial condition.
 *
 * \returns Returns <tt>true</tt> if the stepper->setInitialCondition(...) was
 * called and returns <tt>false</tt> otherwise.
 */
template<class Scalar>
bool setDefaultInitialConditionFromNominalValues(
  const Thyra::ModelEvaluator<Scalar>& model,
  const Ptr<StepperBase<Scalar> >& stepper
  );

/** \brief Restart a time stepper.
 *
 * This simple helper function just grabs the state out of the
 * <tt>*stepper</tt> object and then resets it on itself as an initial
 * condition.  This is generally used to restart a stepper when passing over a
 * breakpoint where the model is expected to be discontinuous in some way.
 */
template<class Scalar>
void restart( StepperBase<Scalar> *stepper );

template<class Scalar>
void eval_model_explicit(
    const Thyra::ModelEvaluator<Scalar> &model,
    Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const VectorBase<Scalar>& x_in,
    const typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t_in,
    const Ptr<VectorBase<Scalar> >& f_out
    );


#ifdef HAVE_THYRA_ME_POLYNOMIAL


template<class Scalar>
void eval_model_explicit_poly(
    const Thyra::ModelEvaluator<Scalar> &model,
    Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const Teuchos::Polynomial< VectorBase<Scalar> > &x_poly,
    const typename Thyra::ModelEvaluatorBase::InArgs<Scalar>::ScalarMag &t,
    const Ptr<Teuchos::Polynomial<VectorBase<Scalar> > >& f_poly
    );


#endif // HAVE_THYRA_ME_POLYNOMIAL


// This function simply returns the boundary points if they're asked for.  Otherwise it throws.
template<class Scalar>
void defaultGetPoints(
    const Scalar& t_old, // required inArg
    const Ptr<const VectorBase<Scalar> >& x_old, // optional inArg
    const Ptr<const VectorBase<Scalar> >& xdot_old, // optional inArg
    const Scalar& t, // required inArg
    const Ptr<const VectorBase<Scalar> >& x, // optional inArg
    const Ptr<const VectorBase<Scalar> >& xdot, // optional inArg
    const Array<Scalar>& time_vec, // required inArg
    const Ptr<Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > >& x_vec, // optional outArg
    const Ptr<Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > >& xdot_vec, // optional outArg
    const Ptr<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> >& accuracy_vec, // optional outArg
    const Ptr<InterpolatorBase<Scalar> > interpolator // optional inArg (note:  not const)
    );

// This function sets a model on a stepper by creating the appropriate
// ConstNonconstObjectContainer object.
template<class Scalar>
  void setStepperModel(
      const Ptr<StepperBase<Scalar> >& stepper,
      const RCP<const Thyra::ModelEvaluator<Scalar> >& model
      );

template<class Scalar>
  void setStepperModel(
      const Ptr<StepperBase<Scalar> >& stepper,
      const RCP<Thyra::ModelEvaluator<Scalar> >& model
      );

template<class Scalar>
  void setStepperModel(
      const Ptr<StepperBase<Scalar> >& stepper,
      Teuchos::ConstNonconstObjectContainer<Thyra::ModelEvaluator<Scalar> >& 
        model
      );


} // namespace Rythmos


#endif // RYTHMOS_STEPPER_HELPERS_DECL_HPP
