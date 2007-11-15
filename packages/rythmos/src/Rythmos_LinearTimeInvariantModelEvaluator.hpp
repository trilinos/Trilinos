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


#ifndef RYTHMOS_LINEAR_TIME_INVARIANT_MODEL_EVALUATOR_HPP
#define RYTHMOS_LINEAR_TIME_INVARIANT_MODEL_EVALUATOR_HPP


#include "Thyra_StateFuncModelEvaluatorBase.hpp"


namespace Rythmos {


/** \brief .
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class LinearTimeInvariantModelEvaluator
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  LinearTimeInvariantModelEvaluator();

  //@}

  /** \name Overridden from LinearTimeInvariantModelEvaluatorBase */
  //@{

  /** \brief . */
  void initialize(
    const RCP<const Thyra::LinearOpBase<Scalar> > &M,
    const RCP<const Thyra::LinearOpBase<Scalar> > &K,
    const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_factory
    );

  //@}

  /** \name Public functions overridden from ModelEvaluator */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;
  /** \breif . */
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaluatorDefaultBase */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
    ) const;

  //@}

private:

  RCP<const Thyra::LinearOpBase<Scalar> > M_;
  RCP<const Thyra::LinearOpBase<Scalar> > K_;
  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

};


/** \brief Non-member constructor.
 *
 * \relates LinearTimeInvariantModelEvaluator.
 */
template<class Scalar>
RCP<LinearTimeInvariantModelEvaluator<Scalar> >
linearTimeInvariantModelEvaluator(
  const RCP<const Thyra::LinearOpBase<Scalar> > &M,
  const RCP<const Thyra::LinearOpBase<Scalar> > &K,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_factory
  )
{
  RCP<LinearTimeInvariantModelEvaluator<Scalar> >
    model(new LinearTimeInvariantModelEvaluator<Scalar>());
  model->initialize(M,K,W_factory);
  return model;
}


// ///////////////////////
// Definition


// Constructors/initializers/accessors


template<class Scalar>
LinearTimeInvariantModelEvaluator<Scalar>::LinearTimeInvariantModelEvaluator()
{}


// Overridden from LinearTimeInvariantModelEvaluatorBase


template<class Scalar>
void LinearTimeInvariantModelEvaluator<Scalar>::initialize(
  const RCP<const Thyra::LinearOpBase<Scalar> > &M,
  const RCP<const Thyra::LinearOpBase<Scalar> > &K,
  const RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > &W_factory
  )
{
  TEST_FOR_EXCEPT(is_null(M));
  TEST_FOR_EXCEPT(is_null(K));
  TEST_FOR_EXCEPT(is_null(W_factory));
  THYRA_ASSERT_LINEAR_OP_PLUS_LINEAR_OP_SPACES_NAMES(
    "LinearTimeInvariantModelEvaluator<Scalar>::initialize(...)",
    *M, Thyra::NOTRANS, ("M"), *K, Thyra::NOTRANS, ("K")
    );
  M_ = M;
  K_ = K;
  W_factory_ = W_factory;
}


// Public functions overridden from ModelEvaluator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
LinearTimeInvariantModelEvaluator<Scalar>::get_x_space() const
{
  if (is_null(M_))
    return Teuchos::null;
  return M_->domain();
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
LinearTimeInvariantModelEvaluator<Scalar>::get_f_space() const
{
  if (is_null(M_))
    return Teuchos::null;
  return M_->range();
}


template<class Scalar>
RCP<Thyra::LinearOpBase<Scalar> >
LinearTimeInvariantModelEvaluator<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
LinearTimeInvariantModelEvaluator<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
LinearTimeInvariantModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  TEST_FOR_EXCEPT(true);
  return MEB::InArgs<Scalar>();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
LinearTimeInvariantModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x_dot);
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.setSupports(MEB::IN_ARG_t);
  inArgs.setSupports(MEB::IN_ARG_alpha);
  inArgs.setSupports(MEB::IN_ARG_beta);
  return inArgs;
}


// Private functions overridden from ModelEvaluatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
LinearTimeInvariantModelEvaluator<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.set_W_properties(
    MEB::DerivativeProperties(
      MEB::DERIV_LINEARITY_CONST,
      MEB::DERIV_RANK_FULL,
      true // We support adjoints?
      )
    );
  return outArgs;
}


template<class Scalar>
void LinearTimeInvariantModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs
  ) const
{
  TEST_FOR_EXCEPT(true);
}


} // namespace Rythmos


#endif // RYTHMOS_LINEAR_TIME_INVARIANT_MODEL_EVALUATOR_HPP
