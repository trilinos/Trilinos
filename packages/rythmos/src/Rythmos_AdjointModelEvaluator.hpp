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

#ifndef RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP
#define RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP


#include "Rythmos_IntegratorBase.hpp"
#include "Thyra_ModelEvaluator.hpp" // Interface
#include "Thyra_StateFuncModelEvaluatorBase.hpp" // Implementation
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_Assert.hpp"


namespace Rythmos {


/** \brief Adjoint ModelEvaluator.
 *
 * ToDo: Finish Documentation!
 */
template<class Scalar>
class AdjointModelEvaluator
  : virtual public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /** \name Constructors/Intializers/Accessors */
  //@{

  /** \brief . */
  AdjointModelEvaluator();

  /** \brief . */
  void setFwdStateModel(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel );
  
  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  // /////////////////////////
  // Private data members

  RCP<const Thyra::ModelEvaluator<Scalar> > fwdStateModel_;
};


/** \brief Nonmember constructor.
 *
 * \relates AdjointModelEvaluator
 */
template<class Scalar>
RCP<AdjointModelEvaluator<Scalar> >
adjointModelEvaluator(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel )
{
  RCP<AdjointModelEvaluator<Scalar> >
    adjointModel = Teuchos::rcp(new AdjointModelEvaluator<Scalar>);
  adjointModel->setFwdStateModel(fwdStateModel);
  return adjointModel;
}


// /////////////////////////////////
// Implementations


// Constructors/Intializers/Accessors


template<class Scalar>
AdjointModelEvaluator<Scalar>::AdjointModelEvaluator()
{}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::setFwdStateModel(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &fwdStateModel )
{
  // ToDo: Validate the form of fwdStateModel
  TEST_FOR_EXCEPT(is_null(fwdStateModel));
  fwdStateModel_ = fwdStateModel;
}


// Public functions overridden from ModelEvaulator


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointModelEvaluator<Scalar>::get_x_space() const
{
  return fwdStateModel_->get_f_space();
}


template<class Scalar>
RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointModelEvaluator<Scalar>::get_f_space() const
{
  return fwdStateModel_->get_x_space();
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointModelEvaluator<Scalar>::getNominalValues() const
{
  TEST_FOR_EXCEPT(true);
  return Thyra::ModelEvaluatorBase::InArgs<Scalar>();
}


template<class Scalar>
RCP<Thyra::LinearOpWithSolveBase<Scalar> >
AdjointModelEvaluator<Scalar>::create_W() const
{
  //return Thyra::nonconstAdjoint<Scalar>(fwdStateModel_->create_W());
  TEST_FOR_EXCEPT(true);
  return Teuchos::null;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointModelEvaluator<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> stateModelInArgs = fwdStateModel_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports( MEB::IN_ARG_x_dot );
  inArgs.setSupports( MEB::IN_ARG_x );
  inArgs.setSupports( MEB::IN_ARG_t );
  inArgs.setSupports( MEB::IN_ARG_alpha );
  inArgs.setSupports( MEB::IN_ARG_beta );
  return inArgs;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
AdjointModelEvaluator<Scalar>::createOutArgsImpl() const
{

  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgs<Scalar> stateModelOutArgs = fwdStateModel_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;

  outArgs.setModelEvalDescription(this->description());

  outArgs.setSupports(MEB::OUT_ARG_f);

  if (stateModelOutArgs.supports(MEB::OUT_ARG_W) ) {
    outArgs.setSupports(MEB::OUT_ARG_W);
    outArgs.set_W_properties(stateModelOutArgs.get_W_properties());
  }

  return outArgs;

}


template<class Scalar>
void AdjointModelEvaluator<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::rcp_dynamic_cast;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::VerboseObjectTempState<Thyra::ModelEvaluatorBase> VOTSME;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_GEN_BEGIN(
    "AdjointModelEvaluator", inArgs, outArgs, Teuchos::null );

  TEST_FOR_EXCEPT(true);

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Rythmos


#endif // RYTHMOS_ADJOINT_MODEL_EVALUATOR_HPP
