// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and computes certain
 * derivatives using finite differences.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultFiniteDifferenceModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief Utility object that computes directional finite differences */
  STANDARD_COMPOSITION_MEMBERS(
    DirectionalFiniteDiffCalculator<Scalar>, direcFiniteDiffCalculator );

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator();

  /** \brief . */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel
    ,const RCP<DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}
 
};


/** \brief Nonmember constructor.
 *
 * \relates DefaultFiniteDifferenceModelEvaluator
 */
template<class Scalar>
RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> >
defaultFiniteDifferenceModelEvaluator()
{
  return Teuchos::rcp(new DefaultFiniteDifferenceModelEvaluator<Scalar>());
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultFiniteDifferenceModelEvaluator
 */
template<class Scalar>
RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> >
defaultFiniteDifferenceModelEvaluator(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
  )
{
  RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> > fdModel =
    defaultFiniteDifferenceModelEvaluator<Scalar>();
  fdModel->initialize(thyraModel, direcFiniteDiffCalculator);
  return fdModel;
}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP
