// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_Time.hpp"


//#define THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_DUMP_ALL


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and allows the client to
 * capture the final point that is returned by a client.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultFinalPointCaptureModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultFinalPointCaptureModelEvaluator();

  /** \brief . */
  DefaultFinalPointCaptureModelEvaluator(
    const Teuchos::RCP<ModelEvaluator<Scalar> >  &thyraModel
    );

  /** \brief . */
  const ModelEvaluatorBase::InArgs<Scalar>& getFinalPoint() const;

  /** \brief . */
  bool finalPointWasSolved() const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase */
  //@{

  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private:

  ModelEvaluatorBase::InArgs<Scalar> finalPoint_;
  bool finalPointWasSolved_;
  
};


// /////////////////////////////////
// Implementations


// Constructors/initializers/accessors/utilities


template<class Scalar>
DefaultFinalPointCaptureModelEvaluator<Scalar>::DefaultFinalPointCaptureModelEvaluator()
  :finalPointWasSolved_(false)
{}


template<class Scalar>
DefaultFinalPointCaptureModelEvaluator<Scalar>::DefaultFinalPointCaptureModelEvaluator(
  const Teuchos::RCP<ModelEvaluator<Scalar> >                     &thyraModel
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  finalPoint_ = thyraModel->createInArgs();
  finalPoint_.setArgs(thyraModel->getNominalValues());
  finalPointWasSolved_ = false;
}


template<class Scalar>
const ModelEvaluatorBase::InArgs<Scalar>&
DefaultFinalPointCaptureModelEvaluator<Scalar>::getFinalPoint() const
{
#ifdef THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_DUMP_ALL
  *Teuchos::VerboseObjectBase::getDefaultOStream()
    << "\nDefaultFinalPointCaptureModelEvaluator<Scalar>::getFinalPoint():"
    << " finalPoint =\n" << Teuchos::describe(finalPoint_,Teuchos::VERB_EXTREME);
#endif  
  return finalPoint_;
}


template<class Scalar>
bool DefaultFinalPointCaptureModelEvaluator<Scalar>::finalPointWasSolved() const
{
  return finalPointWasSolved_;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultFinalPointCaptureModelEvaluator<Scalar>::description() const
{
  const Teuchos::RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultFinalPointCaptureModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Overridden from ModelEvaulator.


template<class Scalar>
void DefaultFinalPointCaptureModelEvaluator<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
  const bool wasSolved
  )
{
  finalPoint_.setArgs(finalPoint);
  finalPointWasSolved_ = wasSolved;
  if(!this->isUnderlyingModelConst())
    this->getNonconstUnderlyingModel()->reportFinalPoint(finalPoint,wasSolved);
#ifdef THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_DUMP_ALL
  *Teuchos::VerboseObjectBase::getDefaultOStream()
    << "\nDefaultFinalPointCaptureModelEvaluator<Scalar>::reportFinalPoint(...):"
    << " finalPoint =\n" << Teuchos::describe(finalPoint_,Teuchos::VERB_EXTREME);
#endif  
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
void DefaultFinalPointCaptureModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultFinalPointCaptureModelEvaluator",inArgs,outArgs
    );

  thyraModel->evalModel(inArgs,outArgs);

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();

}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINAL_POINT_CAPTURE_MODEL_EVALUATOR_HPP
