// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_NOMINAL_BOUNDS_OVERRIDE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_NOMINAL_BOUNDS_OVERRIDE_MODEL_EVALUATOR_HPP


#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_Time.hpp"


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and allows the client to
 * overide the state contained in the nominal values and the upper and lower
 * bounds.
 *
 * To only override selected nominal values and bounds, do the following:
 *
 \code
 
   template<class Scalar>
   RCP<DefaultNominalBoundsOverrideModelEvaluator<Scalar> >
   override(
     const RCP<ModelEvaluator<Scalar> >   &thyraModel
     ...
     )
   {

     using Teuchos::rcp;
     typedef Thyra::ModelEvaluatorBase MEB;

     // Get the defaults
     RCP<MEB::InArgs<Scalar> >
       nominalValues = clone(thyraModel->getNominalValues()),
       lowerBounds = clone(thyraModel->getLowerBounds()),
       upperBounds = clone(thyraModel->getUpperBounds());

     // Override selected components ...
     ...

     // Initialize the overridden state
     RCP<DefaultNominalBoundsOverrideModelEvaluator<Scalar> >
       defaultOverridder = rcp(
         new DefaultNominalBoundsOverrideModelEvaluator<Scalar>(thyraModel));
     defaultOverridder->setNominalValues(nominalValues);
     defaultOverridder->setLowerBounds(lowerBounds);
     defaultOverridder->setUpperBounds(upperBounds);
     
     return defaultOverridder;

   }
 
 \endcode
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultNominalBoundsOverrideModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultNominalBoundsOverrideModelEvaluator();

  /** \brief . */
  DefaultNominalBoundsOverrideModelEvaluator(
    const RCP<ModelEvaluator<Scalar> > &thyraModel,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &nominalValues,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &lowerBounds = Teuchos::null,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &upperBounds = Teuchos::null
    );

  /** \brief Initalize.
   *
   * \param thyraModel [in] Model being wrapped.
   *
   * \param nominalValues [in] Completely overrides
   * thyraModel->getNominalValues()
   *
   * \param lowerBounds [in] If non-null, completely overrides
   * thyraModel->getLowerBounds()
   *
   * \param upperBounds [in] If non-null, completely overrides
   * thyraModel->getUpperBounds()
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>thyraModel.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getUnderlyingModel.get() == thyraModel.get()</tt>
   * <li>[nominalValues.get()] <tt>this->getNominalValues()</tt> returns <tt>*nominalValues</tt>
   * <li>[lowerBounds.get()] <tt>this->getLowerBounds()</tt> returns <tt>*lowerBounds</tt>
   * <li>[upperBounds.get()] <tt>this->getUpperBounds()</tt> returns <tt>*upperBounds</tt>
   * </ul>
   */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &nominalValues,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &lowerBounds = Teuchos::null,
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &upperBounds = Teuchos::null
    );
  
  /** \brief Set only nominal values. */
  void setNominalValues(
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
    );
  
  /** \brief Set only lower bounds. */
  void setLowerBounds(
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds
    );
  
  /** \brief Set only upper bounds. */
  void setUpperBounds(
    const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds
    );
  
  // ToDo: Add functions to reset lower and upper bounds when needed!

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;

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

  RCP<const ModelEvaluatorBase::InArgs<Scalar> >  nominalValues_;
  RCP<const ModelEvaluatorBase::InArgs<Scalar> >  lowerBounds_;
  RCP<const ModelEvaluatorBase::InArgs<Scalar> >  upperBounds_;
  
};


// /////////////////////////////////
// Implementations


// Constructors/initializers/accessors/utilities


template<class Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::DefaultNominalBoundsOverrideModelEvaluator()
{}


template<class Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::DefaultNominalBoundsOverrideModelEvaluator(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &nominalValues,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &lowerBounds,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &upperBounds
  )
{
  initialize(thyraModel,nominalValues,lowerBounds,upperBounds);
}


template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &nominalValues,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &lowerBounds,
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> > &upperBounds
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  nominalValues_ = nominalValues;
  lowerBounds_ = lowerBounds;
  upperBounds_ = upperBounds;
}


template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::setNominalValues(
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
  )
{
  nominalValues_ = nominalValues;
}


template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::setLowerBounds(
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds
  )
{
  lowerBounds_ = lowerBounds;
}


template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::setUpperBounds(
  const RCP<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds
  )
{
  upperBounds_ = upperBounds;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultNominalBoundsOverrideModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultNominalBoundsOverrideModelEvaluator{";
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
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getNominalValues() const
{
  if(nominalValues_.get())
    return *nominalValues_;
  return this->getUnderlyingModel()->getNominalValues();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getLowerBounds() const
{
  if(lowerBounds_.get())
    return *lowerBounds_;
  return this->getUnderlyingModel()->getLowerBounds();
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getUpperBounds() const
{
  if(upperBounds_.get())
    return *upperBounds_;
  return this->getUnderlyingModel()->getUpperBounds();
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef ModelEvaluatorBase MEB;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultNominalBoundsOverrideModelEvaluator",inArgs,outArgs
    );

  // First set the inArgs to what was overridden
  MEB::InArgs<Scalar>
    wrappedInArgs = ( !is_null(nominalValues_) ? *nominalValues_ : this->createInArgs() );

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\nwrappedInArgs after assigning to nominalValues =\n" << Teuchos::describe(wrappedInArgs,verbLevel);

  // Reset those not at their nominal values
  wrappedInArgs.setArgs(inArgs);

  // This is a special exception: see evalModel() in Thyra::ME
  // documentation.  If inArgs() supports x_dot (or x_dot_dot) but the
  // evaluate call passes in a null value, then we need to make sure
  // the null value gets passed on instead of the nominal value.
  if (wrappedInArgs.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot)) {
    if (is_null(inArgs.get_x_dot())) {
      wrappedInArgs.set_x_dot(Teuchos::null);
    }
  }
  if (wrappedInArgs.supports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot)) {
    if (is_null(inArgs.get_x_dot_dot())) {
      wrappedInArgs.set_x_dot_dot(Teuchos::null);
    }
  }

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\nwrappedInArgs after setting input values =\n" << Teuchos::describe(wrappedInArgs,verbLevel);

  thyraModel->evalModel(wrappedInArgs,outArgs);

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


} // namespace Thyra


#endif // THYRA_DEFAULT_NOMINAL_BOUNDS_OVERRIDE_MODEL_EVALUATOR_HPP
