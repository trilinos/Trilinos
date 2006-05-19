// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
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
 * Hint: To only overide some of the nominal values and bounds you can do:
 *
 \code
 
   template<class Scalar>
   Teuchos::RefCountPtr<DefaultNominalBoundsOverrideModelEvaluator<Scalar> >
   override(
     const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >   &thyraModel
     ...
     )
   {
     // Get the defaults
     typedef ModelEvaluatorBase MEB;
     MEB::InArgs<Scalar>  nominalValues = thyraModel->getNominalValues();
     MEB::InArgs<Scalar>  lowerBounds   = thyraModel->getLowerBounds();
     MEB::InArgs<Scalar>  upperBounds   = thyraModel->getUpperBounds();
     // Override selected components ...
     ...
     // Initialize the overridden state
     return Teuchos::rcp(
       new DefaultNominalBoundsOverrideModelEvaluator<Scalar>(
         thyraModel,nominalValues,lowerBounds,upperBounds
         )
       );
   }
 
 \endcode
 *
 * ToDo: Finish documentation!
 *
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
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                     &thyraModel
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds = Teuchos::null
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds = Teuchos::null
    );

  /** \brief Initalize.
   *
   * \param  thyraModel     [in] Model being wrapped.
   * \param  nominalValues  [in] Completely overrides thyraModel->getNominalValues()
   * \param  lowerBounds    [in] If non-null, completely overrides thyraModel->getLowerBounds()
   * \param  upperBounds    [in] If non-null, completely overrides thyraModel->getUpperBounds()
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
    const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                     &thyraModel
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds = Teuchos::null
    ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds = Teuchos::null
    );
  
  /** \brief Set only nominal values. */
  void setNominalValues(
    const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
    );
  
  // ToDo: Add functions to reset lower and upper bounds when needed!

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  nominalValues_;
  Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  lowerBounds_;
  Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  upperBounds_;
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::DefaultNominalBoundsOverrideModelEvaluator()
{}

template<class Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::DefaultNominalBoundsOverrideModelEvaluator(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                     &thyraModel
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds
  )
{
  initialize(thyraModel,nominalValues,lowerBounds,upperBounds);
}

template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::initialize(
  const Teuchos::RefCountPtr<ModelEvaluator<Scalar> >                     &thyraModel
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &lowerBounds
  ,const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &upperBounds
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  nominalValues_ = nominalValues;
  lowerBounds_ = lowerBounds;
  upperBounds_ = upperBounds;
}

template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::setNominalValues(
  const Teuchos::RefCountPtr<const ModelEvaluatorBase::InArgs<Scalar> >  &nominalValues
  )
{
  nominalValues_ = nominalValues;
}

// Overridden from ModelEvaulator.

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getNominalValues() const
{
  if(nominalValues_.get()) return *nominalValues_;
  return this->getUnderlyingModel()->getNominalValues();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getLowerBounds() const
{
  if(lowerBounds_.get()) return *lowerBounds_;
  return this->getUnderlyingModel()->getLowerBounds();
}

template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
DefaultNominalBoundsOverrideModelEvaluator<Scalar>::getUpperBounds() const
{
  if(upperBounds_.get()) return *upperBounds_;
  return this->getUnderlyingModel()->getUpperBounds();
}

template<class Scalar>
void DefaultNominalBoundsOverrideModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;

  Teuchos::Time totalTimer(""), timer("");
  totalTimer.start(true);

  const Teuchos::RefCountPtr<Teuchos::FancyOStream> out       = this->getOStream();
  const Teuchos::EVerbosityLevel                    verbLevel = this->getVerbLevel();
  Teuchos::OSTab tab(out);
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nEntering Thyra::DefaultNominalBoundsOverrideModelEvaluator<Scalar>::evalModel(...) ...\n";

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\ninArgs =\n" << Teuchos::describe(inArgs,verbLevel)
      << "\noutArgs on input =\n" << Teuchos::describe(outArgs,Teuchos::VERB_LOW);

  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();

  typedef Teuchos::VerboseObjectTempState<ModelEvaluatorBase> VOTSME;
  VOTSME thyraModel_outputTempState(thyraModel,out,verbLevel);

  // First set the inArgs to what was overridden
  MEB::InArgs<Scalar> wrappedInArgs = *nominalValues_;

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\nwrappedInArgs after assigning to nominalValues =\n" << Teuchos::describe(wrappedInArgs,verbLevel);

  // Reset those not at their nominal values
  wrappedInArgs.setArgs(inArgs);

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\nwrappedInArgs after setting input values =\n" << Teuchos::describe(wrappedInArgs,verbLevel);

  thyraModel->evalModel(wrappedInArgs,outArgs);

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_EXTREME))
    *out
      << "\noutArgs on output =\n" << Teuchos::describe(outArgs,verbLevel);

  totalTimer.stop();
  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out
      << "\nLeaving Thyra::DefaultNominalBoundsOverrideModelEvaluator<Scalar>::evalModel(...) ...\n";
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultNominalBoundsOverrideModelEvaluator<Scalar>::description() const
{
  const Teuchos::RefCountPtr<const ModelEvaluator<Scalar> >
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

} // namespace Thyra

#endif // THYRA_DEFAULT_NOMINAL_BOUNDS_OVERRIDE_MODEL_EVALUATOR_HPP
