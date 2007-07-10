//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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

#ifndef Rythmos_SIMPLE_INTEGRATOR_H
#define Rythmos_SIMPLE_INTEGRATOR_H

#include "Rythmos_IntegratorBase.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_as.hpp"


namespace Rythmos {


/** \brief A very simple concrete subclass for <tt>IntegratorBase</tt> that
 * allows just for simple fixed steps or variable steps.
 */
template<class Scalar> 
class SimpleIntegrator : virtual public IntegratorBase<Scalar>
{
public:

  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  SimpleIntegrator();

  //@}

  /** \name Overridden from IntegratorBase */
  //@{
  
  /** \brief . */
  void setInterpolationBuffer(
    const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    );

  /** \brief . */
  void setStepper(
    const RCP<StepperBase<Scalar> > &stepper
    );

  /** \brief . */
  bool getFwdPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    );

  /** \brief . */
  TimeRange<Scalar> getFwdTimeRange() const;

  //@}

  /** \name Overridden from InterpolationBufferBase */
  //@{
    
  /** \brief . */
  bool setPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
    const Array<ScalarMag> & accuracy_vec 
    );

  /** \brief . */
  bool getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  bool setRange(
    const TimeRange<Scalar>& range,
    const InterpolationBufferBase<Scalar>& interpBuffer
    );

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  bool getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  bool removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}

  /** \name Overridden from Teuchos::Describable */
  //@{

  /** \brief . */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<ParameterList> getParameterList();

  /** \brief . */
  RCP<ParameterList> unsetParameterList();

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

private:

  RCP<StepperBase<Scalar> > stepper_;
  RCP<ParameterList> paramList_;
  bool takeVariableSteps_;
  Scalar fixed_dt_;
  Scalar finalTime_;

};


/** \brief .
 *
 * \relates SimpleIntegrator
 */
template<class Scalar> 
RCP<SimpleIntegrator<Scalar> > simpleIntegrator()
{
  return Teuchos::rcp(new SimpleIntegrator<Scalar>());
}


// ////////////////////////////
// Defintions


// Constructors, Initializers, Misc


template<class Scalar> 
SimpleIntegrator<Scalar>::SimpleIntegrator()
  :takeVariableSteps_(true),
   fixed_dt_(ScalarTraits<Scalar>::zero()),
   finalTime_(ScalarTraits<Scalar>::zero())
{}


// Overridden from IntegratorBase


template<class Scalar> 
void SimpleIntegrator<Scalar>::setInterpolationBuffer(
    const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    )
{
  TEST_FOR_EXCEPT(true); // We will never implement this function in this class!
}


template<class Scalar> 
void SimpleIntegrator<Scalar>::setStepper(
  const RCP<StepperBase<Scalar> > &stepper
  )
{
  TEST_FOR_EXCEPT(is_null(stepper));
  // ToDo: Validate state of the stepper
  stepper_ = stepper;
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::getFwdPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  )
{
  TEST_FOR_EXCEPT("ToDo: Implement the forward integration method");
  return false;
}


template<class Scalar> 
TimeRange<Scalar>
SimpleIntegrator<Scalar>::getFwdTimeRange() const
{
  return timeRange(
    stepper_->getTimeRange().lower(),
    finalTime_
    );
}


// Overridden from InterpolationBufferBase


template<class Scalar> 
bool SimpleIntegrator<Scalar>::setPoints(
  const Array<Scalar>& time_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
  const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec,
  const Array<ScalarMag> & accuracy_vec 
  )
{
  TEST_FOR_EXCEPT(true);
  return false;
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::getPoints(
  const Array<Scalar>& time_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
  Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
  Array<ScalarMag>* accuracy_vec
  ) const
{
  return stepper_->getPoints(time_vec,x_vec,xdot_vec,accuracy_vec);
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::setRange(
  const TimeRange<Scalar>& range,
  const InterpolationBufferBase<Scalar>& interpBuffer
  )
{
  return stepper_->setRange(range,interpBuffer);
}


template<class Scalar> 
TimeRange<Scalar> SimpleIntegrator<Scalar>::getTimeRange() const
{
  return stepper_->getTimeRange();
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::getNodes(Array<Scalar>* time_vec) const
{
  return stepper_->getNodes(time_vec);
}


template<class Scalar> 
bool SimpleIntegrator<Scalar>::removeNodes(Array<Scalar>& time_vec)
{
  return stepper_->removeNodes(time_vec);
}


template<class Scalar> 
int SimpleIntegrator<Scalar>::getOrder() const
{
  return stepper_->getOrder();
}


// Overridden from Teuchos::Describable


template<class Scalar> 
void SimpleIntegrator<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  TEST_FOR_EXCEPT(true);
}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void SimpleIntegrator<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(true);
}


template<class Scalar> 
RCP<ParameterList>
SimpleIntegrator<Scalar>::getParameterList()
{
  return paramList_;
}


template<class Scalar> 
RCP<ParameterList>
SimpleIntegrator<Scalar>::unsetParameterList()
{
  RCP<ParameterList> tempParamList = paramList_;
  paramList_ = Teuchos::null;
  return tempParamList;
}


template<class Scalar> 
RCP<const ParameterList>
SimpleIntegrator<Scalar>::getValidParameters() const
{
  return paramList_;
}


} // namespace Rythmos


#endif //Rythmos_SIMPLE_INTEGRATOR_H
