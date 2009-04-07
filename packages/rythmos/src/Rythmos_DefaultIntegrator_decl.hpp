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

#ifndef RYTHMOS_DEFAULT_INTEGRATOR_DECL_HPP
#define RYTHMOS_DEFAULT_INTEGRATOR_DECL_HPP


#include "Rythmos_IntegrationControlStrategyAcceptingIntegratorBase.hpp"
#include "Rythmos_InterpolationBufferAppenderAcceptingIntegratorBase.hpp"
#include "Rythmos_TrailingInterpolationBufferAcceptingIntegratorBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_StepControlInfo.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"




namespace Rythmos {


/** \brief A concrete subclass for <tt>IntegratorBase</tt> that allows a good
 * deal of customization.
 */
template<class Scalar> 
class DefaultIntegrator
  : virtual public IntegrationControlStrategyAcceptingIntegratorBase<Scalar>,
    virtual public InterpolationBufferAppenderAcceptingIntegratorBase<Scalar>,
    virtual public TrailingInterpolationBufferAcceptingIntegratorBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:
  
  /** \brief . */
  typedef typename ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \name Constructors, Initializers, Misc */
  //@{
  
  /** \brief . */
  DefaultIntegrator();

  /** \brief . */
  void setIntegrationObserver(
    const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
    );

  /** \name Overridden from InterpolationBufferAppenderAcceptingIntegratorBase */
  //@{

  /** \brief . */
  void setInterpolationBufferAppender(
    const RCP<InterpolationBufferAppenderBase<Scalar> > &interpBufferAppender
    );

  /** \brief . */
  RCP<const InterpolationBufferAppenderBase<Scalar> >
    getInterpolationBufferAppender();

  /** \brief . */
  RCP<InterpolationBufferAppenderBase<Scalar> >
    getNonconstInterpolationBufferAppender();

  /** \brief . */
  RCP<InterpolationBufferAppenderBase<Scalar> >
    unSetInterpolationBufferAppender();

  //@}
  
  /** \name Overridden from IntegrationControlStrategyAcceptingIntegratorBase */
  //@{
  
  /** \brief . */
  void setIntegrationControlStrategy(
    const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
    );

  /** \brief . */
  RCP<IntegrationControlStrategyBase<Scalar> > 
    getNonconstIntegrationControlStrategy();

  /** \brief . */
  RCP<const IntegrationControlStrategyBase<Scalar> > 
    getIntegrationControlStrategy() const;

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from IntegratorBase */
  //@{

  /** \brief . */
  RCP<IntegratorBase<Scalar> > cloneIntegrator() const;
  
  /** \brief . */
  void setStepper(
    const RCP<StepperBase<Scalar> > &stepper,
    const Scalar &finalTime,
    const bool landOnFinalTime = true
    );

  /** \brief . */
  RCP<StepperBase<Scalar> > unSetStepper();

  /** \brief . */
  RCP<const StepperBase<Scalar> > getStepper() const;

  /** \name Overridden from TrailingInterpolationBufferAcceptingIntegratorBase */
  //@{
  
  /** \brief . */
  void setTrailingInterpolationBuffer(
    const RCP<InterpolationBufferBase<Scalar> > &trailingInterpBuffer
    );

  /** \brief . */
  RCP<InterpolationBufferBase<Scalar> >
    getNonconstTrailingInterpolationBuffer();

  /** \brief . */
  RCP<const InterpolationBufferBase<Scalar> >
    getTrailingInterpolationBuffer() const;

  /** \brief . */
  RCP<InterpolationBufferBase<Scalar> >
    unSetTrailingInterpolationBuffer();

  //@}

  /** \brief . */
  void getFwdPoints(
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
  RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;
    
  /** \brief . */
  void addPoints(
    const Array<Scalar>& time_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& x_vec,
    const Array<RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
    );

  /** \brief . */
  void getPoints(
    const Array<Scalar>& time_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* x_vec,
    Array<RCP<const Thyra::VectorBase<Scalar> > >* xdot_vec,
    Array<ScalarMag>* accuracy_vec
    ) const;

  /** \brief . */
  TimeRange<Scalar> getTimeRange() const;

  /** \brief . */
  void getNodes(Array<Scalar>* time_vec) const;

  /** \brief . */
  void removeNodes(Array<Scalar>& time_vec);

  /** \brief . */
  int getOrder() const;

  //@}

private:

  // ////////////////////////
  // Private data members

  RCP<IntegrationControlStrategyBase<Scalar> > integrationControlStrategy_;
  RCP<IntegrationObserverBase<Scalar> > integrationObserver_;

  RCP<InterpolationBufferBase<Scalar> > trailingInterpBuffer_;
  RCP<InterpolationBufferAppenderBase<Scalar> > interpBufferAppender_;
  
  RCP<StepperBase<Scalar> > stepper_;
  TimeRange<Scalar> integrationTimeDomain_;
  bool landOnFinalTime_;

  int maxNumTimeSteps_;

  int currTimeStepIndex_;
  StepControlInfo<Scalar> stepCtrlInfoLast_;

  static const std::string maxNumTimeSteps_name_;
  static const int maxNumTimeSteps_default_;

  // /////////////////////////
  // Private member functions

  void finalizeSetup();

  bool advanceStepperToTime( const Scalar& t );

};


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
defaultIntegrator();


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
defaultIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy,
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  );


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
controlledDefaultIntegrator(
  const RCP<IntegrationControlStrategyBase<Scalar> > &integrationControlStrategy
  );


/** \brief .
 *
 * \relates DefaultIntegrator
 */
template<class Scalar> 
RCP<DefaultIntegrator<Scalar> >
observedDefaultIntegrator(
  const RCP<IntegrationObserverBase<Scalar> > &integrationObserver
  );

// 2007/08/30: rabartl: Above, note that I had to name the nonmember
// constructors taking an single RCP argument different names from each other
// in order to get around the classic ambiguity problem with implicit
// conversions of smart pointers.


} // namespace Rythmos


#endif //RYTHMOS_DEFAULT_INTEGRATOR_DECL_HPP
