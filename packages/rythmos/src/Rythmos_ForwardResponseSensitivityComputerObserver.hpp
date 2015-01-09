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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_OBSERVER_HPP
#define RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_OBSERVER_HPP


#include "Rythmos_IntegrationObserverBase.hpp"
#include "Rythmos_ForwardResponseSensitivityComputer.hpp"
#include "Rythmos_ResponseAndFwdSensPoint.hpp"
#include "Rythmos_extractStateAndSens.hpp"


namespace Rythmos {


/** \brief Observer class that computes sensitivities at the end of each time
 * step.
 *
 * ToDo: Finish Documentation
 */
template<class Scalar>
class ForwardResponseSensitivityComputerObserver
  : public IntegrationObserverBase<Scalar>
{
public:

  /** \name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  ForwardResponseSensitivityComputerObserver();

  /** \brief . */
  void initialize(
    const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
    const int p_index,
    const int g_index
    );

  /** \brief . */
  const Array<ResponseAndFwdSensPoint<Scalar> >& responseAndFwdSensPoints() const;

  //@}

  /** \name Overridden from IntegrationObserverBase */
  //@{

  /** \brief . */
  virtual RCP<IntegrationObserverBase<Scalar> >
  cloneIntegrationObserver() const;

  /** \brief . */
  virtual void resetIntegrationObserver(
    const TimeRange<Scalar> &integrationTimeDomain
    );

  /** \brief . */
  virtual void observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  //@}

private:

  ForwardResponseSensitivityComputer<Scalar> forwardResponseSensitivityComputer_; 
  Array<ResponseAndFwdSensPoint<Scalar> > responseAndFwdSensPoints_;

  RCP<Thyra::VectorBase<Scalar> > g_hat_;
  RCP<Thyra::MultiVectorBase<Scalar> > D_g_hat_D_p_;

};


/** \brief Non-member constructor.
 *
 * \relates ForwardResponseSensitivityComputerObserver
 */
template<class Scalar>
RCP<ForwardResponseSensitivityComputerObserver<Scalar> >
forwardResponseSensitivityComputerObserver()
{
  RCP<ForwardResponseSensitivityComputerObserver<Scalar> >
    frsco(new ForwardResponseSensitivityComputerObserver<Scalar>());
  return frsco;
}


/** \brief Non-member constructor.
 *
 * \relates ForwardResponseSensitivityComputerObserver
 */
template<class Scalar>
RCP<ForwardResponseSensitivityComputerObserver<Scalar> >
forwardResponseSensitivityComputerObserver(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const int p_index,
  const int g_index
  )
{
  RCP<ForwardResponseSensitivityComputerObserver<Scalar> >
    frsco = Rythmos::forwardResponseSensitivityComputerObserver<Scalar>();
  frsco->initialize(responseFunc,basePoint,p_index,g_index);
  return frsco;
}


//
// Implementations
//


// Constructors/Initializers/Accessors


template<class Scalar>
ForwardResponseSensitivityComputerObserver<Scalar>::ForwardResponseSensitivityComputerObserver()
{}


template<class Scalar>
void ForwardResponseSensitivityComputerObserver<Scalar>::initialize(
  const RCP<const Thyra::ModelEvaluator<Scalar> > &responseFunc,
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &basePoint,
  const int p_index,
  const int g_index
  )
{
  forwardResponseSensitivityComputer_.setResponseFunction(
    responseFunc, basePoint, p_index, g_index );
  g_hat_ = forwardResponseSensitivityComputer_.create_g_hat();
  D_g_hat_D_p_ = forwardResponseSensitivityComputer_.create_D_g_hat_D_p();
}


template<class Scalar>
const Array<ResponseAndFwdSensPoint<Scalar> >&
ForwardResponseSensitivityComputerObserver<Scalar>::responseAndFwdSensPoints() const
{
  return responseAndFwdSensPoints_;
}


// Overridden from IntegrationObserverBase


template<class Scalar>
RCP<IntegrationObserverBase<Scalar> >
ForwardResponseSensitivityComputerObserver<Scalar>::cloneIntegrationObserver() const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return Teuchos::null;                       
}


template<class Scalar>
void ForwardResponseSensitivityComputerObserver<Scalar>::resetIntegrationObserver(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  responseAndFwdSensPoints_.clear();
}


template<class Scalar>
void ForwardResponseSensitivityComputerObserver<Scalar>::observeCompletedTimeStep(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfo,
  const int timeStepIter
  )
{

  using Teuchos::OSTab;
  using Teuchos::includesVerbLevel;

  // Setup the output info

  const RCP<FancyOStream> out = this->getOStream();
  const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();

  const bool trace = 
    ( !is_null(out) && includesVerbLevel(verbLevel,Teuchos::VERB_LOW) );

  forwardResponseSensitivityComputer_.setOStream(out);
  forwardResponseSensitivityComputer_.setVerbLevel(verbLevel);

  OSTab tab(out);

  if (trace)
    *out << "\nEntering ForwardResponseSensitivityComputerObserver<Scalar>::observeCompletedTimeStep(...) ...\n";

  // A) Get x_bar and x_dot_bar for this time step

  const Scalar t = stepper.getStepStatus().time;

  RCP<const Thyra::VectorBase<Scalar> > x_bar, x_bar_dot;

  get_x_and_x_dot( stepper, t, &x_bar, &x_bar_dot );

  RCP<const Thyra::VectorBase<Scalar> > x;
  RCP<const Thyra::MultiVectorBase<Scalar> > S;
  RCP<const Thyra::VectorBase<Scalar> >  x_dot;
  RCP<const Thyra::MultiVectorBase<Scalar> > S_dot;

  extractStateAndSens( x_bar, x_bar_dot, &x, &S, &x_dot, &S_dot );

  // B) Compute and assemble the response and reduced sensitivity

  forwardResponseSensitivityComputer_.computeResponseAndSensitivity(
    x_dot.get(), S_dot.get(), *x, *S, t, &*g_hat_, &*D_g_hat_D_p_
    );
  
  // C) Store this point

  responseAndFwdSensPoints_.push_back(
    ResponseAndFwdSensPoint<Scalar>(
      t, g_hat_->clone_v(), D_g_hat_D_p_->clone_mv()
      )
    );
  
  if (trace)
    *out << "\nEntering ForwardResponseSensitivityComputerObserver<Scalar>::observeCompletedTimeStep(...) ...\n";

}


} // namespace Rythmos


#endif //RYTHMOS_FORWARD_RESPONSE_SENSITIVITY_COMPUTER_OBSERVER_HPP
