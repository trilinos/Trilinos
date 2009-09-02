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

#include "Rythmos_Charon_IntegrationControlAndObserver.hpp"
#include "Rythmos_Charon_Solver.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"


namespace RythmosCharon {


CharonIntegrationControlAndObserver::CharonIntegrationControlAndObserver() { }
CharonIntegrationControlAndObserver::~CharonIntegrationControlAndObserver() { }
// Overridden from Rythmos::IntegrationControlStrategyBase
Teuchos::RCP<Rythmos::IntegrationControlStrategyBase<double> > CharonIntegrationControlAndObserver::cloneIntegrationControlStrategy() const
{
  return Teuchos::null;
}
void CharonIntegrationControlAndObserver::resetIntegrationControlStrategy(
  const Rythmos::TimeRange<double> &integrationTimeDomain
  )
{ }
Rythmos::StepControlInfo<double> CharonIntegrationControlAndObserver::getNextStepControlInfo(
  const Rythmos::StepperBase<double> &stepper,
  const Rythmos::StepControlInfo<double> &stepCtrlInfoLast,
  const int timeStepIter
  )
{
  Rythmos::StepControlInfo<double> sci;
  return sci;
}
// Overridden from Rythmos::IntegrationObserverBase
Teuchos::RCP<Rythmos::IntegrationObserverBase<double> > CharonIntegrationControlAndObserver::cloneIntegrationObserver() const
{
  return Teuchos::null;
}
void CharonIntegrationControlAndObserver::resetIntegrationObserver(
  const Rythmos::TimeRange<double> &integrationTimeDomain
  )
{ }
void CharonIntegrationControlAndObserver::observeCompletedTimeStep(
  const Rythmos::StepperBase<double> &stepper,
  const Rythmos::StepControlInfo<double> &stepCtrlInfo,
  const int timeStepIter
  )
{ }
// Overridden from Teuchos::ParameterListAcceptorDefaultBase
void CharonIntegrationControlAndObserver::setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList )
{ }
Teuchos::RCP<const Teuchos::ParameterList> CharonIntegrationControlAndObserver::getValidParameters() const
{
  return Teuchos::parameterList();
}

} // namespace RythmosCharon


