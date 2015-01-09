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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_CHARON_INTEGRATION_CONTROL_AND_OBSERVER_HPP
#define RYTHMOS_CHARON_INTEGRATION_CONTROL_AND_OBSERVER_HPP

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Rythmos_StepControlInfo.hpp"
#include "Rythmos_StepperSupportTypes.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"

#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Rythmos_StepControlInfo.hpp"


namespace RythmosCharon {


class CharonIntegrationControlAndObserver
  : virtual public Rythmos::IntegrationControlStrategyBase<double>,
    virtual public Rythmos::IntegrationObserverBase<double>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
  CharonIntegrationControlAndObserver();
  virtual ~CharonIntegrationControlAndObserver();
  // Overridden from Rythmos::IntegrationControlStrategyBase
  Teuchos::RCP<Rythmos::IntegrationControlStrategyBase<double> > cloneIntegrationControlStrategy() const;
  void resetIntegrationControlStrategy(
    const Rythmos::TimeRange<double> &integrationTimeDomain
    );
  Rythmos::StepControlInfo<double> getNextStepControlInfo(
    const Rythmos::StepperBase<double> &stepper,
    const Rythmos::StepControlInfo<double> &stepCtrlInfoLast,
    const int timeStepIter
    );
  // Overridden from Rythmos::IntegrationObserverBase
  Teuchos::RCP<Rythmos::IntegrationObserverBase<double> > cloneIntegrationObserver() const;
  void resetIntegrationObserver(
    const Rythmos::TimeRange<double> &integrationTimeDomain
    );
  void observeCompletedTimeStep(
    const Rythmos::StepperBase<double> &stepper,
    const Rythmos::StepControlInfo<double> &stepCtrlInfo,
    const int timeStepIter
    );
  // Overridden from Teuchos::ParameterListAcceptorDefaultBase
  void setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList );
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
};

} // namespace RythmosCharon

#endif // RYTHMOS_CHARON_INTEGRATION_CONTROL_AND_OBSERVER_HPP

