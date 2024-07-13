// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

#include "LOCA_Bifurcation_Factory.H"
#include "LOCA_TurningPoint_MooreSpence_ExtendedGroup.H"
#include "LOCA_TurningPoint_MooreSpence_AbstractGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Pitchfork_MooreSpence_ExtendedGroup.H"
#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_ExtendedGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_ExtendedGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_PhaseTransition_AbstractGroup.H"
#include "LOCA_PhaseTransition_ExtendedGroup.H"

LOCA::Bifurcation::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::Bifurcation::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
LOCA::Bifurcation::Factory::create(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp)
{
  std::string methodName = "LOCA::Bifurcation::Factory::create()";
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*bifurcationParams);

  if (name == "None")
    strategy = grp;

  else if (name == "Turning Point:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RCP<LOCA::TurningPoint::MooreSpence::AbstractGroup> msg =
      Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
            methodName,
            std::string("Underlying group must be derived from ") +
            std::string("LOCA::TurningPoint::MooreSpence::AbstractGroup ") +
            std::string("for Moore-Spence turning point continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedGroup(
                               globalData,
                               topParams,
                               bifurcationParams,
                               msg));
  }
  else if (name == "Turning Point:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup> mag =
      Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
        methodName,
        std::string("Underlying group must be derived from ") +
        std::string("LOCA::TurningPoint::MinimallyAugmented::AbstractGroup ") +
        std::string("for minimally augmented turning point continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup(
                               globalData,
                               topParams,
                               bifurcationParams,
                               mag));
  }
  else if (name == "Pitchfork:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RCP<LOCA::Pitchfork::MooreSpence::AbstractGroup> msg =
      Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
            methodName,
            std::string("Underlying group must be derived from ") +
            std::string("LOCA::Pitchfork::MooreSpence::AbstractGroup ") +
            std::string("for Moore-Spence pitchfork continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedGroup(
                               globalData,
                               topParams,
                               bifurcationParams,
                               msg));
  }
  else if (name == "Pitchfork:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RCP<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup> mag =
      Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
        methodName,
        std::string("Underlying group must be derived from ") +
        std::string("LOCA::Pitchfork::MinimallyAugmented::AbstractGroup ") +
        std::string("for minimally augmented pitchfork continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup(
                               globalData,
                               topParams,
                               bifurcationParams,
                               mag));
  }
  else if (name == "Hopf:  Moore-Spence") {

    // Cast group to MooreSpence group
    Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup> msg =
      Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::AbstractGroup>(grp);
    if (msg.get() == NULL)
      globalData->locaErrorCheck->throwError(
            methodName,
            std::string("Underlying group must be derived from ") +
            std::string("LOCA::Hopf::MooreSpence::AbstractGroup ") +
            std::string("for Moore-Spence Hopf continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedGroup(
                                 globalData,
                                 topParams,
                                 bifurcationParams,
                                 msg));
  }
  else if (name == "Hopf:  Minimally Augmented") {

    // Cast group to MinimallyAugmented group
    Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup> mag =
      Teuchos::rcp_dynamic_cast<LOCA::Hopf::MinimallyAugmented::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
        methodName,
        std::string("Underlying group must be derived from ") +
        std::string("LOCA::Hopf::MinimallyAugmented::AbstractGroup ") +
        std::string("for minimally augmented Hopf continuation!"));

    strategy =
      Teuchos::rcp(new LOCA::Hopf::MinimallyAugmented::ExtendedGroup(
                               globalData,
                               topParams,
                               bifurcationParams,
                               mag));
  }
  else if (name == "Phase Transition") {

    // Cast group to MinimallyAugmented group
    Teuchos::RCP<LOCA::PhaseTransition::AbstractGroup> mag =
      Teuchos::rcp_dynamic_cast<LOCA::PhaseTransition::AbstractGroup>(grp);
    if (mag.get() == NULL)
      globalData->locaErrorCheck->throwError(
        methodName,
        std::string("Underlying group must be derived from ") +
        std::string("LOCA::PhaseTransition::AbstractGroup ") +
        std::string("for Phase Transition tracking!"));

    strategy =
      Teuchos::rcp(new LOCA::PhaseTransition::ExtendedGroup(
                           globalData,
                           bifurcationParams,
                           mag));
  }
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName = bifurcationParams->get(
                             "User-Defined Name",
                             "???");
    if ((*bifurcationParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> >(userDefinedName))
      strategy = (*bifurcationParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                       methodName,
                       "Cannot find user-defined strategy: " +
                       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid bifurcation method: " +
                      name);

  return strategy;
}

std::string
LOCA::Bifurcation::Factory::strategyName(
                Teuchos::ParameterList& bifurcationParams) const
{
  // Get bifurcation type
  std::string bif_type =  bifurcationParams.get("Type", "None");

  // Get the formulation
  if (bif_type == "Turning Point" || bif_type == "Pitchfork" ||
      bif_type == "Hopf") {
    std::string formulation =
      bifurcationParams.get("Formulation", "Moore-Spence");
    bif_type += ":  " + formulation;
  }

  return bif_type;
}
