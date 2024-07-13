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

#include "LOCA_AnasaziOperator_Factory.H"
#include "LOCA_AnasaziOperator_AbstractStrategy.H"
#include "LOCA_AnasaziOperator_JacobianInverse.H"
#include "LOCA_AnasaziOperator_ShiftInvert.H"
#include "LOCA_AnasaziOperator_ShiftInvert2Matrix.H"
#include "LOCA_AnasaziOperator_Cayley.H"
#include "LOCA_AnasaziOperator_Cayley2Matrix.H"

LOCA::AnasaziOperator::Factory::Factory(
            const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
}

LOCA::AnasaziOperator::Factory::~Factory()
{
}

Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>
LOCA::AnasaziOperator::Factory::create(
       const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
       const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
       const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
       const Teuchos::RCP<NOX::Abstract::Group>& grp)
{
  std::string methodName = "LOCA::AnasaziOperator::Factory::create()";
  Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy> strategy;

  // Get name of strategy
  const std::string& name = strategyName(*eigenParams);

  if (name == "Jacobian Inverse")
    strategy =
      Teuchos::rcp(new LOCA::AnasaziOperator::JacobianInverse(globalData,
                                  topParams,
                                  eigenParams,
                                  solverParams,
                                  grp));
  else if (name == "Shift-Invert") {
    Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp =
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
    methodName,
    std::string("Group argument for Shift-Invert Anasazi operator ") +
    std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy =
      Teuchos::rcp(new LOCA::AnasaziOperator::ShiftInvert(globalData,
                              topParams,
                              eigenParams,
                              solverParams,
                              tdGrp));
  }
  else if (name == "Shift-Invert 2 Matrix") {
    Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp =
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
    methodName,
    std::string("Group argument for Shift-Invert 2 Matrix Anasazi operator ") +
    std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy =
      Teuchos::rcp(new LOCA::AnasaziOperator::ShiftInvert2Matrix(globalData,
                              topParams,
                              eigenParams,
                              solverParams,
                              tdGrp));
  }
  else if (name == "Cayley") {
    Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp =
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
    methodName,
    std::string("Group argument for Cayley Anasazi operator ") +
    std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy =
      Teuchos::rcp(new LOCA::AnasaziOperator::Cayley(globalData,
                             topParams,
                             eigenParams,
                             solverParams,
                             tdGrp));
  }
  else if (name == "Cayley 2 Matrix") {
    Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp =
      Teuchos::rcp_dynamic_cast<LOCA::TimeDependent::AbstractGroup>(grp);
    if (tdGrp == Teuchos::null)
      globalData->locaErrorCheck->throwError(
    methodName,
    std::string("Group argument for Cayley 2 Matrix Anasazi operator ") +
    std::string("strategy must be a LOCA::TimeDependent::AbstractGroup."));
    strategy =
      Teuchos::rcp(new LOCA::AnasaziOperator::Cayley2Matrix(globalData,
                             topParams,
                             eigenParams,
                             solverParams,
                             tdGrp));
  }
  else if (name == "User-Defined") {

    // Get name of user-defined strategy
    std::string userDefinedName =
      eigenParams->get("Operator User-Defined Name", "???");
    if ((*eigenParams).INVALID_TEMPLATE_QUALIFIER
    isType< Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy> >(userDefinedName))
      strategy = (*eigenParams).INVALID_TEMPLATE_QUALIFIER
    get< Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy> >(userDefinedName);
    else
       globalData->locaErrorCheck->throwError(
                       methodName,
                       "Cannot find user-defined strategy: " +
                       userDefinedName);
  }
  else
    globalData->locaErrorCheck->throwError(
                      methodName,
                      "Invalid Anasazi operator strategy: " +
                      name);

  return strategy;
}

const std::string&
LOCA::AnasaziOperator::Factory::strategyName(
                  Teuchos::ParameterList& eigenParams) const
{
  return eigenParams.get("Operator", "Jacobian Inverse");
}
