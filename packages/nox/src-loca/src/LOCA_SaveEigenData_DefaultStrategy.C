// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_SaveEigenData_DefaultStrategy.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"


LOCA::SaveEigenData::DefaultStrategy::DefaultStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) :
  globalData(global_data)
{
}

LOCA::SaveEigenData::DefaultStrategy::~DefaultStrategy()
{
}

NOX::Abstract::Group::ReturnType
LOCA::SaveEigenData::DefaultStrategy::save(
         Teuchos::RCP< std::vector<double> >& /* evals_r */,
         Teuchos::RCP< std::vector<double> >& /* evals_i */,
         Teuchos::RCP< NOX::Abstract::MultiVector >& /* evecs_r */,
             Teuchos::RCP< NOX::Abstract::MultiVector >& /* evecs_i */)
{
  return NOX::Abstract::Group::Ok;
}
