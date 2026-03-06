//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
#include "Tempus_PhiEvaluatorLeja_decl.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DetachedVectorView.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorLeja.hpp"

#include "../TestModels/SinCosModel.hpp"


namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::sublist;

using Tempus::LejaPoint;
using Tempus::LpType;

TEUCHOS_UNIT_TEST(PhiEvaluator, Leja_SinCos)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_PhiEvaluator_SinCos.xml");

  // Setup the SinCosModel as mock model
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model = rcp(new SinCosModel<double>(scm_pl));

  // Setup the PhiEvaluator
  RCP<ParameterList> phi_pl = sublist(pList, "PhiEvaluator");
  auto phiEvaluator = Tempus::createPhiEvaluatorLeja<double>(phi_pl);
  phiEvaluator->setModel(model);
  phiEvaluator->initialize();
  phiEvaluator->setLejaEllipse(-1.0, 0.0, 0.5);

  // Check the first leja points
  LejaPoint lp = phiEvaluator->getLpSc(0);
  TEST_ASSERT(lp.lpt == LpType::LPREAL);
  TEST_FLOATING_EQUALITY(lp.get().at(0).real(), 0.0, 1e-6);
  TEST_FLOATING_EQUALITY(lp.get().at(0).imag(), 0.0, 1e-6);
  lp = phiEvaluator->getLpSc(1);
  TEST_ASSERT(lp.lpt == LpType::LPREAL);
  TEST_FLOATING_EQUALITY(lp.get().at(0).real(), -1.0, 1e-6);
  TEST_FLOATING_EQUALITY(lp.get().at(0).imag(), 0.0, 1e-6);
  lp = phiEvaluator->getLpSc(2);
  TEST_ASSERT(lp.lpt == LpType::LPCONJ);
  TEST_FLOATING_EQUALITY(lp.get().at(0).real(), 0.5, 1e-6);
  TEST_FLOATING_EQUALITY(lp.get().at(0).imag(), 0.5, 1e-6);
  TEST_FLOATING_EQUALITY(lp.get().at(1).real(), 0.5, 1e-6);
  TEST_FLOATING_EQUALITY(lp.get().at(1).imag(), -0.5, 1e-6);

  // Check the first divided diffs
  auto lp_dd = phiEvaluator->getDividedDiffs(0, 1.0);
}

}  // namespace Tempus_Test
