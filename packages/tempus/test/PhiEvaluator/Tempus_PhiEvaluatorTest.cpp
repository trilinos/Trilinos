//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
#include "Tempus_PhiEvaluatorLeja_decl.hpp"
#include "Tempus_PhiEvaluatorTaylor_decl.hpp"
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
#include "Thyra_VectorStdOps_decl.hpp"

#include <cmath>

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
  RCP<ParameterList> pListTay =
      getParametersFromXmlFile("Tempus_PhiEvaluator_SinCos.xml");

  // Setup the SinCosModel as mock model
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model = rcp(new SinCosModel<double>(scm_pl));
  // auto model_tay = rcp(new SinCosModel<double>(scm_pl));

  // Setup the PhiEvaluator
  RCP<ParameterList> phi_pl = sublist(pList, "PhiEvaluator");
  auto phiEvaluator = Tempus::createPhiEvaluatorLeja<double>(phi_pl);
  phiEvaluator->setModel(model);
  phiEvaluator->initialize();

  // Setup Taylor PhiEvaluator for comparison
  RCP<ParameterList> phi_pl_tay = sublist(pListTay, "PhiEvaluator");
  auto phiEvaluatorTay = Tempus::createPhiEvaluatorTaylor<double>(phi_pl_tay);
  phiEvaluatorTay->setModel(model);
  phiEvaluatorTay->initialize();

  double leja_a = -1.0;
  double leja_b = 0.0;
  double leja_c = 0.5;
  phiEvaluator->setLejaEllipse(leja_a, leja_b, leja_c);

  // Check the first leja points with scaling
  LejaPoint lp = phiEvaluator->getLpSc(0);
  TEST_ASSERT(lp.lpt == LpType::LPREAL);
  TEST_FLOATING_EQUALITY(lp.lp.real(), 0.0, 1e-6);
  TEST_FLOATING_EQUALITY(lp.lp.imag(), 0.0, 1e-6);
  lp = phiEvaluator->getLpSc(1);
  TEST_ASSERT(lp.lpt == LpType::LPREAL);
  TEST_FLOATING_EQUALITY(lp.lp.real(), leja_a, 1e-6);
  TEST_FLOATING_EQUALITY(lp.lp.imag(), 0.0, 1e-6);
  lp = phiEvaluator->getLpSc(2);
  TEST_ASSERT(lp.lpt == LpType::LPCONJ);
  TEST_FLOATING_EQUALITY(lp.lp.real(), 0.5 * leja_a, 1e-6);
  TEST_FLOATING_EQUALITY(lp.lp.imag(), leja_c, 1e-6);

  // Check the first divided diffs
  const int exp_order = 4;
  auto lp_dd = phiEvaluator->getDividedDiffs(0, 1.0, exp_order);
  // the first divided diff is: y(x_0) == exp(0.0) == 1.0
  // the second divided diff is: (y(x_1) - y(x_0) / x_1 - x_0)) == (exp(-1.0) - 1.0) / (-1.0 - 0.0)
  std::cout << "lp_dd 0: " << lp_dd[0] << std::endl;
  std::cout << "lp_dd 1: " << lp_dd[1] << std::endl;
  TEST_FLOATING_EQUALITY(lp_dd[0].real(), std::exp(leja_b), 1e-8);
  // TEST_FLOATING_EQUALITY(lp_dd[0].imag(), 0.0, 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[1].real(), 0.316060279414, 1e-8);
  // TEST_FLOATING_EQUALITY(lp_dd[1].imag(), 0.0, 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[2].real(), 0.075829495185, 1e-8);
  // TEST_FLOATING_EQUALITY(lp_dd[2].imag(), 0.012636995600793, 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[3].real(), 0.01263699560, 1e-8);
  // TEST_FLOATING_EQUALITY(lp_dd[3].imag(), 0.0, 1e-8);
  //TODO, do not test imaginary Leja dds, not needed

  leja_a = -1.0e-18;
  leja_c = 1.0;
  phiEvaluator->setLejaEllipse(leja_a, leja_b, leja_c);

  // compute exp(dt*A)*v using PhiEvaluatorLeja
  // make a digonal linop from SinCosModel
  // with A = [[0, -1], [0, 1]]
  auto x_space = model->get_x_space();
  Teuchos::RCP<Thyra::VectorBase<double>> xdot_init = createMember(x_space);
  Teuchos::RCP<Thyra::VectorBase<double>> v = createMember(x_space);
  Teuchos::RCP<Thyra::VectorBase<double>> vend = createMember(x_space);
  Teuchos::RCP<Thyra::VectorBase<double>> vend_tay = createMember(x_space);

  // set initial condition v0 = (-1.0, 0)
  Thyra::assign(v.ptr(), 0.0);
  Thyra::set_ele(0, -1.0, v.ptr());
  Thyra::assign(xdot_init.ptr(), 0.0);

  auto inArgs = model->createInArgs();
  inArgs.set_x(v);
  inArgs.set_x_dot(xdot_init);
  inArgs.set_t(0.0);
  double dt = 1.0;
  phiEvaluator->setLinearizationPoint(inArgs);
  phiEvaluator->computePhi(vend.ptr(), 0, dt, v);

  Teuchos::RCP<Teuchos::FancyOStream> ostream = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  vend->describe(*ostream, Teuchos::VERB_EXTREME);

  // compare to Taylor PhiEvaluator
  phiEvaluatorTay->setLinearizationPoint(inArgs);
  phiEvaluatorTay->computePhi(vend_tay.ptr(), 0, dt, v);
  vend_tay->describe(*ostream, Teuchos::VERB_EXTREME);
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*vend_tay, 0), Thyra::get_ele(*vend, 0), 1e-6);
  TEST_FLOATING_EQUALITY(Thyra::get_ele(*vend_tay, 1), Thyra::get_ele(*vend, 1), 1e-6);

  // compare to analytic solution
  // v(t) = [-cos(t), sin(t)]
  std::vector<double> sol = {-std::cos(dt), std::sin(dt)};
  TEST_FLOATING_EQUALITY(sol[0], Thyra::get_ele(*vend, 0), 1e-6);
  TEST_FLOATING_EQUALITY(sol[1], Thyra::get_ele(*vend, 1), 1e-6);
}

}  // namespace Tempus_Test
