//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER
#include "Tempus_PhiEvaluatorLeja.hpp"
#include "Tempus_PhiEvaluatorTaylor.hpp"
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
#include <chrono>
#include <ctime>
#include <sstream>

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
  const int expansion_order = 20;
  RCP<ParameterList> phi_pl = sublist(pList, "PhiEvaluator");
  // use the H-factorization method by default
  const int dd_method = 2;
  phi_pl->set("Leja DD Method", dd_method);
  phi_pl->set("Expansion Order", expansion_order);
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
  Teuchos::ArrayRCP<double> lp_dd = phiEvaluator->getDividedDiffs(0, 1.0, expansion_order);
  // the first divided diff is: y(x_0) == exp(0.0) == 1.0
  // the second divided diff is: (y(x_1) - y(x_0) / x_1 - x_0)) == (exp(-1.0) - 1.0) / (-1.0 - 0.0)
  std::cout << "lp_dd 0: " << lp_dd[0] << std::endl;
  std::cout << "lp_dd 1: " << lp_dd[1] << std::endl;
  std::cout << "lp_dd 2: " << lp_dd[2] << std::endl;
  std::cout << "lp_dd 3: " << lp_dd[3] << std::endl;
  TEST_FLOATING_EQUALITY(lp_dd[0], std::exp(leja_b), 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[1], 0.316060279414, 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[2], 0.075829495185, 1e-8);
  TEST_FLOATING_EQUALITY(lp_dd[3], 0.01263699560, 1e-8);

  // test large ellipse and runtime
  const int expansion_order_high = 300;
  leja_a                         = -600.0;
  leja_c                         = 300.0;
  // test dt scaling: scale Leja ellipse in one case, use dt as argument in another
  double dt                = 0.5;

  phiEvaluator->setLejaEllipse(leja_a, leja_b, leja_c);
  phiEvaluator->setExpansionOrder(expansion_order_high);

  // update for new ellipse
  std::ostringstream ss;
  Teuchos::RCP<Teuchos::Time> localTimer;
  ss << "running baseline dd_phi method " << dd_method << " with dt=" << dt;
  localTimer = Teuchos::TimeMonitor::getNewCounter(ss.str());
  {
    Teuchos::TimeMonitor localTimeMonitor(*localTimer);
    lp_dd = phiEvaluator->getDividedDiffs(0, dt, expansion_order_high);
  }

  // print the entire dd array
  //std::cout << "LP_DD: " << lp_dd() << std::endl;
  std::cout << "First 10 dd: " << lp_dd(lp_dd.lowerOffset(), 10) << std::endl;
  std::cout << "Last 10 dd: " << lp_dd(lp_dd.upperOffset()-10, 10) << std::endl;

  for (int dd_m = 0; dd_m <= 3; dd_m++)
  {
    // check divided difference calculation against other versions

    phiEvaluator->setLejaEllipse(dt * leja_a, dt * leja_b, dt * leja_c);
    phiEvaluator->setDivideDifferenceMethod(dd_m);
    Teuchos::ArrayRCP<double> lp_dd_alt;

    Teuchos::RCP<Teuchos::Time> localTimer2;
    ss.str("");
    ss << "running dd_phi method " << dd_m << " with dt=1, and scaled ellipse";
    localTimer2 = Teuchos::TimeMonitor::getNewCounter(ss.str());
    {
      Teuchos::TimeMonitor localTimeMonitor(*localTimer2);
      lp_dd_alt = phiEvaluator->getDividedDiffs(0, 1.0, expansion_order_high);
    }

    if (dd_m == 0)
    {
      // compare only first 90 entries due to instability of recurrence relation
      TEST_COMPARE_FLOATING_ARRAYS(lp_dd.view(lp_dd.lowerOffset(), 90), lp_dd_alt.view(lp_dd.lowerOffset(), 90), 1e-10);
    }
    else
    {
      TEST_COMPARE_FLOATING_ARRAYS(lp_dd, lp_dd_alt, 1e-10);
    }
  }

  leja_a = -1.0e-18;
  leja_c = 1.0;
  phiEvaluator->setLejaEllipse(leja_a, leja_b, leja_c);
  phiEvaluator->setExpansionOrder(expansion_order);

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
  dt = 1.0;
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

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
