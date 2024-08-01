//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Tempus_config.hpp"
#include "Tempus_IntegratorPseudoTransientForwardSensitivity.hpp"
#include "Tempus_IntegratorPseudoTransientAdjointSensitivity.hpp"

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "../TestModels/SteadyQuadraticModel.hpp"

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultProductVector.hpp"

#include <vector>
#include <fstream>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
void test_pseudotransient_fsa(const bool use_dfdp_as_tangent,
                              Teuchos::FancyOStream &out, bool &success)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_DIRK_SteadyQuadratic.xml");

  // Setup the SteadyQuadraticModel
  RCP<ParameterList> scm_pl = sublist(pList, "SteadyQuadraticModel", true);
  scm_pl->set("Use DfDp as Tangent", use_dfdp_as_tangent);
  RCP<SteadyQuadraticModel<double> > model =
      Teuchos::rcp(new SteadyQuadraticModel<double>(scm_pl));

  // Setup sensitivities
  RCP<ParameterList> pl  = sublist(pList, "Tempus", true);
  ParameterList &sens_pl = pl->sublist("Sensitivities");
  sens_pl.set("Use DfDp as Tangent", use_dfdp_as_tangent);
  sens_pl.set("Reuse State Linear Solver", true);
  sens_pl.set("Force W Update", true);  // Have to do this because for this
  // model the solver seems to be overwriting the matrix

  // Setup the Integrator
  RCP<Tempus::IntegratorPseudoTransientForwardSensitivity<double> > integrator =
      Tempus::createIntegratorPseudoTransientForwardSensitivity<double>(pl,
                                                                        model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus);

  // Test if at 'Final Time'
  double time      = integrator->getTime();
  double timeFinal = pl->sublist("Default Integrator")
                         .sublist("Time Step Control")
                         .get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution
  RCP<const Thyra::VectorBase<double> > x_vec         = integrator->getX();
  RCP<const Thyra::MultiVectorBase<double> > DxDp_vec = integrator->getDxDp();
  const double x                                      = Thyra::get_ele(*x_vec, 0);
  const double dxdb                                   = Thyra::get_ele(*(DxDp_vec->col(0)), 0);
  const double x_exact                                = model->getSteadyStateSolution();
  const double dxdb_exact                             = model->getSteadyStateSolutionSensitivity();

  TEST_FLOATING_EQUALITY(x, x_exact, 1.0e-6);
  TEST_FLOATING_EQUALITY(dxdb, dxdb_exact, 1.0e-6);
}

TEUCHOS_UNIT_TEST(DIRK, SteadyQuadratic_PseudoTransient_FSA)
{
  test_pseudotransient_fsa(false, out, success);
}

TEUCHOS_UNIT_TEST(DIRK, SteadyQuadratic_PseudoTransient_FSA_Tangent)
{
  test_pseudotransient_fsa(true, out, success);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(DIRK, SteadyQuadratic_PseudoTransient_ASA)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_DIRK_SteadyQuadratic.xml");

  // Setup the SteadyQuadraticModel
  RCP<ParameterList> scm_pl = sublist(pList, "SteadyQuadraticModel", true);
  RCP<SteadyQuadraticModel<double> > model =
      Teuchos::rcp(new SteadyQuadraticModel<double>(scm_pl));

  // Setup sensitivities
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  // ParameterList& sens_pl = pl->sublist("Sensitivities");

  // Setup the Integrator
  RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<double> > integrator =
      Tempus::integratorPseudoTransientAdjointSensitivity<double>(pl, model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus);

  // Test if at 'Final Time'
  double time      = integrator->getTime();
  double timeFinal = pl->sublist("Default Integrator")
                         .sublist("Time Step Control")
                         .get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Time-integrated solution and the exact solution using the fact that g = x
  RCP<const Thyra::VectorBase<double> > x_vec         = integrator->getX();
  RCP<const Thyra::MultiVectorBase<double> > DxDp_vec = integrator->getDgDp();
  const double x                                      = Thyra::get_ele(*x_vec, 0);
  const double dxdb                                   = Thyra::get_ele(*(DxDp_vec->col(0)), 0);
  const double x_exact                                = model->getSteadyStateSolution();
  const double dxdb_exact                             = model->getSteadyStateSolutionSensitivity();

  TEST_FLOATING_EQUALITY(x, x_exact, 1.0e-6);
  TEST_FLOATING_EQUALITY(dxdb, dxdb_exact, 1.0e-6);
}

}  // namespace Tempus_Test
