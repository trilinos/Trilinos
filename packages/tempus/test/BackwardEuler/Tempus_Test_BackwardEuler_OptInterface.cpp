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
#include "Teuchos_DefaultComm.hpp"

#include "Thyra_MultiVectorStdOps.hpp"

#include "Tempus_config.hpp"
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_StepperBackwardEuler.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <vector>
#include <fstream>
#include <sstream>
#include <limits>

namespace Tempus_Test {

using Teuchos::getParametersFromXmlFile;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::sublist;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, OptInterface)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
      getParametersFromXmlFile("Tempus_BackwardEuler_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  auto model                = rcp(new SinCosModel<double>(scm_pl));

  // Setup the Integrator
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
      Tempus::createIntegratorBasic<double>(pl, model);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus);

  // Get solution history
  RCP<const SolutionHistory<double> > solutionHistory =
      integrator->getSolutionHistory();

  // Get the stepper and cast to optimization interface
  RCP<Tempus::Stepper<double> > stepper = integrator->getStepper();
  RCP<Tempus::StepperOptimizationInterface<double> > opt_stepper =
      Teuchos::rcp_dynamic_cast<Tempus::StepperOptimizationInterface<double> >(
          stepper, true);

  // Check stencil length
  TEST_EQUALITY(opt_stepper->stencilLength(), 2);

  // Create needed vectors/multivectors
  Teuchos::Array<RCP<const Thyra::VectorBase<double> > > x(2);
  Teuchos::Array<double> t(2);
  RCP<const Thyra::VectorBase<double> > p = model->getNominalValues().get_p(0);
  RCP<Thyra::VectorBase<double> > x_dot =
      Thyra::createMember(model->get_x_space());
  RCP<Thyra::VectorBase<double> > f = Thyra::createMember(model->get_f_space());
  RCP<Thyra::VectorBase<double> > f2 =
      Thyra::createMember(model->get_f_space());
  RCP<Thyra::LinearOpBase<double> > dfdx  = model->create_W_op();
  RCP<Thyra::LinearOpBase<double> > dfdx2 = model->create_W_op();
  RCP<Thyra::MultiVectorBase<double> > dfdx_mv =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(dfdx, true);
  RCP<Thyra::MultiVectorBase<double> > dfdx_mv2 =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(dfdx2, true);
  const int num_p = p->range()->dim();
  RCP<Thyra::MultiVectorBase<double> > dfdp =
      Thyra::createMembers(model->get_f_space(), num_p);
  RCP<Thyra::MultiVectorBase<double> > dfdp2 =
      Thyra::createMembers(model->get_f_space(), num_p);
  RCP<Thyra::LinearOpWithSolveBase<double> > W  = model->create_W();
  RCP<Thyra::LinearOpWithSolveBase<double> > W2 = model->create_W();
  RCP<Thyra::MultiVectorBase<double> > tmp =
      Thyra::createMembers(model->get_x_space(), num_p);
  RCP<Thyra::MultiVectorBase<double> > tmp2 =
      Thyra::createMembers(model->get_x_space(), num_p);
  std::vector<double> nrms(num_p);
  double err;

  // Loop over states, checking residuals and derivatives
  const int n = solutionHistory->getNumStates();
  for (int i = 1; i < n; ++i) {
    RCP<const SolutionState<double> > state      = (*solutionHistory)[i];
    RCP<const SolutionState<double> > prev_state = (*solutionHistory)[i - 1];

    // Fill x, t stencils
    x[0] = state->getX();
    x[1] = prev_state->getX();
    t[0] = state->getTime();
    t[1] = prev_state->getTime();

    // Compute x_dot
    const double dt = t[0] - t[1];
    Thyra::V_StVpStV(x_dot.ptr(), 1.0 / dt, *(x[0]), -1.0 / dt, *(x[1]));

    // Create model inargs
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<double> in_args   = model->createInArgs();
    MEB::OutArgs<double> out_args = model->createOutArgs();
    in_args.set_x(x[0]);
    in_args.set_x_dot(x_dot);
    in_args.set_t(t[0]);
    in_args.set_p(0, p);

    const double tol = 1.0e-14;

    // Check residual
    opt_stepper->computeStepResidual(*f, x, t, *p, 0);
    out_args.set_f(f2);
    model->evalModel(in_args, out_args);
    out_args.set_f(Teuchos::null);
    Thyra::V_VmV(f.ptr(), *f, *f2);
    err = Thyra::norm(*f);
    TEST_FLOATING_EQUALITY(err, 0.0, tol);

    // Check df/dx_n
    // df/dx_n = df/dx_dot * dx_dot/dx_n + df/dx_n = 1/dt*df/dx_dot + df/dx_n
    opt_stepper->computeStepJacobian(*dfdx, x, t, *p, 0, 0);
    out_args.set_W_op(dfdx2);
    in_args.set_alpha(1.0 / dt);
    in_args.set_beta(1.0);
    model->evalModel(in_args, out_args);
    out_args.set_W_op(Teuchos::null);
    Thyra::V_VmV(dfdx_mv.ptr(), *dfdx_mv, *dfdx_mv2);
    Thyra::norms(*dfdx_mv, Teuchos::arrayViewFromVector(nrms));
    err = 0.0;
    for (auto nrm : nrms) err += nrm;
    TEST_FLOATING_EQUALITY(err, 0.0, tol);

    // Check df/dx_{n-1}
    // df/dx_{n-1} = df/dx_dot * dx_dot/dx_{n-1} = -1/dt*df/dx_dot
    opt_stepper->computeStepJacobian(*dfdx, x, t, *p, 0, 1);
    out_args.set_W_op(dfdx2);
    in_args.set_alpha(-1.0 / dt);
    in_args.set_beta(0.0);
    model->evalModel(in_args, out_args);
    out_args.set_W_op(Teuchos::null);
    Thyra::V_VmV(dfdx_mv.ptr(), *dfdx_mv, *dfdx_mv2);
    Thyra::norms(*dfdx_mv, Teuchos::arrayViewFromVector(nrms));
    err = 0.0;
    for (auto nrm : nrms) err += nrm;
    TEST_FLOATING_EQUALITY(err, 0.0, tol);

    // Check df/dp
    opt_stepper->computeStepParamDeriv(*dfdp, x, t, *p, 0);
    out_args.set_DfDp(
        0, MEB::Derivative<double>(dfdp2, MEB::DERIV_MV_JACOBIAN_FORM));
    model->evalModel(in_args, out_args);
    out_args.set_DfDp(0, MEB::Derivative<double>());
    Thyra::V_VmV(dfdp.ptr(), *dfdp, *dfdp2);
    Thyra::norms(*dfdp, Teuchos::arrayViewFromVector(nrms));
    err = 0.0;
    for (auto nrm : nrms) err += nrm;
    TEST_FLOATING_EQUALITY(err, 0.0, tol);

    // Check W
    opt_stepper->computeStepSolver(*W, x, t, *p, 0);
    out_args.set_W(W2);
    in_args.set_alpha(1.0 / dt);
    in_args.set_beta(1.0);
    model->evalModel(in_args, out_args);
    out_args.set_W(Teuchos::null);
    // note:  dfdp overwritten above so dfdp != dfdp2
    Thyra::solve(*W, Thyra::NOTRANS, *dfdp2, tmp.ptr());
    Thyra::solve(*W2, Thyra::NOTRANS, *dfdp2, tmp2.ptr());
    Thyra::V_VmV(tmp.ptr(), *tmp, *tmp2);
    Thyra::norms(*tmp, Teuchos::arrayViewFromVector(nrms));
    err = 0.0;
    for (auto nrm : nrms) err += nrm;
    TEST_FLOATING_EQUALITY(err, 0.0, tol);
  }

  Teuchos::TimeMonitor::summarize();
}

}  // namespace Tempus_Test
