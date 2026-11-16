//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_AdjointAuxSensitivityModelEvaluator.hpp"
#include "Tempus_SolutionHistory.hpp"

#include "../TestModels/SinCosModel.hpp"

namespace Tempus_Unit_Test {

using Tempus_Test::SinCosModel;
using Tempus_Test::SinCosModelAdjoint;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

typedef Thyra::ModelEvaluatorBase MEB;
typedef Thyra::VectorBase<double> Vector;
typedef Tempus::AdjointAuxSensitivityModelEvaluator<double> AuxME;

// Model that counts the evaluations of df/dx (W with alpha = 0, beta = 1)
// and of df/dp of the underlying model
class CountingModel : public Thyra::ModelEvaluatorDelegatorBase<double> {
 public:
  CountingModel(const RCP<const Thyra::ModelEvaluator<double> >& model)
    : num_dfdx(0), num_dfdp(0)
  {
    this->initialize(model);
  }
  mutable int num_dfdx;
  mutable int num_dfdp;

 private:
  void evalModelImpl(const MEB::InArgs<double>& inArgs,
                     const MEB::OutArgs<double>& outArgs) const
  {
    if (nonnull(outArgs.get_W_op()) && inArgs.get_alpha() == 0.0 &&
        inArgs.get_beta() == 1.0)
      ++num_dfdx;
    if (!outArgs.supports(MEB::OUT_ARG_DfDp, 0).none() &&
        !outArgs.get_DfDp(0).isEmpty())
      ++num_dfdp;
    this->getUnderlyingModel()->evalModel(inArgs, outArgs);
  }
};

// SinCos model, its adjoint, and a forward solution history with the exact
// solution at t = 0, 0.5, 1
struct SinCosSetup {
  RCP<SinCosModel<double> > sincos;
  RCP<CountingModel> model;
  RCP<CountingModel> adjoint_model;
  RCP<Tempus::SolutionHistory<double> > sh;

  SinCosSetup()
  {
    RCP<ParameterList> scm_pl = Teuchos::parameterList();
    scm_pl->set("Accept model parameters", true);
    scm_pl->set("Provide nominal values", true);
    scm_pl->set("Coeff a", 0.5);
    scm_pl->set("Coeff f", 2.0);
    scm_pl->set("Coeff L", 1.5);
    sincos        = rcp(new SinCosModel<double>(scm_pl));
    model         = rcp(new CountingModel(sincos));
    adjoint_model = rcp(new CountingModel(rcp(new SinCosModelAdjoint<double>(scm_pl))));

    sh = rcp(new Tempus::SolutionHistory<double>());
    sh->setStorageTypeString("Unlimited");
    const double times[3] = {0.0, 0.5, 1.0};
    for (int i = 0; i < 3; ++i) {
      MEB::InArgs<double> exact = sincos->getExactSolution(times[i]);
      RCP<Tempus::SolutionState<double> > state =
          Tempus::createSolutionStateX<double>(exact.get_x()->clone_v(),
                                               exact.get_x_dot()->clone_v());
      state->setTime(times[i]);
      state->setIndex(i);
      sh->addState(state);
    }
    sh->initialize();
  }

  RCP<AuxME> createAuxME(const std::string& responseType) const
  {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set("Response Type", responseType);
    RCP<AuxME> aux = rcp(new AuxME(model, adjoint_model, 0.0, 1.0, pl));
    aux->setForwardSolutionHistory(sh);
    return aux;
  }
};

// The inArgs of the adjoint model at the adjoint time t, with random adjoint
// state (y, z) and time derivative (x_dot null if !set_x_dot)
MEB::InArgs<double> randomInArgs(const AuxME& aux, const double t,
                                 const double alpha, const double beta,
                                 const bool set_x_dot = true)
{
  MEB::InArgs<double> inArgs = aux.getNominalValues();
  RCP<Vector> x              = Thyra::createMember(aux.get_x_space());
  Thyra::randomize(-1.0, 1.0, x.ptr());
  inArgs.set_x(x);
  if (set_x_dot) {
    RCP<Vector> x_dot = Thyra::createMember(aux.get_x_space());
    Thyra::randomize(-1.0, 1.0, x_dot.ptr());
    inArgs.set_x_dot(x_dot);
  }
  else
    inArgs.set_x_dot(Teuchos::null);
  inArgs.set_t(t);
  inArgs.set_alpha(alpha);
  inArgs.set_beta(beta);
  return inArgs;
}

RCP<Vector> evalResidual(const AuxME& aux, const MEB::InArgs<double>& inArgs)
{
  RCP<Vector> f                = Thyra::createMember(aux.get_f_space());
  MEB::OutArgs<double> outArgs = aux.createOutArgs();
  outArgs.set_f(f);
  aux.evalModel(inArgs, outArgs);
  return f;
}

RCP<Thyra::LinearOpBase<double> > evalW(const AuxME& aux,
                                        const MEB::InArgs<double>& inArgs)
{
  RCP<Thyra::LinearOpBase<double> > W = aux.create_W_op();
  MEB::OutArgs<double> outArgs        = aux.createOutArgs();
  outArgs.set_W_op(W);
  aux.evalModel(inArgs, outArgs);
  return W;
}

// ************************************************************
// df/dx and df/dp are evaluated once for each time, parameters and form
// (x_dot given or not), and the residual does not change.
TEUCHOS_UNIT_TEST(AdjointAuxSensitivityModelEvaluator, DerivativesReused)
{
  SinCosSetup setup;
  RCP<AuxME> aux    = setup.createAuxME("Time Integral");
  CountingModel& fm = *setup.model;
  CountingModel& am = *setup.adjoint_model;

  // Residual and W at t = 0.5
  MEB::InArgs<double> inArgs = randomInArgs(*aux, 0.5, 3.0, 0.7);
  RCP<Vector> F0             = evalResidual(*aux, inArgs);
  TEST_EQUALITY(am.num_dfdx, 1);
  TEST_EQUALITY(fm.num_dfdp, 1);
  MEB::InArgs<double> inArgs2 = randomInArgs(*aux, 0.5, 3.0, 0.7);
  evalResidual(*aux, inArgs2);
  evalW(*aux, inArgs2);
  TEST_EQUALITY(am.num_dfdx, 1);
  TEST_EQUALITY(fm.num_dfdp, 1);

  // Another time
  MEB::InArgs<double> inArgs3 = randomInArgs(*aux, 0.0, 3.0, 0.7);
  evalW(*aux, inArgs3);
  evalResidual(*aux, inArgs3);
  TEST_EQUALITY(am.num_dfdx, 2);
  TEST_EQUALITY(fm.num_dfdp, 2);

  // Other parameter values (changed in place), then the nominal ones
  RCP<Vector> p = inArgs3.get_p(0)->clone_v();
  inArgs3.set_p(0, p);
  Thyra::scale(1.1, p.ptr());
  evalResidual(*aux, inArgs3);
  TEST_EQUALITY(am.num_dfdx, 3);
  TEST_EQUALITY(fm.num_dfdp, 3);
  Thyra::assign(p.ptr(), *inArgs.get_p(0));
  evalResidual(*aux, inArgs3);
  TEST_EQUALITY(am.num_dfdx, 4);
  TEST_EQUALITY(fm.num_dfdp, 4);

  // Explicit form (no x_dot)
  MEB::InArgs<double> inArgs4 = randomInArgs(*aux, 0.0, 3.0, 0.7, false);
  evalResidual(*aux, inArgs4);
  TEST_EQUALITY(am.num_dfdx, 5);
  TEST_EQUALITY(fm.num_dfdp, 5);

  // The residual is the one of a new model evaluator
  RCP<Vector> F1   = evalResidual(*aux, inArgs);
  RCP<Vector> F2   = evalResidual(*setup.createAuxME("Time Integral"), inArgs);
  RCP<Vector> diff = F1->clone_v();
  Thyra::Vp_StV(diff.ptr(), -1.0, *F0);
  TEST_COMPARE(Thyra::norm_2(*diff), <=, 1.0e-14 * Thyra::norm_2(*F0));
  Thyra::V_VmV(diff.ptr(), *F2, *F0);
  TEST_COMPARE(Thyra::norm_2(*diff), <=, 1.0e-14 * Thyra::norm_2(*F0));
}

}  // namespace Tempus_Unit_Test
