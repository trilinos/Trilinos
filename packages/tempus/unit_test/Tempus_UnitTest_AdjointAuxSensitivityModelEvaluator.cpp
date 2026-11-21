//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_AdjointAuxCouplingLinearOp.hpp"
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
typedef Thyra::MultiVectorBase<double> MultiVector;
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

// The residual is affine in (x, x_dot), so W*v = alpha*dF/dx_dot*v +
// beta*dF/dx*v = alpha*(F(x, x_dot+v) - F(x, x_dot)) + beta*(F(x+v, x_dot) -
// F(x, x_dot)) up to rounding.
void checkWIsJacobian(const std::string& responseType,
                      Teuchos::FancyOStream& out, bool& success)
{
  SinCosSetup setup;
  RCP<AuxME> aux = setup.createAuxME(responseType);

  const double alpha = 3.0, beta = 0.7;
  for (const double t : {0.0, 0.5}) {
    MEB::InArgs<double> inArgs = randomInArgs(*aux, t, alpha, beta);
    RCP<Vector> v              = Thyra::createMember(aux->get_x_space());
    Thyra::randomize(-1.0, 1.0, v.ptr());

    RCP<Vector> F0               = evalResidual(*aux, inArgs);
    MEB::InArgs<double> inArgs_x = inArgs;
    RCP<Vector> x_plus_v         = inArgs.get_x()->clone_v();
    Thyra::Vp_V(x_plus_v.ptr(), *v);
    inArgs_x.set_x(x_plus_v);
    RCP<Vector> Fx                  = evalResidual(*aux, inArgs_x);
    MEB::InArgs<double> inArgs_xdot = inArgs;
    RCP<Vector> x_dot_plus_v        = inArgs.get_x_dot()->clone_v();
    Thyra::Vp_V(x_dot_plus_v.ptr(), *v);
    inArgs_xdot.set_x_dot(x_dot_plus_v);
    RCP<Vector> Fxdot = evalResidual(*aux, inArgs_xdot);

    // Jv = beta*(Fx - F0) + alpha*(Fxdot - F0)
    RCP<Vector> Jv = Thyra::createMember(aux->get_f_space());
    Thyra::V_StVpStV(Jv.ptr(), beta, *Fx, alpha, *Fxdot);
    Thyra::Vp_StV(Jv.ptr(), -(alpha + beta), *F0);

    RCP<Thyra::LinearOpBase<double> > W = evalW(*aux, inArgs);
    RCP<Vector> Wv                      = Thyra::createMember(aux->get_f_space());
    Thyra::apply(*W, Thyra::NOTRANS, *v, Wv.ptr());

    RCP<Vector> diff = Wv->clone_v();
    Thyra::Vp_StV(diff.ptr(), -1.0, *Jv);
    out << "t = " << t << ": |W*v| = " << Thyra::norm_2(*Wv)
        << ", |W*v - J*v| = " << Thyra::norm_2(*diff) << std::endl;
    TEST_COMPARE(Thyra::norm_2(*diff), <=, 1.0e-12 * Thyra::norm_2(*Jv));
  }
}

// ************************************************************
TEUCHOS_UNIT_TEST(AdjointAuxSensitivityModelEvaluator, W_IsJacobian_FinalTime)
{
  checkWIsJacobian("Final Time", out, success);
}

// ************************************************************
TEUCHOS_UNIT_TEST(AdjointAuxSensitivityModelEvaluator, W_IsJacobian_TimeIntegral)
{
  checkWIsJacobian("Time Integral", out, success);
}

// ************************************************************
// The solver of W from get_W_factory() solves the coupled system.
TEUCHOS_UNIT_TEST(AdjointAuxSensitivityModelEvaluator, W_Solve)
{
  SinCosSetup setup;
  RCP<AuxME> aux                      = setup.createAuxME("Final Time");
  MEB::InArgs<double> inArgs          = randomInArgs(*aux, 0.5, 3.0, 0.7);
  RCP<Thyra::LinearOpBase<double> > W = evalW(*aux, inArgs);
  RCP<Thyra::LinearOpWithSolveBase<double> > W_lows =
      Thyra::linearOpWithSolve<double>(*aux->get_W_factory(), W);

  RCP<Vector> b = Thyra::createMember(aux->get_f_space());
  Thyra::randomize(-1.0, 1.0, b.ptr());
  RCP<Vector> x = Thyra::createMember(aux->get_x_space());
  Thyra::assign(x.ptr(), 0.0);
  Thyra::solve<double>(*W_lows, Thyra::NOTRANS, *b, x.ptr());

  RCP<Vector> r = b->clone_v();
  Thyra::apply(*W, Thyra::NOTRANS, *x, r.ptr(), 1.0, -1.0);
  out << "|W*x - b| = " << Thyra::norm_2(*r) << ", |b| = " << Thyra::norm_2(*b)
      << std::endl;
  TEST_COMPARE(Thyra::norm_2(*r), <=, 1.0e-12 * Thyra::norm_2(*b));
}

// ************************************************************
// The factory of W rejects, when reinitializing, an operator that is not
// block triangular or has another number of blocks, and leaves the solver
// unchanged (it checks the operator before reinitializing the diagonal
// blocks).
TEUCHOS_UNIT_TEST(BlockedTriangularLinearOpWithSolveFactory,
                  ReinitializeRejected)
{
  typedef Thyra::LinearOpBase<double> LinearOp;
  typedef Thyra::PhysicallyBlockedLinearOpBase<double> PBLOB;
  SinCosSetup setup;
  RCP<AuxME> aux = setup.createAuxME("Time Integral");
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > lowsf =
      aux->get_W_factory();
  // The two W differ in their diagonal blocks
  RCP<LinearOp> W  = evalW(*aux, randomInArgs(*aux, 0.5, 3.0, 0.7));
  RCP<LinearOp> W2 = evalW(*aux, randomInArgs(*aux, 0.5, 1.5, 0.4));
  RCP<Thyra::LinearOpWithSolveBase<double> > W_lows =
      Thyra::linearOpWithSolve<double>(*lowsf, W);

  // Relative residual of the solve of op*x = b with W_lows
  auto solveResidual = [&](const LinearOp& op) {
    RCP<Vector> b = Thyra::createMember(op.range());
    Thyra::randomize(-1.0, 1.0, b.ptr());
    RCP<Vector> x = Thyra::createMember(op.domain());
    Thyra::assign(x.ptr(), 0.0);
    Thyra::solve<double>(*W_lows, Thyra::NOTRANS, *b, x.ptr());
    RCP<Vector> r = b->clone_v();
    Thyra::apply(op, Thyra::NOTRANS, *x, r.ptr(), 1.0, -1.0);
    return Thyra::norm_2(*r) / Thyra::norm_2(*b);
  };

  // W2 with a block above the diagonal, and its first diagonal block alone
  RCP<const PBLOB> W2_blo   = Teuchos::rcp_dynamic_cast<const PBLOB>(W2, true);
  RCP<const LinearOp> W2_10 = W2_blo->getBlock(1, 0).assert_not_null();
  RCP<const LinearOp> notTriangular =
      Thyra::block2x2<double>(W2_blo->getBlock(0, 0), Thyra::adjoint(W2_10),
                              W2_10, W2_blo->getBlock(1, 1));
  RCP<const LinearOp> oneBlock = Thyra::block1x1<double>(W2_blo->getBlock(0, 0));
  TEST_THROW(Thyra::initializeOp<double>(*lowsf, notTriangular, W_lows.ptr()),
             std::logic_error);
  TEST_THROW(Thyra::initializeOp<double>(*lowsf, oneBlock, W_lows.ptr()),
             std::logic_error);

  // W_lows still solves W, and can be reinitialized with W2
  TEST_COMPARE(solveResidual(*W), <=, 1.0e-12);
  Thyra::initializeOp<double>(*lowsf, W2, W_lows.ptr());
  TEST_COMPARE(solveResidual(*W2), <=, 1.0e-12);
  TEST_COMPARE(solveResidual(*W), >, 1.0e-6);
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
  evalW(*aux, inArgs3);  // df/dp for the coupling block of W
  TEST_EQUALITY(fm.num_dfdp, 2);
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

// ************************************************************
// A W keeps the df/dp of its evaluation point: later evaluations of another W
// or of the residual at other times (the df/dp of SinCos depends on the
// forward state) do not change it.
TEUCHOS_UNIT_TEST(AdjointAuxSensitivityModelEvaluator, W_KeepsItsEvaluationPoint)
{
  SinCosSetup setup;
  RCP<AuxME> aux = setup.createAuxME("Final Time");
  RCP<Vector> v  = Thyra::createMember(aux->get_x_space());
  Thyra::randomize(-1.0, 1.0, v.ptr());
  auto applyW = [&](const Thyra::LinearOpBase<double>& W) {
    RCP<Vector> Wv = Thyra::createMember(aux->get_f_space());
    Thyra::apply(W, Thyra::NOTRANS, *v, Wv.ptr());
    return Wv;
  };
  auto relDiff = [](const Vector& a, const Vector& b) {
    RCP<Vector> d = a.clone_v();
    Thyra::Vp_StV(d.ptr(), -1.0, b);
    return Thyra::norm_2(*d) / Thyra::norm_2(b);
  };

  // Same alpha and beta: the two W differ in their coupling blocks
  RCP<Thyra::LinearOpBase<double> > W1 =
      evalW(*aux, randomInArgs(*aux, 0.0, 3.0, 0.7));
  RCP<Vector> W1v = applyW(*W1);
  RCP<Thyra::LinearOpBase<double> > W2 =
      evalW(*aux, randomInArgs(*aux, 0.5, 3.0, 0.7));
  RCP<Vector> W2v = applyW(*W2);
  evalResidual(*aux, randomInArgs(*aux, 1.0, 3.0, 0.7));
  out << "|W2*v - W1*v|/|W1*v| = " << relDiff(*W2v, *W1v) << std::endl;
  TEST_COMPARE(relDiff(*W2v, *W1v), >, 1.0e-3);
  TEST_COMPARE(relDiff(*applyW(*W1), *W1v), <=, 1.0e-14);
  TEST_COMPARE(relDiff(*applyW(*W2), *W2v), <=, 1.0e-14);

  // df/dp is still evaluated once per evaluation point
  const int num_dfdp = setup.model->num_dfdp;
  evalW(*aux, randomInArgs(*aux, 1.0, 1.5, 0.4));
  evalResidual(*aux, randomInArgs(*aux, 1.0, 3.0, 0.7));
  TEST_EQUALITY(setup.model->num_dfdp, num_dfdp);
}

// ************************************************************
// The coupling operator s*df/dp^T, with left and right scalings.
TEUCHOS_UNIT_TEST(AdjointAuxCouplingLinearOp, ApplyScalingRowStats)
{
  typedef Tempus::AdjointAuxCouplingLinearOp<double> CouplingOp;
  RCP<const Thyra::VectorSpaceBase<double> > f_space =
      Thyra::defaultSpmdVectorSpace<double>(5);
  RCP<const Thyra::VectorSpaceBase<double> > p_space =
      Thyra::defaultSpmdVectorSpace<double>(3);
  RCP<MultiVector> dfdp = Thyra::createMembers(f_space, 3);
  Thyra::seed_randomize<double>(41);
  Thyra::randomize(-1.0, 1.0, dfdp.ptr());

  RCP<CouplingOp> op = rcp(new CouplingOp(p_space, f_space));
  RCP<Vector> v      = Thyra::createMember(f_space);
  Thyra::randomize(-1.0, 1.0, v.ptr());
  RCP<Vector> y = Thyra::createMember(p_space);

  // Zero before initialize()
  Thyra::assign(y.ptr(), 1.0);
  Thyra::apply(*op, Thyra::NOTRANS, *v, y.ptr());
  TEST_EQUALITY(Thyra::norm_inf(*y), 0.0);

  // s*df/dp^T*v, scaled left and right
  const double s = -0.7;
  op->initialize(dfdp, Thyra::CONJTRANS, s);
  RCP<Vector> dl = Thyra::createMember(p_space);
  RCP<Vector> dr = Thyra::createMember(f_space);
  Thyra::randomize(0.5, 2.0, dl.ptr());
  Thyra::randomize(0.5, 2.0, dr.ptr());
  op->scaleLeft(*dl);
  op->scaleRight(*dr);
  op->scaleLeft(*dl);

  RCP<Vector> drv = v->clone_v();
  Thyra::ele_wise_scale(*dr, drv.ptr());
  RCP<Vector> y_ref = Thyra::createMember(p_space);
  Thyra::apply(*dfdp, Thyra::CONJTRANS, *drv, y_ref.ptr(), s);
  Thyra::ele_wise_scale(*dl, y_ref.ptr());
  Thyra::ele_wise_scale(*dl, y_ref.ptr());
  Thyra::apply(*op, Thyra::NOTRANS, *v, y.ptr());
  Thyra::Vp_StV(y_ref.ptr(), -1.0, *y);
  TEST_COMPARE(Thyra::norm_inf(*y_ref), <=, 1.0e-14 * Thyra::norm_inf(*y));

  // Consistency of the apply and its transpose
  Thyra::LinearOpTester<double> tester;
  TEST_ASSERT(tester.check(*op, Teuchos::inOutArg(out)));

  // The row statistics are zero
  RCP<Vector> rs = Thyra::createMember(p_space);
  op->getRowStat(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_ROW_SUM, rs.ptr());
  TEST_EQUALITY(Thyra::norm_inf(*rs), 0.0);
  TEST_ASSERT(!op->rowStatIsSupported(Thyra::RowStatLinearOpBaseUtils::ROW_STAT_INV_ROW_SUM));

  // initialize() removes the scalings
  op->initialize(dfdp, Thyra::CONJTRANS, s);
  Thyra::apply(*dfdp, Thyra::CONJTRANS, *v, y_ref.ptr(), s);
  Thyra::apply(*op, Thyra::NOTRANS, *v, y.ptr());
  Thyra::Vp_StV(y_ref.ptr(), -1.0, *y);
  TEST_COMPARE(Thyra::norm_inf(*y_ref), <=, 1.0e-14 * Thyra::norm_inf(*y));
}

}  // namespace Tempus_Unit_Test
