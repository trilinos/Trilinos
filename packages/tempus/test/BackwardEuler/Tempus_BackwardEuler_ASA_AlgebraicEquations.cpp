//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

// Adjoint sensitivities of a semi-explicit index-1 DAE ("Algebraic Equations
// Mask" option of IntegratorAdjointSensitivity), for the "Final Time", "Sum
// Over Times" and "Time Integral" responses.

#include <cmath>
#include <string>
#include <vector>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_DAEStructure.hpp"
#include "Tempus_IntegratorAdjointSensitivity.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

namespace Tempus_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

/** \brief Semi-explicit index-1 DAE with x = (u, h) and parameters p = (p, q):
 *
 *    f_0 = u - kappa h - p          = 0   (algebraic equation),
 *    f_1 = h_dot + h - u - 2 q      = 0   (differential equation),
 *
 *  i.e., u = kappa h + p and h_dot = a h + p + 2 q with a = kappa - 1.  The
 *  response is g = x (two components).  With s(t) = (exp(a t) - 1)/a and
 *  h(0) independent of p:
 *    dh/dp = s,            dh/dq = 2 s,
 *    du/dp = kappa s + 1,  du/dq = 2 kappa s.
 *  The term 1 of du/dp comes from the algebraic component of the adjoint jump
 *  of the "Final Time" and "Sum Over Times" responses.
 */
class AlgebraicDAEModel : public Thyra::StateFuncModelEvaluatorBase<double> {
 public:
  typedef Thyra::ModelEvaluatorBase MEB;

  AlgebraicDAEModel(const double kappa = 0.5, const double h0 = 1.0,
                    const double p0 = 0.3, const double q0 = 0.1)
    : kappa_(kappa), h0_(h0), p0_(p0), q0_(q0)
  {
    space_ = Thyra::defaultSpmdVectorSpace<double>(2);
    p_names_.push_back("p");
    p_names_.push_back("q");
    g_names_.push_back("g = x");

    // The descriptions must be this->description(), which is checked by
    // Thyra::ModelEvaluatorDefaultBase in debug builds
    MEB::InArgsSetup<double> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_t);
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs.setSupports(MEB::IN_ARG_x_dot);
    inArgs.setSupports(MEB::IN_ARG_alpha);
    inArgs.setSupports(MEB::IN_ARG_beta);
    inArgs.set_Np(1);
    inArgs_ = inArgs;

    MEB::OutArgsSetup<double> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    outArgs.set_Np_Ng(1, 1);
    outArgs.setSupports(MEB::OUT_ARG_DfDp, 0, MEB::DERIV_MV_JACOBIAN_FORM);
    outArgs.setSupports(MEB::OUT_ARG_DgDx, 0, MEB::DERIV_MV_GRADIENT_FORM);
    outArgs.setSupports(MEB::OUT_ARG_DgDp, 0, 0, MEB::DERIV_MV_JACOBIAN_FORM);
    outArgs_ = outArgs;

    nominal_ = inArgs_;
    nominal_.set_t(0.0);
    RCP<Thyra::VectorBase<double> > x     = Thyra::createMember(space_);
    RCP<Thyra::VectorBase<double> > x_dot = Thyra::createMember(space_);
    RCP<Thyra::VectorBase<double> > p     = Thyra::createMember(space_);
    {
      Thyra::DetachedVectorView<double> xv(*x), xdv(*x_dot), pv(*p);
      xv[0]  = kappa_ * h0_ + p0_;
      xv[1]  = h0_;
      xdv[0] = 0.0;
      xdv[1] = a() * h0_ + p0_ + 2.0 * q0_;
      pv[0]  = p0_;
      pv[1]  = q0_;
    }
    nominal_.set_x(x);
    nominal_.set_x_dot(x_dot);
    nominal_.set_p(0, p);
  }

  double a() const { return kappa_ - 1.0; }
  double kappa() const { return kappa_; }

  /// Consistent initial condition x(0) = (kappa h0 + p, h0)
  RCP<const Thyra::VectorBase<double> > getInitialCondition() const
  {
    return nominal_.get_x();
  }

  /// dx(0)/dp (2 x 2): u(0) depends on p
  RCP<Thyra::MultiVectorBase<double> > getInitialSensitivity() const
  {
    RCP<Thyra::MultiVectorBase<double> > dxdp = Thyra::createMembers(space_, 2);
    Thyra::assign(dxdp.ptr(), 0.0);
    Thyra::DetachedMultiVectorView<double> v(*dxdp);
    v(0, 0) = 1.0;
    return dxdp;
  }

  /// Exact dx/dp(t), rows: (u, h), columns: (p, q)
  void exactSensitivity(const double t, double s[2][2]) const
  {
    const double st = (std::exp(a() * t) - 1.0) / a();
    s[0][0]         = kappa_ * st + 1.0;
    s[0][1]         = 2.0 * kappa_ * st;
    s[1][0]         = st;
    s[1][1]         = 2.0 * st;
  }

  /// Exact int_0^t dx/dp
  void exactIntegratedSensitivity(const double t, double s[2][2]) const
  {
    const double St = ((std::exp(a() * t) - 1.0) / a() - t) / a();
    s[0][0]         = kappa_ * St + t;
    s[0][1]         = 2.0 * kappa_ * St;
    s[1][0]         = St;
    s[1][1]         = 2.0 * St;
  }

  RCP<const Thyra::VectorSpaceBase<double> > get_x_space() const
  {
    return space_;
  }
  RCP<const Thyra::VectorSpaceBase<double> > get_f_space() const
  {
    return space_;
  }
  RCP<const Thyra::VectorSpaceBase<double> > get_p_space(int l) const
  {
    TEUCHOS_ASSERT_EQUALITY(l, 0);
    return space_;
  }
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const
  {
    TEUCHOS_ASSERT_EQUALITY(l, 0);
    return Teuchos::rcpFromRef(p_names_);
  }
  RCP<const Thyra::VectorSpaceBase<double> > get_g_space(int j) const
  {
    TEUCHOS_ASSERT_EQUALITY(j, 0);
    return space_;
  }
  Teuchos::ArrayView<const std::string> get_g_names(int j) const
  {
    TEUCHOS_ASSERT_EQUALITY(j, 0);
    return g_names_();
  }
  MEB::InArgs<double> getNominalValues() const { return nominal_; }
  MEB::InArgs<double> createInArgs() const { return inArgs_; }
  RCP<Thyra::LinearOpBase<double> > create_W_op() const
  {
    return Thyra::createMembers(space_, 2);
  }
  RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > get_W_factory()
      const
  {
    return Thyra::defaultSerialDenseLinearOpWithSolveFactory<double>();
  }

 private:
  MEB::OutArgs<double> createOutArgsImpl() const { return outArgs_; }

  void evalModelImpl(const MEB::InArgs<double>& inArgs,
                     const MEB::OutArgs<double>& outArgs) const
  {
    Thyra::ConstDetachedVectorView<double> x(*inArgs.get_x().assert_not_null());
    RCP<const Thyra::VectorBase<double> > p_in = inArgs.get_p(0);
    if (p_in == Teuchos::null) p_in = nominal_.get_p(0);
    Thyra::ConstDetachedVectorView<double> p(*p_in);
    double h_dot = 0.0;
    if (inArgs.get_x_dot() != Teuchos::null)
      h_dot = Thyra::get_ele(*inArgs.get_x_dot(), 1);

    RCP<Thyra::VectorBase<double> > f = outArgs.get_f();
    if (f != Teuchos::null) {
      Thyra::DetachedVectorView<double> fv(*f);
      fv[0] = x[0] - kappa_ * x[1] - p[0];
      fv[1] = h_dot + x[1] - x[0] - 2.0 * p[1];
    }

    // W = alpha df/dx_dot + beta df/dx, with df/dx_dot = [0 0; 0 1]
    RCP<Thyra::LinearOpBase<double> > W_op = outArgs.get_W_op();
    if (W_op != Teuchos::null) {
      const double alpha = inArgs.get_alpha();
      const double beta  = inArgs.get_beta();
      RCP<Thyra::MultiVectorBase<double> > W =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<double> >(W_op, true);
      Thyra::DetachedMultiVectorView<double> Wv(*W);
      Wv(0, 0) = beta;
      Wv(0, 1) = -kappa_ * beta;
      Wv(1, 0) = -beta;
      Wv(1, 1) = alpha + beta;
    }

    RCP<Thyra::MultiVectorBase<double> > dfdp = outArgs.get_DfDp(0).getMultiVector();
    if (dfdp != Teuchos::null) {
      Thyra::DetachedMultiVectorView<double> v(*dfdp);
      v(0, 0) = -1.0;
      v(0, 1) = 0.0;
      v(1, 0) = 0.0;
      v(1, 1) = -2.0;
    }

    RCP<Thyra::VectorBase<double> > g = outArgs.get_g(0);
    if (g != Teuchos::null) Thyra::assign(g.ptr(), *inArgs.get_x());

    RCP<Thyra::MultiVectorBase<double> > dgdx = outArgs.get_DgDx(0).getMultiVector();
    if (dgdx != Teuchos::null) {
      Thyra::DetachedMultiVectorView<double> v(*dgdx);
      v(0, 0) = 1.0;
      v(0, 1) = 0.0;
      v(1, 0) = 0.0;
      v(1, 1) = 1.0;
    }

    RCP<Thyra::MultiVectorBase<double> > dgdp = outArgs.get_DgDp(0, 0).getMultiVector();
    if (dgdp != Teuchos::null) Thyra::assign(dgdp.ptr(), 0.0);
  }

  double kappa_, h0_, p0_, q0_;
  RCP<const Thyra::VectorSpaceBase<double> > space_;
  Teuchos::Array<std::string> p_names_;
  Teuchos::Array<std::string> g_names_;
  MEB::InArgs<double> inArgs_;
  MEB::OutArgs<double> outArgs_;
  MEB::InArgs<double> nominal_;
};

/// Tempus parameters: Backward Euler with a constant time step dt on [0, 1]
RCP<ParameterList> getTempusParameters(const double dt)
{
  RCP<ParameterList> pl = Teuchos::getParametersFromXmlString(
      "<ParameterList name=\"Tempus\">"
      "  <Parameter name=\"Integrator Name\" type=\"string\""
      "             value=\"Default Integrator\"/>"
      "  <ParameterList name=\"Default Integrator\">"
      "    <Parameter name=\"Integrator Type\" type=\"string\""
      "               value=\"Integrator Basic\"/>"
      "    <Parameter name=\"Screen Output Index Interval\" type=\"int\""
      "               value=\"100000\"/>"
      "    <Parameter name=\"Stepper Name\" type=\"string\""
      "               value=\"Default Stepper\"/>"
      "    <ParameterList name=\"Solution History\">"
      "      <Parameter name=\"Storage Type\" type=\"string\""
      "                 value=\"Unlimited\"/>"
      "      <Parameter name=\"Storage Limit\" type=\"int\" value=\"20\"/>"
      "      <ParameterList name=\"Interpolator\">"
      "        <Parameter name=\"Interpolator Type\" type=\"string\""
      "                   value=\"Lagrange\"/>"
      "        <Parameter name=\"Order\" type=\"int\" value=\"0\"/>"
      "      </ParameterList>"
      "    </ParameterList>"
      "    <ParameterList name=\"Time Step Control\">"
      "      <Parameter name=\"Initial Time\" type=\"double\" value=\"0.0\"/>"
      "      <Parameter name=\"Final Time\" type=\"double\" value=\"1.0\"/>"
      "      <Parameter name=\"Final Time Index\" type=\"int\""
      "                 value=\"100000\"/>"
      "      <Parameter name=\"Minimum Time Step\" type=\"double\""
      "                 value=\"1.0e-6\"/>"
      "      <Parameter name=\"Maximum Time Step\" type=\"double\""
      "                 value=\"1.0\"/>"
      "    </ParameterList>"
      "  </ParameterList>"
      "  <ParameterList name=\"Default Stepper\">"
      "    <Parameter name=\"Stepper Type\" type=\"string\""
      "               value=\"Backward Euler\"/>"
      "    <Parameter name=\"Use FSAL\" type=\"bool\" value=\"false\"/>"
      "    <Parameter name=\"Initial Condition Consistency\" type=\"string\""
      "               value=\"None\"/>"
      "    <Parameter name=\"Initial Condition Consistency Check\""
      "               type=\"bool\" value=\"false\"/>"
      "    <Parameter name=\"Solver Name\" type=\"string\""
      "               value=\"Default Solver\"/>"
      "    <Parameter name=\"Zero Initial Guess\" type=\"bool\" value=\"0\"/>"
      "    <ParameterList name=\"Default Solver\">"
      "      <ParameterList name=\"NOX\">"
      "        <ParameterList name=\"Direction\">"
      "          <Parameter name=\"Method\" type=\"string\" value=\"Newton\"/>"
      "        </ParameterList>"
      "        <ParameterList name=\"Line Search\">"
      "          <Parameter name=\"Method\" type=\"string\""
      "                     value=\"Full Step\"/>"
      "        </ParameterList>"
      "        <Parameter name=\"Nonlinear Solver\" type=\"string\""
      "                   value=\"Line Search Based\"/>"
      "        <ParameterList name=\"Status Tests\">"
      "          <Parameter name=\"Test Type\" type=\"string\" value=\"Combo\"/>"
      "          <Parameter name=\"Combo Type\" type=\"string\" value=\"OR\"/>"
      "          <Parameter name=\"Number of Tests\" type=\"int\" value=\"2\"/>"
      "          <ParameterList name=\"Test 0\">"
      "            <Parameter name=\"Test Type\" type=\"string\""
      "                       value=\"NormF\"/>"
      "            <Parameter name=\"Tolerance\" type=\"double\""
      "                       value=\"1.0e-12\"/>"
      "          </ParameterList>"
      "          <ParameterList name=\"Test 1\">"
      "            <Parameter name=\"Test Type\" type=\"string\""
      "                       value=\"MaxIters\"/>"
      "            <Parameter name=\"Maximum Iterations\" type=\"int\""
      "                       value=\"10\"/>"
      "          </ParameterList>"
      "        </ParameterList>"
      "      </ParameterList>"
      "    </ParameterList>"
      "  </ParameterList>"
      "</ParameterList>");
  pl->sublist("Default Integrator")
      .sublist("Time Step Control")
      .set("Initial Time Step", dt);
  return pl;
}

/// Mask of the algebraic equations
RCP<const Thyra::VectorBase<double> > createMask(
    const AlgebraicDAEModel& model, const double m0, const double m1)
{
  RCP<Thyra::VectorBase<double> > mask = Thyra::createMember(model.get_f_space());
  Thyra::DetachedVectorView<double> v(*mask);
  v[0] = m0;
  v[1] = m1;
  return mask;
}

/// Run the adjoint sensitivity analysis; returns dg/dp (Jacobian form)
RCP<Thyra::MultiVectorBase<double> > runAdjoint(
    const RCP<AlgebraicDAEModel>& model, const std::string& responseType,
    const std::string& responseTimes, const double dt,
    const RCP<const Thyra::VectorBase<double> >& mask, bool& status)
{
  RCP<ParameterList> pl  = getTempusParameters(dt);
  ParameterList& sens_pl = pl->sublist("Sensitivities");
  sens_pl.set("Mass Matrix Is Identity", false);
  sens_pl.set("Response Type", responseType);
  if (!responseTimes.empty()) sens_pl.set("Response Time List", responseTimes);
  sens_pl.set(Tempus::DAEStructure<double>::algebraicMaskName(), mask);
  sens_pl.sublist("Algebraic Jump Solver").set("Relative Tolerance", 1.0e-12);

  RCP<Tempus::IntegratorAdjointSensitivity<double> > integrator =
      Tempus::createIntegratorAdjointSensitivity<double>(pl, model);
  integrator->initializeSolutionHistory(
      0.0, model->getInitialCondition(), Teuchos::null, Teuchos::null,
      model->getInitialSensitivity(), Teuchos::null, Teuchos::null);
  status = integrator->advanceTime();

  // getDgDp() is in gradient form: (p, g)
  RCP<const Thyra::MultiVectorBase<double> > DgDp = integrator->getDgDp();
  RCP<Thyra::MultiVectorBase<double> > dgdp       = Thyra::createMembers(model->get_g_space(0), 2);
  Thyra::ConstDetachedMultiVectorView<double> grad(*DgDp);
  Thyra::DetachedMultiVectorView<double> jac(*dgdp);
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j) jac(i, j) = grad(j, i);
  return dgdp;
}

/// Convergence of the adjoint sensitivities to the exact ones (order 1)
void checkConvergence(const std::string& responseType,
                      Teuchos::FancyOStream& out, bool& success)
{
  RCP<AlgebraicDAEModel> model    = rcp(new AlgebraicDAEModel());
  const std::vector<double> times = {0.25, 0.5, 1.0};
  const std::string timeList      = "0.25, 0.5, 1.0";

  // Exact sensitivity of the response
  double exact[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  if (responseType == "Final Time")
    model->exactSensitivity(1.0, exact);
  else if (responseType == "Time Integral")
    model->exactIntegratedSensitivity(1.0, exact);
  else {
    for (const double t : times) {
      double s[2][2];
      model->exactSensitivity(t, s);
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) exact[i][j] += s[i][j];
    }
  }

  std::vector<double> stepSize, errorNorm;
  for (int n = 0; n < 5; ++n) {
    const double dt = 1.0 / (16 << n);  // exact in binary
    bool status     = false;
    RCP<Thyra::MultiVectorBase<double> > dgdp =
        runAdjoint(model, responseType,
                   responseType == "Sum Over Times" ? timeList : "", dt,
                   createMask(*model, 1.0, 0.0), status);
    TEST_ASSERT(status);

    // Backward Euler sensitivity s_k of h at t_k = k dt (exact for the
    // discrete solution, for comparison only)
    double s_BE = 0.0, sum_BE = 0.0;
    const int N = 16 << n;
    for (int k = 1; k <= N; ++k) {
      s_BE = (s_BE + dt) / (1.0 - model->a() * dt);
      for (const double t : times)
        if (std::abs(k * dt - t) < 1.0e-12) sum_BE += s_BE;
    }
    const double s_disc = responseType == "Final Time" ? s_BE : sum_BE;
    const double ns     = responseType == "Final Time" ? 1.0 : 3.0;

    Thyra::ConstDetachedMultiVectorView<double> v(*dgdp);
    double err2 = 0.0;
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        err2 += (v(i, j) - exact[i][j]) * (v(i, j) - exact[i][j]);
    stepSize.push_back(dt);
    errorNorm.push_back(std::sqrt(err2));
    out << "  " << responseType << ", dt = " << dt << ": du/dp = " << v(0, 0)
        << ", du/dq = " << v(0, 1) << ", dh/dp = " << v(1, 0)
        << ", dh/dq = " << v(1, 1) << ", |error| = " << std::sqrt(err2);
    if (responseType != "Time Integral")
      out << ", du/dp - (discrete BE) = "
          << v(0, 0) - (model->kappa() * s_disc + ns);
    out << std::endl;
  }

  const double slope = computeLinearRegressionLogLog<double>(stepSize, errorNorm);
  out << "  " << responseType << ": order = " << slope << std::endl;
  TEST_FLOATING_EQUALITY(slope, 1.0, 0.05);
  TEST_COMPARE(errorNorm.back(), <, 1.0e-2);
}

// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_FinalTime)
{
  checkConvergence("Final Time", out, success);
}

// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_SumOverTimes)
{
  checkConvergence("Sum Over Times", out, success);
}

// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_TimeIntegral)
{
  checkConvergence("Time Integral", out, success);
}

// ************************************************************
// The algebraic component of the jump: for g = u at the final time, the
// sensitivity to p jumps by 1 (the ODE formula would need df/dx_dot^{-T},
// which does not exist).
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_AlgebraicJump)
{
  RCP<AlgebraicDAEModel> model              = rcp(new AlgebraicDAEModel());
  bool status                               = false;
  const double dt                           = 1.0 / 64;
  RCP<Thyra::MultiVectorBase<double> > dgdp = runAdjoint(
      model, "Final Time", "", dt, createMask(*model, 1.0, 0.0), status);
  TEST_ASSERT(status);
  Thyra::ConstDetachedMultiVectorView<double> v(*dgdp);
  // du/dp = kappa dh/dp + 1 and du/dq = kappa dh/dq hold exactly for the
  // discrete solution
  TEST_FLOATING_EQUALITY(v(0, 0), model->kappa() * v(1, 0) + 1.0, 1.0e-10);
  TEST_FLOATING_EQUALITY(v(0, 1), model->kappa() * v(1, 1), 1.0e-10);
}

// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_InvalidMask)
{
  RCP<AlgebraicDAEModel> model = rcp(new AlgebraicDAEModel());
  bool status                  = false;

  // Entries must be 0 or 1
  TEST_THROW(runAdjoint(model, "Final Time", "", 0.1,
                        createMask(*model, 0.5, 0.0), status),
             std::invalid_argument);

  // The differential equation depends on x_dot and cannot be flagged
  bool massMatrixError = false;
  try {
    runAdjoint(model, "Final Time", "", 0.1, createMask(*model, 1.0, 1.0),
               status);
  }
  catch (const std::logic_error& e) {
    massMatrixError = std::string(e.what()).find(
                          "the rows of the mass matrix df/dx_dot") !=
                      std::string::npos;
    if (!massMatrixError) out << "Unexpected exception: " << e.what() << std::endl;
  }
  TEST_ASSERT(massMatrixError);
}

// ************************************************************
TEUCHOS_UNIT_TEST(BackwardEuler, ASA_AlgebraicEquations_DAEStructure)
{
  RCP<AlgebraicDAEModel> model = rcp(new AlgebraicDAEModel());
  Tempus::DAEStructure<double> dae(createMask(*model, 1.0, 0.0));
  TEST_EQUALITY(Thyra::get_ele(*dae.getDifferentialMask(), 0), 0.0);
  TEST_EQUALITY(Thyra::get_ele(*dae.getDifferentialMask(), 1), 1.0);

  RCP<Thyra::MultiVectorBase<double> > v = Thyra::createMembers(model->get_f_space(), 2);
  Thyra::assign(v.ptr(), 3.0);
  dae.zeroAlgebraicComponents(*v);
  {
    Thyra::ConstDetachedMultiVectorView<double> w(*v);
    TEST_EQUALITY(w(0, 0), 0.0);
    TEST_EQUALITY(w(0, 1), 0.0);
    TEST_EQUALITY(w(1, 0), 3.0);
    TEST_EQUALITY(w(1, 1), 3.0);
  }
  Thyra::assign(v.ptr(), 3.0);
  dae.zeroDifferentialComponents(*v);
  {
    Thyra::ConstDetachedMultiVectorView<double> w(*v);
    TEST_EQUALITY(w(0, 0), 3.0);
    TEST_EQUALITY(w(1, 1), 0.0);
  }

  // Absent or null mask: no DAE structure
  ParameterList pl;
  TEST_ASSERT(Tempus::getDAEStructure<double>(pl) == Teuchos::null);
  pl.set(Tempus::DAEStructure<double>::algebraicMaskName(),
         createMask(*model, 1.0, 0.0));
  TEST_ASSERT(Tempus::getDAEStructure<double>(pl) != Teuchos::null);
}

}  // namespace Tempus_Test
