// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_ConfigDefs.hpp"

#if defined(HAVE_PIRO_TEMPUS) && defined(HAVE_PIRO_ROL)

#include "Piro_ThyraProductME_ROL_DynamicObjective.hpp"

#include "MassSpringDamperModel.hpp"

#include "ROL_ThyraVector.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_Core.hpp"

#include <array>
#include <cmath>
#include <vector>

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;

namespace {

typedef Piro::ThyraProductME_ROL_DynamicObjective<double> DynamicObjective;

const double t0 = 0.5, dt = 0.1;
const int nt = 3;

// Time stamps as in Piro::PerformROLTransientAnalysis: time stamp 0 ends at the initial
// time and is not integrated, the time stamps 1, ..., nt cover [t0, t0 + nt*dt]
std::vector<ROL::TimeStamp<double>> timeStamps()
{
  std::vector<ROL::TimeStamp<double>> ts(nt+1);
  for (int k = 0; k <= nt; ++k) {
    ts[k].k = k;
    ts[k].t = {t0 + (k-1)*dt, t0 + k*dt};
  }
  return ts;
}

RCP<Thyra::VectorBase<double>> vector2(const RCP<const Thyra::VectorSpaceBase<double>>& space, double a, double b)
{
  const RCP<Thyra::VectorBase<double>> v = Thyra::createMember(space);
  Thyra::set_ele(0, a, v.ptr());
  Thyra::set_ele(1, b, v.ptr());
  return v;
}

RCP<DynamicObjective> dynamicObjective(const RCP<Thyra::ModelEvaluator<double>>& model, ParameterList& piroParams,
                                       const bool onlyFinalTime)
{
  const RCP<Tempus::Integrator<double>> noIntegrator;
  const RCP<Thyra::ModelEvaluator<double>> noAdjointModel;
  return rcp(new DynamicObjective(model, noIntegrator, noIntegrator, noAdjointModel, 0, piroParams,
                                  nt+1, onlyFinalTime, true, Teuchos::VERB_NONE));
}

// Values of the response g(x, x_dot, p, t) at the two ends of a time stamp, as in the
// documentation of Piro::ThyraProductME_ROL_DynamicObjective
struct EndResponses {
  double g_old, g_new;
};

EndResponses endResponses(const DynamicObjective& obj, const RCP<const Thyra::VectorBase<double>>& uo,
                          const RCP<const Thyra::VectorBase<double>>& un, const RCP<const Thyra::VectorBase<double>>& z,
                          const ROL::TimeStamp<double>& ts)
{
  const double h = ts.t.back() - ts.t.front();
  const RCP<Thyra::VectorBase<double>> u_dot = un->clone_v();
  Thyra::V_StVpStV(u_dot.ptr(), 1.0/h, *un, -1.0/h, *uo);
  EndResponses g;
  g.g_old = obj.value(uo, u_dot, z, ts.t.front());
  g.g_new = obj.value(un, u_dot, z, ts.t.back());
  return g;
}

// The gradients of the contribution of a time stamp with respect to u_old, u_new and z
// are the derivatives of its value (central differences, exact up to rounding for the
// quadratic response of the mass-spring-damper model)
void checkGradients(const DynamicObjective& obj, const ROL::TimeStamp<double>& ts,
                    const RCP<const Thyra::VectorSpaceBase<double>>& x_space,
                    const RCP<const Thyra::VectorSpaceBase<double>>& p_space,
                    Teuchos::FancyOStream& out, bool& success)
{
  const RCP<Thyra::VectorBase<double>> uo = vector2(x_space, 0.3, -0.2);
  const RCP<Thyra::VectorBase<double>> un = vector2(x_space, 0.5, 0.4);
  const RCP<Thyra::VectorBase<double>> z = vector2(p_space, 1.2, 0.8);
  const ROL::ThyraVector<double> rol_uo(uo), rol_un(un), rol_z(z);
  const std::array<RCP<Thyra::VectorBase<double>>, 3> args = {{uo, un, z}};
  const char* names[3] = {"u_old", "u_new", "z"};

  const double h = 1.0e-4;
  for (int i = 0; i < 3; ++i) {
    const RCP<Thyra::VectorBase<double>> grad = args[i]->clone_v();
    ROL::ThyraVector<double> rol_grad(grad);
    if (i == 0) obj.gradient_uo(rol_grad, rol_uo, rol_un, rol_z, ts);
    else if (i == 1) obj.gradient_un(rol_grad, rol_uo, rol_un, rol_z, ts);
    else obj.gradient_z(rol_grad, rol_uo, rol_un, rol_z, ts);

    for (int j = 0; j < 2; ++j) {
      const double a = Thyra::get_ele(*args[i], j);
      Thyra::set_ele(j, a + h, args[i].ptr());
      const double value_plus = obj.value(rol_uo, rol_un, rol_z, ts);
      Thyra::set_ele(j, a - h, args[i].ptr());
      const double value_minus = obj.value(rol_uo, rol_un, rol_z, ts);
      Thyra::set_ele(j, a, args[i].ptr());

      const double fd = (value_plus - value_minus)/(2*h);
      const double derivative = Thyra::get_ele(*grad, j);
      out << "time stamp " << ts.k << ", d/d" << names[i] << "[" << j << "]: " << derivative
          << " (finite difference " << fd << ")\n";
      TEST_COMPARE(std::abs(derivative - fd), <=, 1.0e-7*(1.0 + std::abs(fd)));
    }
  }
}

} // namespace

// The gradients are the derivatives of the values for the three types of responses, with
// and without dependence of the response on x_dot. For "Sum Over Times", the response
// times include the initial time (start of time stamp 1) and a repeated time.
TEUCHOS_UNIT_TEST(Piro_ROL_DynamicObjective, Gradients)
{
  const auto comm = Tpetra::getDefaultComm();
  const std::vector<ROL::TimeStamp<double>> ts = timeStamps();
  const std::vector<double> responseTimes = {t0, t0 + dt, t0 + dt, t0 + 3*dt};

  for (const bool use_x_dot_in_g : {false, true}) {
    const RCP<Thyra::ModelEvaluator<double>> model =
      rcp(new MassSpringDamperModel(comm, false, Teuchos::null, false, use_x_dot_in_g));
    for (const std::string responseType : {"Final Time", "Time Integral", "Sum Over Times"}) {
      out << "Response type \"" << responseType << "\", response depends on x_dot: " << use_x_dot_in_g << "\n";
      ParameterList piroParams("Piro");
      const RCP<DynamicObjective> obj = dynamicObjective(model, piroParams, responseType == "Final Time");
      if (responseType == "Sum Over Times") obj->setResponseTimes(responseTimes, ts);
      for (int k = 1; k <= nt; ++k)
        checkGradients(*obj, ts[k], model->get_x_space(), model->get_p_space(0), out, success);
    }
  }
}

// Values of the time stamp contributions
TEUCHOS_UNIT_TEST(Piro_ROL_DynamicObjective, Values)
{
  const auto comm = Tpetra::getDefaultComm();
  const std::vector<ROL::TimeStamp<double>> ts = timeStamps();
  const RCP<Thyra::ModelEvaluator<double>> model =
    rcp(new MassSpringDamperModel(comm, false, Teuchos::null, false, true));
  const RCP<Thyra::VectorBase<double>> uo = vector2(model->get_x_space(), 0.3, -0.2);
  const RCP<Thyra::VectorBase<double>> un = vector2(model->get_x_space(), 0.5, 0.4);
  const RCP<Thyra::VectorBase<double>> z = vector2(model->get_p_space(0), 1.2, 0.8);
  const ROL::ThyraVector<double> rol_uo(uo), rol_un(un), rol_z(z);
  const double tol = 1.0e-14;

  ParameterList finalParams("Piro"), integralParams("Piro"), sumParams("Piro");
  const RCP<DynamicObjective> finalTime = dynamicObjective(model, finalParams, true);
  const RCP<DynamicObjective> timeIntegral = dynamicObjective(model, integralParams, false);
  const RCP<DynamicObjective> sumOverTimes = dynamicObjective(model, sumParams, true);
  sumOverTimes->setResponseTimes({t0, t0 + dt, t0 + dt, t0 + 3*dt}, ts);

  for (int k = 1; k <= nt; ++k) {
    const EndResponses g = endResponses(*finalTime, uo, un, z, ts[k]);
    out << "time stamp " << k << ": g_old = " << g.g_old << ", g_new = " << g.g_new << "\n";
    TEST_FLOATING_EQUALITY(finalTime->value(rol_uo, rol_un, rol_z, ts[k]), (k == nt ? g.g_new : 0.0), tol);
    TEST_FLOATING_EQUALITY(timeIntegral->value(rol_uo, rol_un, rol_z, ts[k]), dt/2*(g.g_old + g.g_new), tol);
    // Initial time once, t0 + dt twice, t0 + 2 dt not in the list, t0 + 3 dt once
    const double sum = (k == 1) ? g.g_old + 2*g.g_new : ((k == 3) ? g.g_new : 0.0);
    TEST_FLOATING_EQUALITY(sumOverTimes->value(rol_uo, rol_un, rol_z, ts[k]), sum, tol);
  }
}

// The response times must be on the time grid: the initial time or the end of a time stamp
TEUCHOS_UNIT_TEST(Piro_ROL_DynamicObjective, ResponseTimesOnTimeGrid)
{
  const auto comm = Tpetra::getDefaultComm();
  const std::vector<ROL::TimeStamp<double>> ts = timeStamps();
  const RCP<Thyra::ModelEvaluator<double>> model = rcp(new MassSpringDamperModel(comm));
  ParameterList piroParams("Piro");
  const RCP<DynamicObjective> obj = dynamicObjective(model, piroParams, true);

  TEST_NOTHROW(obj->setResponseTimes({t0, t0 + 2*dt, t0 + 3*dt}, ts));
  TEST_THROW(obj->setResponseTimes({t0 + dt/2}, ts), Teuchos::Exceptions::InvalidParameter);
  TEST_THROW(obj->setResponseTimes({t0 - dt}, ts), Teuchos::Exceptions::InvalidParameter);
  TEST_THROW(obj->setResponseTimes({t0 + 4*dt}, ts), Teuchos::Exceptions::InvalidParameter);
  TEST_THROW(obj->setResponseTimes({}, ts), Teuchos::Exceptions::InvalidParameter);
  const std::vector<ROL::TimeStamp<double>> missingTimeStamp(ts.begin(), ts.end() - 1);
  TEST_THROW(obj->setResponseTimes({t0}, missingTimeStamp), Teuchos::Exceptions::InvalidParameter);
}

#endif /* HAVE_PIRO_TEMPUS && HAVE_PIRO_ROL */
