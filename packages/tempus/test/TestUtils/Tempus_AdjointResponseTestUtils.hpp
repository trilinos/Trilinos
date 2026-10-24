//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_ADJOINT_RESPONSE_TEST_UTILS_HPP
#define TEMPUS_ADJOINT_RESPONSE_TEST_UTILS_HPP

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorAdjointSensitivity.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "Tempus_ConvergenceTestUtils.hpp"

namespace Tempus_Test {

/** \brief Exact sensitivity dJ/dp (Jacobian form, x_space x num_p) of the
 *  "Time Integral" or "Sum Over Times" response of the SinCos model, whose
 *  response is g(x) = x.
 *
 *  For "Time Integral", J = int_{t0}^{tf} x(t) dt and dJ/dp is computed with
 *  composite 5-point Gauss-Legendre quadrature of the exact sensitivities.
 *  For "Sum Over Times", J = sum_k x(t_k) and dJ/dp = sum_k dx/dp(t_k).
 */
inline Teuchos::RCP<Thyra::MultiVectorBase<double> >
sinCosExactResponseSensitivity(const SinCosModel<double>& model,
                               const std::string& responseType,
                               const double t0, const double tf,
                               const std::vector<double>& responseTimes)
{
  const int num_param = model.get_p_space(0)->dim();
  Teuchos::RCP<Thyra::MultiVectorBase<double> > dJdp =
      Thyra::createMembers(model.get_x_space(), num_param);
  Thyra::assign(dJdp.ptr(), 0.0);

  std::vector<double> times, weights;
  if (responseType == "Time Integral") {
    const double nodes[5]    = {-0.9061798459386640, -0.5384693101056831, 0.0,
                                0.5384693101056831, 0.9061798459386640};
    const double gweights[5] = {0.2369268850561891, 0.4786286704993665,
                                0.5688888888888889, 0.4786286704993665,
                                0.2369268850561891};
    const int numIntervals   = 100;
    const double h           = (tf - t0) / numIntervals;
    for (int k = 0; k < numIntervals; ++k) {
      const double tm = t0 + (k + 0.5) * h;
      for (int q = 0; q < 5; ++q) {
        times.push_back(tm + 0.5 * h * nodes[q]);
        weights.push_back(0.5 * h * gweights[q]);
      }
    }
  }
  else {
    times = responseTimes;
    weights.assign(times.size(), 1.0);
  }

  for (std::size_t k = 0; k < times.size(); ++k)
    for (int i = 0; i < num_param; ++i)
      Thyra::Vp_StV(dJdp->col(i).ptr(), weights[k],
                    *(model.getExactSensSolution(i, times[k]).get_x()));
  return dJdp;
}

/** \brief Solve the SinCos adjoint problem with the "Time Integral" or
 *  "Sum Over Times" response for a sequence of time step sizes, and return
 *  the time step sizes and the errors of dJ/dp.
 *
 *  \param xmlFile        Input file with "SinCosModel" and "Tempus" sublists
 *  \param configure      Called with the "Tempus" sublist to set the stepper,
 *                        interpolator, ...
 *  \param responseType   "Time Integral" or "Sum Over Times"
 *  \param responseTimes  Times of the "Sum Over Times" response
 *  \param dt0            Largest time step size
 *  \param useAdjointModel Use SinCosModelAdjoint as explicit adjoint model
 */
inline void computeSinCosAdjointResponseErrors(
    const std::string& xmlFile,
    const std::function<void(Teuchos::ParameterList&)>& configure,
    const std::string& responseType, const std::vector<double>& responseTimes,
    const int nTimeStepSizes, const double dt0, const bool useAdjointModel,
    std::vector<double>& stepSizes, std::vector<double>& errors,
    double& stepperOrder, Teuchos::FancyOStream& out, bool& success)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::sublist;

  std::ostringstream timeList;
  timeList.precision(17);
  for (std::size_t k = 0; k < responseTimes.size(); ++k)
    timeList << (k > 0 ? ", " : "") << responseTimes[k];

  stepSizes.clear();
  errors.clear();
  double dt = dt0;
  for (int n = 0; n < nTimeStepSizes; ++n, dt /= 2) {
    RCP<ParameterList> pList  = Teuchos::getParametersFromXmlFile(xmlFile);
    RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
    RCP<SinCosModel<double> > model =
        Teuchos::rcp(new SinCosModel<double>(scm_pl));
    RCP<Thyra::ModelEvaluator<double> > adjoint_model;
    if (useAdjointModel)
      adjoint_model = Teuchos::rcp(new SinCosModelAdjoint<double>(scm_pl));

    RCP<ParameterList> pl = sublist(pList, "Tempus", true);
    const std::string integratorName =
        pl->get<std::string>("Integrator Name");
    ParameterList& tscPL =
        pl->sublist(integratorName).sublist("Time Step Control");
    tscPL.set("Initial Time Step", dt);

    ParameterList& sens_pl = pl->sublist("Sensitivities");
    sens_pl.set("Response Type", responseType);
    if (!responseTimes.empty())
      sens_pl.set("Response Time List", timeList.str());

    configure(*pl);

    RCP<Tempus::IntegratorAdjointSensitivity<double> > integrator =
        Tempus::createIntegratorAdjointSensitivity<double>(pl, model,
                                                           adjoint_model);
    stepperOrder = integrator->getStepper()->getOrder();

    // Initial conditions and their sensitivities from the exact solution
    const double t0 = tscPL.get<double>("Initial Time");
    RCP<const Thyra::VectorBase<double> > x0 =
        model->getExactSolution(t0).get_x();
    const int num_param = model->get_p_space(0)->dim();
    RCP<Thyra::MultiVectorBase<double> > DxDp0 =
        Thyra::createMembers(model->get_x_space(), num_param);
    for (int i = 0; i < num_param; ++i)
      Thyra::assign(DxDp0->col(i).ptr(),
                    *(model->getExactSensSolution(i, t0).get_x()));
    integrator->initializeSolutionHistory(t0, x0, Teuchos::null, Teuchos::null,
                                          DxDp0, Teuchos::null, Teuchos::null);

    const bool integratorStatus = integrator->advanceTime();
    TEST_ASSERT(integratorStatus);
    if (!integratorStatus) return;

    const double tf = integrator->getTime();
    TEST_FLOATING_EQUALITY(tf, tscPL.get<double>("Final Time"), 1.0e-12);

    // The integrator returns dJ/dp in gradient form (num_param x num_g)
    RCP<const Thyra::MultiVectorBase<double> > DgDp = integrator->getDgDp();
    RCP<Thyra::MultiVectorBase<double> > dJdp_exact =
        sinCosExactResponseSensitivity(*model, responseType, t0, tf,
                                       responseTimes);
    double error2 = 0.0;
    {
      Thyra::ConstDetachedMultiVectorView<double> dgdp_view(*DgDp);
      Thyra::ConstDetachedMultiVectorView<double> exact_view(*dJdp_exact);
      const int num_g = DgDp->domain()->dim();
      for (int i = 0; i < num_g; ++i)
        for (int j = 0; j < num_param; ++j) {
          const double diff = dgdp_view(j, i) - exact_view(i, j);
          error2 += diff * diff;
        }
    }
    stepSizes.push_back(dt);
    errors.push_back(std::sqrt(error2));
    out << "  dt = " << dt << "  error = " << errors.back() << std::endl;
  }
}

}  // namespace Tempus_Test

#endif  // TEMPUS_ADJOINT_RESPONSE_TEST_UTILS_HPP
