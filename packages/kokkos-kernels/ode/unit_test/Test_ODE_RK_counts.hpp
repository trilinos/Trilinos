// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include "KokkosKernels_TestUtils.hpp"

#include "KokkosODE_RungeKutta.hpp"
#include "Test_ODE_TestProblems.hpp"

namespace Test {

std::string RK_type_to_name(const KokkosODE::Experimental::RK_type RK) {
  std::string name;

  switch (RK) {
    case KokkosODE::Experimental::RK_type::RKFE: name = "Forward-Euler"; break;
    case KokkosODE::Experimental::RK_type::RKEH: name = "Euler-Heun"; break;
    case KokkosODE::Experimental::RK_type::RKF12: name = "Fehlberg 1-2"; break;
    case KokkosODE::Experimental::RK_type::RKBS: name = "Bogacki-Shampine"; break;
    case KokkosODE::Experimental::RK_type::RK4: name = "Classic RK order 4"; break;
    case KokkosODE::Experimental::RK_type::RKF45: name = "Fehlberg 4-5"; break;
    case KokkosODE::Experimental::RK_type::RKCK: name = "Cash-Karp"; break;
    case KokkosODE::Experimental::RK_type::RKDP: name = "Dormand-Prince"; break;
    default: name = "Unknown Runge-Kutta method";
  }

  return name;
}

template <KokkosODE::Experimental::RK_type RK, class Device, class OdeType>
void RK_Count(const Device, const OdeType myODE, const double relTol, const double absTol,
              const int /*expected_count*/) {
  using execution_space = typename Device::execution_space;
  using vec_type        = Kokkos::View<double*, Device>;
  using mv_type         = Kokkos::View<double**, Device>;
  using count_type      = Kokkos::View<int*, execution_space>;

  constexpr int neqs = myODE.neqs;

  constexpr double tstart = myODE.tstart(), tend = myODE.tend();
  constexpr int num_steps      = myODE.numsteps();
  constexpr int maxSteps       = 1e6;
  constexpr double minStepSize = (tend - tstart) / (100 * maxSteps);
  KokkosODE::Experimental::ODE_params params(
      num_steps, maxSteps, 1.0e-12, (RK == KokkosODE::Experimental::RK_type::RKF12) ? 1.0e-8 : 1.0e-6, minStepSize);

  vec_type y("solution", neqs), f("function", neqs);
  vec_type y_new("y new", neqs), y_old("y old", neqs);
  count_type count("time step count", 1);

  auto y_h                                    = Kokkos::create_mirror_view(y);
  typename vec_type::host_mirror_type y_old_h = Kokkos::create_mirror(y_old);
  auto y_ref_h                                = Kokkos::create_mirror(y);
  for (int dofIdx = 0; dofIdx < neqs; ++dofIdx) {
    y_h(dofIdx)     = myODE.expected_val(tstart, dofIdx);
    y_old_h(dofIdx) = y_h(dofIdx);
    y_ref_h(dofIdx) = myODE.expected_val(tend, dofIdx);
  }
  Kokkos::deep_copy(y, y_h);

  vec_type tmp("tmp vector", neqs);
  mv_type kstack("k stack", KokkosODE::Experimental::RungeKutta<RK>::num_stages(), neqs);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  Kokkos::deep_copy(y_old, y_old_h);
  Kokkos::deep_copy(y_new, y_old_h);
  RKSolve_wrapper<OdeType, RK, vec_type, mv_type, double, count_type> solve_wrapper(myODE, params, tstart, tend, y_old,
                                                                                    y_new, tmp, kstack, count);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto y_new_h = Kokkos::create_mirror(y_new);
  Kokkos::deep_copy(y_new_h, y_new);

  typename count_type::host_mirror_type count_h = Kokkos::create_mirror_view(count);
  Kokkos::deep_copy(count_h, count);

  double error = 0.0;
  for (int eqIdx = 0; eqIdx < neqs; ++eqIdx) {
    error += Kokkos::pow(y_ref_h(eqIdx) - y_new_h(eqIdx), 2.0) /
             Kokkos::pow(absTol + relTol * Kokkos::abs(y_new_h(eqIdx)), 2.0);
  }
  error = Kokkos::sqrt(error / neqs);

  std::string msg = std::string(OdeType::name) + ", " + RK_type_to_name(RK);
  EXPECT_LE(error, 1.0) << msg.c_str();
  // EXPECT_LE(count_h(0), expected_count);
}  // RK_Count

}  // namespace Test

template <KokkosODE::Experimental::RK_type RK>
void test_RK_count() {
  //    RK_Count    (Device,       OdeType,                      relTol, absTol, /*expected_count*/)
  Test::RK_Count<RK>(TestDevice(), TestProblem::DegreeOnePoly(), 1.0e-6, 1e-12, 2);
  Test::RK_Count<RK>(TestDevice(), TestProblem::DegreeTwoPoly(), 1.0e-6, 1e-12, 2);
  Test::RK_Count<RK>(TestDevice(), TestProblem::DegreeThreePoly(), 1.0e-6, 1e-12, 2);
  Test::RK_Count<RK>(TestDevice(), TestProblem::DegreeFivePoly(), 1.0e-6, 1e-12, 5);
  Test::RK_Count<RK>(TestDevice(), TestProblem::Exponential(0.7), 2.0e-6, 1e-12, 4);
  Test::RK_Count<RK>(TestDevice(), TestProblem::SpringMassDamper(1001., 1000.), 1.0e-4, 0.0, 272);
  Test::RK_Count<RK>(TestDevice(), TestProblem::CosExp(-10., 2., 1.), 5.3e-5, 0.0, 25);
  if constexpr (RK == KokkosODE::Experimental::RK_type::RKF12) {
    Test::RK_Count<RK>(TestDevice(), TestProblem::StiffChemicalDecayProcess(1e4, 1.), 4e-9, 1e-9, 2786);
  } else {
    Test::RK_Count<RK>(TestDevice(), TestProblem::StiffChemicalDecayProcess(1e4, 1.), 4e-9, 1.8e-10, 2786);
  }
  Test::RK_Count<RK>(TestDevice(), TestProblem::Tracer(10.0), 0.0, 1e-3, 10);
  Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightB5(), 1.3e-2, 0.0, 90);
  if constexpr (RK == KokkosODE::Experimental::RK_type::RKF12) {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightC1(), 1.e-4, 1e-14, 90);
  } else {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightC1(), 1.e-5, 1e-14, 90);
  }
  if constexpr (RK == KokkosODE::Experimental::RK_type::RKF12) {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightC5(), 1.e-4, 1e-14, 97);
  } else {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightC5(), 1.e-5, 1e-14, 97);
  }
  if constexpr (RK == KokkosODE::Experimental::RK_type::RKF12) {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightD2(), 2.e-4, 0.0, 590);
  } else {
    Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightD2(), 1.e-5, 0.0, 590);
  }
  Test::RK_Count<RK>(TestDevice(), TestProblem::EnrightD4(), 1.e-5, 1.e-9, 932);
#if defined(KOKKOS_ENABLE_SYCL)
  if constexpr ((RK != KokkosODE::Experimental::RK_type::RKF12) &&
                !std::is_same_v<typename TestDevice::execution_space, Kokkos::Experimental::SYCL>) {
#else
  if constexpr (RK != KokkosODE::Experimental::RK_type::RKF12) {
#endif
    Test::RK_Count<RK>(TestDevice(), TestProblem::KKStiffChemistry(), 1e-5, 0.0, 1);
  }
}

void test_count() {
  using RK_type = KokkosODE::Experimental::RK_type;

  // test_RK_count<RK_type::RKEH>();
  test_RK_count<RK_type::RKF12>();
  test_RK_count<RK_type::RKBS>();
  test_RK_count<RK_type::RK4>();
  test_RK_count<RK_type::RKF45>();
  test_RK_count<RK_type::RKCK>();
  test_RK_count<RK_type::RKDP>();
  // test_RK_count<RK_type::VER56>();
}

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, RK_Count) { test_count(); }
#endif
