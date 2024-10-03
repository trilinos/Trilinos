//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <gtest/gtest.h>
#include "KokkosKernels_TestUtils.hpp"

#include "KokkosODE_RungeKutta.hpp"

namespace Test {

// damped harmonic undriven oscillator
// m y'' + c y' + k y = 0
// solution: y=A * exp(-xi * omega_0 * t) * sin(sqrt(1-xi^2) * omega_0 * t +
// phi) omega_0 = sqrt(k/m); xi = c / sqrt(4*m*k) A and phi depend on y(0) and
// y'(0); Change of variables: x(t) = y(t)*exp(-c/(2m)*t) = y(t)*exp(-xi *
// omega_0 * t) Change of variables: X = [x ]
//                          [x']
// Leads to X' = A*X  with A = [ 0  1]
//                             [-d  0]
// with d = k/m - (c/(2m)^2) = (1 - xi^2)*omega_0^2
struct duho {
  constexpr static int neqs = 2;
  const double m, c, k, d;
  const double a11 = 0, a12 = 1, a21, a22;

  duho(const double m_, const double c_, const double k_)
      : m(m_), c(c_), k(k_), d(k_ / m_ - (c_ * c_) / (4 * m_ * m_)), a21(-k / m), a22(-c / m){};

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    f(0) = a11 * y(0) + a12 * y(1);
    f(1) = a21 * y(0) + a22 * y(1);
  }

  template <class vec_type>
  KOKKOS_FUNCTION void solution(const double t, const vec_type& y0, const vec_type& y) const {
    using KAT = Kokkos::ArithTraits<double>;

    const double gamma = c / (2 * m);
    const double omega = KAT::sqrt(k / m - gamma * gamma);
    const double phi   = KAT::atan((y0(1) + gamma * y0(0)) / (y0(0) * omega));
    const double A     = y0(0) / KAT::cos(phi);

    y(0) = A * KAT::cos(omega * t - phi) * KAT::exp(-t * gamma);
    y(1) = -y(0) * gamma - omega * A * KAT::sin(omega * t - phi) * KAT::exp(-t * gamma);
  }

};  // duho

template <class ode_type, class vec_type, class scalar_type>
struct solution_wrapper {
  ode_type ode;
  scalar_type t;
  vec_type y_old, y_ref;

  solution_wrapper(const ode_type& ode_, const scalar_type t_, const vec_type& y_old_, const vec_type& y_ref_)
      : ode(ode_), t(t_), y_old(y_old_), y_ref(y_ref_){};

  KOKKOS_FUNCTION
  void operator()(const int /*idx*/) const { ode.solution(t, y_old, y_ref); }
};

template <class ode_type, KokkosODE::Experimental::RK_type rk_type, class vec_type, class mv_type, class scalar_type>
struct RKSolve_wrapper {
  using ode_params = KokkosODE::Experimental::ODE_params;

  ode_type my_ode;
  ode_params params;
  scalar_type tstart, tend;
  int max_steps;
  vec_type y_old, y_new, tmp;
  mv_type kstack;

  RKSolve_wrapper(const ode_type& my_ode_, const ode_params& params_, const scalar_type tstart_,
                  const scalar_type tend_, const vec_type& y_old_, const vec_type& y_new_, const vec_type& tmp_,
                  const mv_type& kstack_)
      : my_ode(my_ode_),
        params(params_),
        tstart(tstart_),
        tend(tend_),
        y_old(y_old_),
        y_new(y_new_),
        tmp(tmp_),
        kstack(kstack_) {}

  KOKKOS_FUNCTION
  void operator()(const int /*idx*/) const {
    KokkosODE::Experimental::RungeKutta<rk_type>::Solve(my_ode, params, tstart, tend, y_old, y_new, tmp, kstack);
  }
};

template <class ode_type, KokkosODE::Experimental::RK_type rk_type, class vec_type, class mv_type, class scalar_type>
void test_method(const std::string label, ode_type& my_ode, const scalar_type& tstart, const scalar_type& tend,
                 const int num_steps, vec_type& y_old, vec_type& y_new, const int order, const int num_stages,
                 const Kokkos::View<double**, Kokkos::HostSpace>& ks,
                 const Kokkos::View<double*, Kokkos::HostSpace>& sol, typename vec_type::HostMirror y_ref_h) {
  using execution_space = typename vec_type::execution_space;
  using solver_type     = KokkosODE::Experimental::RungeKutta<rk_type>;

  KokkosODE::Experimental::ODE_params params(num_steps);
  vec_type tmp("tmp vector", my_ode.neqs);
  mv_type kstack("k stack", solver_type::num_stages(), my_ode.neqs);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  RKSolve_wrapper<ode_type, rk_type, vec_type, mv_type, scalar_type> solve_wrapper(my_ode, params, tstart, tend, y_old,
                                                                                   y_new, tmp, kstack);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto y_new_h = Kokkos::create_mirror_view(y_new);
  Kokkos::deep_copy(y_new_h, y_new);
  auto kstack_h = Kokkos::create_mirror_view(kstack);
  Kokkos::deep_copy(kstack_h, kstack);

  EXPECT_EQ(solver_type::order(), order);
  EXPECT_EQ(solver_type::num_stages(), num_stages);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\n" << label << std::endl;
  std::cout << "  order: " << solver_type::order() << std::endl;
  std::cout << "  number of stages: " << solver_type::num_stages() << std::endl;
#else
  (void)label;
#endif
  for (int stageIdx = 0; stageIdx < solver_type::num_stages(); ++stageIdx) {
    EXPECT_NEAR_KK(ks(0, stageIdx), kstack_h(stageIdx, 0), 1e-8);
    EXPECT_NEAR_KK(ks(1, stageIdx), kstack_h(stageIdx, 1), 1e-8);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    std::cout << "  k" << stageIdx << "={" << kstack_h(stageIdx, 0) << ", " << kstack_h(stageIdx, 1) << "}"
              << std::endl;
#endif
  }
  EXPECT_NEAR_KK(sol(0), y_new_h(0), 1e-8);
  EXPECT_NEAR_KK(sol(1), y_new_h(1), 1e-8);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "  y={" << y_new_h(0) << ", " << y_new_h(1) << "}" << std::endl;
  std::cout << "  error={" << Kokkos::abs(y_new_h(0) - y_ref_h(0)) / Kokkos::abs(y_ref_h(0)) << ", "
            << Kokkos::abs(y_new_h(1) - y_ref_h(1)) / Kokkos::abs(y_ref_h(1)) << "}" << std::endl;
#else
  (void)y_ref_h;
#endif

}  // test_method

template <class Device>
void test_RK() {
  using execution_space = typename Device::execution_space;
  using RK_type         = KokkosODE::Experimental::RK_type;
  using vec_type        = Kokkos::View<double*, Device>;
  using mv_type         = Kokkos::View<double**, Device>;

  duho my_oscillator(1, 1, 4);
  const int neqs = my_oscillator.neqs;

  vec_type y("solution", neqs), f("function", neqs);
  auto y_h = Kokkos::create_mirror(y);
  y_h(0)   = 1;
  y_h(1)   = 0;
  Kokkos::deep_copy(y, y_h);

  constexpr double tstart = 0, tend = 0.01;
  constexpr int num_steps = 1000;
  double dt               = (tend - tstart) / num_steps;
  vec_type y_new("y new", neqs), y_old("y old", neqs);

  // Since y_old_h will be reused to set initial conditions
  // for each method tested we do not want to use
  // create_mirror_view which would not do a copy
  // when y_old is in HostSpace.
  typename vec_type::HostMirror y_old_h = Kokkos::create_mirror(y_old);
  y_old_h(0)                            = 1;
  y_old_h(1)                            = 0;

  // First compute analytical solution as reference
  // and to evaluate the error from each RK method.
  vec_type y_ref("reference value", neqs);
  auto y_ref_h = Kokkos::create_mirror(y_ref);
  {
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    solution_wrapper wrapper(my_oscillator, tstart + dt, y_old, y_ref);
    Kokkos::parallel_for(my_policy, wrapper);

    Kokkos::deep_copy(y_ref_h, y_ref);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    std::cout << "\nAnalytical solution" << std::endl;
    std::cout << "  y={" << y_ref_h(0) << ", " << y_ref_h(1) << "}" << std::endl;
#endif
  }

  // We perform a single step using a RK method
  // and check the values for ki and y_new against
  // expected values.
  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[2] = {0, -4};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 1);
    double sol_raw[2] = {1, -0.04};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKFE, vec_type, mv_type, double>("Euler-Forward", my_oscillator, tstart, tend, 1, y_old,
                                                                y_new, 1, 1, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[4] = {0, -0.04, -4, -3.96};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 2);
    double sol_raw[2] = {0.9998, -0.0398};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKEH, vec_type, mv_type, double>("Euler-Heun", my_oscillator, tstart, tend, 1, y_old,
                                                                y_new, 2, 2, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[6] = {0, -0.02, -0.03980078, -4, -3.98, -3.95940234};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 3);
    double sol_raw[2] = {0.9998, -0.03979999};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKF12, vec_type, mv_type, double>("RKF-12", my_oscillator, tstart, tend, 1, y_old, y_new,
                                                                 2, 3, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[8] = {0, -0.02, -0.02985, -0.039798, -4, -3.98, -3.96955, -3.95940467};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 4);
    double sol_raw[2] = {0.99980067, -0.039798};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKBS, vec_type, mv_type, double>("RKBS", my_oscillator, tstart, tend, 1, y_old, y_new, 3,
                                                                4, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[12] = {0,  -0.01, -0.01497188, -0.03674986, -0.03979499, -0.0199505,
                         -4, -3.99, -3.98491562, -3.96257222, -3.95941166, -3.97984883};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 6);
    double sol_raw[2] = {0.99980067, -0.03979801};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKF45, vec_type, mv_type, double>("RKF-45", my_oscillator, tstart, tend, 1, y_old, y_new,
                                                                 5, 6, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[12] = {0,  -0.008, -0.011982, -0.02392735, -0.03979862, -0.03484563,
                         -4, -3.992, -3.987946, -3.97578551, -3.95940328, -3.96454357};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 6);
    double sol_raw[2] = {0.99980067, -0.03979801};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKCK, vec_type, mv_type, double>("Cash-Karp", my_oscillator, tstart, tend, 1, y_old,
                                                                y_new, 5, 6, ks, sol, y_ref_h);
  }

  {
    Kokkos::deep_copy(y_old, y_old_h);
    double ks_raw[14] = {0,  -0.008, -0.011982, -0.03187008, -0.03539333, -0.0397954,  -0.03979801,
                         -4, -3.992, -3.987946, -3.96762048, -3.96398013, -3.95941068, -3.95940467};
    Kokkos::View<double**, Kokkos::HostSpace> ks(ks_raw, 2, 7);
    double sol_raw[2] = {0.99980067, -0.03979801};
    Kokkos::View<double*, Kokkos::HostSpace> sol(sol_raw, 2);
    test_method<duho, RK_type::RKDP, vec_type, mv_type, double>("Dormand-Prince", my_oscillator, tstart, tend, 1, y_old,
                                                                y_new, 5, 7, ks, sol, y_ref_h);
  }

}  // test_RK

template <class ode_type, KokkosODE::Experimental::RK_type rk_type, class vec_type, class mv_type, class scalar_type>
void test_rate(ode_type& my_ode, const scalar_type& tstart, const scalar_type& tend,
               Kokkos::View<int*, Kokkos::HostSpace> num_steps, typename vec_type::HostMirror& y_old_h,
               typename vec_type::HostMirror& y_ref_h, typename vec_type::HostMirror& error) {
  using execution_space = typename vec_type::execution_space;
  using solver_type     = KokkosODE::Experimental::RungeKutta<rk_type>;

  vec_type tmp("tmp vector", my_ode.neqs);
  mv_type kstack("k stack", solver_type::num_stages(), my_ode.neqs);

  vec_type y_new("solution", my_ode.neqs);
  vec_type y_old("intial conditions", my_ode.neqs);
  auto y_new_h = Kokkos::create_mirror(y_new);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  for (int idx = 0; idx < num_steps.extent_int(0); ++idx) {
    KokkosODE::Experimental::ODE_params params(num_steps(idx));
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::deep_copy(y_new, y_old_h);
    RKSolve_wrapper<ode_type, rk_type, vec_type, mv_type, scalar_type> solve_wrapper(my_ode, params, tstart, tend,
                                                                                     y_old, y_new, tmp, kstack);
    Kokkos::parallel_for(my_policy, solve_wrapper);

    Kokkos::deep_copy(y_new_h, y_new);
    error(idx) = Kokkos::abs(y_new_h(0) - y_ref_h(0)) / Kokkos::abs(y_ref_h(0));

#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    scalar_type dt = (tend - tstart) / num_steps(idx);
    std::cout << "dt=" << dt << ", error=" << error(idx) << ", solution: {" << y_new_h(0) << ", " << y_new_h(1) << "}"
              << std::endl;
#endif
  }

}  // test_method

template <class Device>
void test_convergence_rate() {
  using execution_space = typename Device::execution_space;
  using RK_type         = KokkosODE::Experimental::RK_type;
  using vec_type        = Kokkos::View<double*, Device>;
  using mv_type         = Kokkos::View<double**, Device>;

  duho my_oscillator(1, 1, 4);
  const int neqs = my_oscillator.neqs;

  vec_type y("solution", neqs), f("function", neqs);
  auto y_h = Kokkos::create_mirror(y);
  y_h(0)   = 1;
  y_h(1)   = 0;
  Kokkos::deep_copy(y, y_h);

  constexpr double tstart = 0, tend = 1.024;
  Kokkos::View<int*, Kokkos::HostSpace> num_steps("Max Steps", 8);
  num_steps(0) = 512;
  num_steps(1) = 256;
  num_steps(2) = 128;
  num_steps(3) = 64;
  num_steps(4) = 32;
  num_steps(5) = 16;
  num_steps(6) = 8;
  num_steps(7) = 4;
  vec_type y_new("y new", neqs), y_old("y old", neqs);

  // Since y_old_h will be reused to set initial conditions
  // for each method tested we do not want to use
  // create_mirror_view which would not do a copy
  // when y_old is in HostSpace.
  typename vec_type::HostMirror y_old_h = Kokkos::create_mirror(y_old);
  y_old_h(0)                            = 1;
  y_old_h(1)                            = 0;

  // First compute analytical solution as reference
  // and to evaluate the error from each RK method.
  vec_type y_ref("reference value", neqs);
  auto y_ref_h = Kokkos::create_mirror(y_ref);
  {
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    solution_wrapper wrapper(my_oscillator, tend, y_old, y_ref);
    Kokkos::parallel_for(my_policy, wrapper);

    Kokkos::deep_copy(y_ref_h, y_ref);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    std::cout << "\nAnalytical solution" << std::endl;
    std::cout << "  y={" << y_ref_h(0) << ", " << y_ref_h(1) << "}" << std::endl;
#endif
  }

  typename vec_type::HostMirror error("error", num_steps.extent(0));
  test_rate<duho, RK_type::RKEH, vec_type, mv_type, double>(my_oscillator, tstart, tend, num_steps, y_old_h, y_ref_h,
                                                            error);

  for (int idx = 1; idx < num_steps.extent_int(0) - 2; ++idx) {
    double expected_ratio =
        Kokkos::pow(num_steps(idx) / num_steps(idx + 1), KokkosODE::Impl::ButcherTableau<1, 1>::order);
    double actual_ratio = error(idx + 1) / error(idx);
    EXPECT_NEAR_KK_REL(actual_ratio, expected_ratio, 0.15);

#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    double rel_ratio_diff = Kokkos::abs(actual_ratio - expected_ratio) / Kokkos::abs(expected_ratio);
    std::cout << "error ratio: " << actual_ratio << ", expected ratio: " << expected_ratio
              << ", rel diff: " << rel_ratio_diff << std::endl;
#endif
  }

  Kokkos::deep_copy(error, 0);
  test_rate<duho, RK_type::RKBS, vec_type, mv_type, double>(my_oscillator, tstart, tend, num_steps, y_old_h, y_ref_h,
                                                            error);

  for (int idx = 1; idx < num_steps.extent_int(0) - 2; ++idx) {
    double expected_ratio =
        Kokkos::pow(num_steps(idx) / num_steps(idx + 1), KokkosODE::Impl::ButcherTableau<2, 3>::order);
    double actual_ratio = error(idx + 1) / error(idx);
    EXPECT_NEAR_KK_REL(actual_ratio, expected_ratio, 0.05);

#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    double rel_ratio_diff = Kokkos::abs(actual_ratio - expected_ratio) / Kokkos::abs(expected_ratio);
    std::cout << "error ratio: " << actual_ratio << ", expected ratio: " << expected_ratio
              << ", rel diff: " << rel_ratio_diff << std::endl;
#endif
  }

  Kokkos::deep_copy(error, 0);
  test_rate<duho, RK_type::RKF45, vec_type, mv_type, double>(my_oscillator, tstart, tend, num_steps, y_old_h, y_ref_h,
                                                             error);

  for (int idx = 1; idx < num_steps.extent_int(0) - 2; ++idx) {
    double expected_ratio =
        Kokkos::pow(num_steps(idx) / num_steps(idx + 1), KokkosODE::Impl::ButcherTableau<4, 5>::order);
    double actual_ratio = error(idx + 1) / error(idx);
    EXPECT_NEAR_KK_REL(actual_ratio, expected_ratio, 0.05);

#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    double rel_ratio_diff = Kokkos::abs(actual_ratio - expected_ratio) / Kokkos::abs(expected_ratio);
    std::cout << "error ratio: " << actual_ratio << ", expected ratio: " << expected_ratio
              << ", rel diff: " << rel_ratio_diff << std::endl;
#endif
  }
}  // test_convergence_rate

template <class Device>
void test_adaptivity() {
  using execution_space = typename Device::execution_space;
  using RK_type         = KokkosODE::Experimental::RK_type;
  using vec_type        = Kokkos::View<double*, Device>;
  using mv_type         = Kokkos::View<double**, Device>;

  duho my_oscillator(1, 1, 4);
  const int neqs = my_oscillator.neqs;

  vec_type y("solution", neqs), f("function", neqs);
  auto y_h = Kokkos::create_mirror(y);
  y_h(0)   = 1;
  y_h(1)   = 0;
  Kokkos::deep_copy(y, y_h);

  constexpr double tstart = 0, tend = 1.024;
  constexpr int maxSteps = 512, numSteps = 128;
  constexpr double absTol = 1e-14, relTol = 1e-8, minStepSize = 0.001;
  vec_type y_new("y new", neqs), y_old("y old", neqs);

  // Since y_old_h will be reused to set initial conditions
  // for each method tested we do not want to use
  // create_mirror_view which would not do a copy
  // when y_old is in HostSpace.
  typename vec_type::HostMirror y_old_h = Kokkos::create_mirror(y_old);
  y_old_h(0)                            = 1;
  y_old_h(1)                            = 0;

  // First compute analytical solution as reference
  // and to evaluate the error from each RK method.
  vec_type y_ref("reference value", neqs);
  auto y_ref_h = Kokkos::create_mirror(y_ref);
  {
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    solution_wrapper wrapper(my_oscillator, tend, y_old, y_ref);
    Kokkos::parallel_for(my_policy, wrapper);

    Kokkos::deep_copy(y_ref_h, y_ref);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    std::cout << "\nAnalytical solution" << std::endl;
    std::cout << "  y={" << y_ref_h(0) << ", " << y_ref_h(1) << "}" << std::endl;
#endif
  }

  vec_type tmp("tmp vector", neqs);
  mv_type kstack("k stack", KokkosODE::Experimental::RungeKutta<RK_type::RKF45>::num_stages(), neqs);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  KokkosODE::Experimental::ODE_params params(numSteps, maxSteps, absTol, relTol, minStepSize);
  Kokkos::deep_copy(y_old, y_old_h);
  Kokkos::deep_copy(y_new, y_old_h);
  RKSolve_wrapper<duho, RK_type::RKF45, vec_type, mv_type, double> solve_wrapper(my_oscillator, params, tstart, tend,
                                                                                 y_old, y_new, tmp, kstack);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto y_new_h = Kokkos::create_mirror(y_new);
  Kokkos::deep_copy(y_new_h, y_new);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "Results: " << std::endl;
  std::cout << "  y_ref={ ";
  for (int idx = 0; idx < y_ref_h.extent_int(0); ++idx) {
    std::cout << y_ref_h(idx) << " ";
  }
  std::cout << "}" << std::endl;
  std::cout << "  y_new={ ";
  for (int idx = 0; idx < y_new_h.extent_int(0); ++idx) {
    std::cout << y_new_h(idx) << " ";
  }
  std::cout << "}" << std::endl;
  std::cout << "  error={ ";
  double error;
#endif

  for (int idx = 0; idx < y_new_h.extent_int(0); ++idx) {
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    error = Kokkos::abs(y_new_h(idx) - y_ref_h(idx)) / Kokkos::abs(y_ref_h(idx));
    std::cout << error << " ";
#endif
    EXPECT_NEAR_KK_REL(y_new_h(idx), y_ref_h(idx), 1e-7);
  }
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "}" << std::endl;
#endif

}  // test_adaptivity

}  // namespace Test

void test_RK() { Test::test_RK<TestDevice>(); }

void test_RK_conv_rate() { Test::test_convergence_rate<TestDevice>(); }

void test_RK_adaptivity() { Test::test_adaptivity<TestDevice>(); }

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, RKSolve_serial) { test_RK(); }
TEST_F(TestCategory, RK_conv_rate) { test_RK_conv_rate(); }
TEST_F(TestCategory, RK_adaptivity) { test_RK_adaptivity(); }
#endif
