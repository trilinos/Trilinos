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

// R1 = 1e-6*1.85e10 * exp(-15618 / T) * (reac) ( 1 – (1- 10^-9) reac)
// d(reac)/dt = -R1
// d(prod)/dt = R1
struct chem_model_1 {
  constexpr static int neqs = 2;
  // constexpr static double alpha = 1e-6*1.85e10;
  constexpr static double alpha = 1.85e10;
  constexpr static double beta  = 15618;
  constexpr static double gamma = 1 - 10e-9;

  const double tstart, tend, T0, T1;

  chem_model_1(const double tstart_ = 0, const double tend_ = 100, const double T0_ = 300, const double T1_ = 800)
      : tstart(tstart_), tend(tend_), T0(T0_), T1(T1_){};

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double t, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    // First compute the temperature
    // using linear ramp from T0 to T1
    // between tstart and tend.
    double T = (T1 - T0) * (t - tstart) / (tend - tstart) + T0;

    // Evaluate the chemical reaction rate
    f(0) = -alpha * Kokkos::exp(-beta / T) * y(0) * (1 - gamma * y(0));
    f(1) = -f(0);
  }
};

struct chem_model_2 {
  constexpr static int neqs      = 7;
  constexpr static double alpha1 = 1e-6 * 3334169440721739.0 * 1500;
  constexpr static double beta1  = 207850000.0 / 8314.0;
  constexpr static double alpha2 = 1e-6 * 49997793980831.89 * 1500;
  constexpr static double beta2  = 207850000.0 / 8314.0;

  const double tstart, tend, T0, T1;

  chem_model_2(const double tstart_ = 0, const double tend_ = 1200, const double T0_ = 300, const double T1_ = 1000)
      : tstart(tstart_), tend(tend_), T0(T0_), T1(T1_){};

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double t, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    // First compute the temperature
    // using linear ramp from T0 to T1
    // between tstart and tend.
    double T = (T1 - T0) * (t - tstart) / (1500 - tstart) + T0;

    // Evaluate the chemical reaction rates
    double R1 = y(0) * alpha1 * Kokkos::exp(-beta1 / T);
    double R2 = y(1) * alpha2 * Kokkos::exp(-beta2 / T);

    // Evaluate the chemical reaction rate
    f(0) = -R1;
    f(1) = -R2;
    f(2) = R1 + 0.08 * R2;
    f(3) = 0.147 * R2;
    f(4) = 0.453 * R2;
    f(5) = 0.187 * R2;
    f(6) = 0.133 * R2;
  }
};

template <class Device>
void test_chem() {
  using execution_space = typename Device::execution_space;
  using vec_type        = Kokkos::View<double*, Device>;
  using mv_type         = Kokkos::View<double**, Device>;
  using RK_type         = KokkosODE::Experimental::RK_type;
  using solver_type     = KokkosODE::Experimental::RungeKutta<RK_type::RKCK>;

  {
    chem_model_1 chem_model;
    const int neqs      = chem_model.neqs;
    const int num_steps = 15000;

    KokkosODE::Experimental::ODE_params params(num_steps);
    vec_type tmp("tmp vector", neqs);
    mv_type kstack("k stack", solver_type::num_stages(), neqs);

    // Set initial conditions
    vec_type y_new("solution", neqs);
    vec_type y_old("initial conditions", neqs);
    auto y_old_h = Kokkos::create_mirror(y_old);
    y_old_h(0)   = 1;
    y_old_h(1)   = 0;
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::deep_copy(y_new, y_old_h);

    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    RKSolve_wrapper<chem_model_1, RK_type::RKCK, vec_type, mv_type, double> solve_wrapper(
        chem_model, params, chem_model.tstart, chem_model.tend, y_old, y_new, tmp, kstack);
    Kokkos::parallel_for(my_policy, solve_wrapper);

    auto y_new_h = Kokkos::create_mirror(y_new);
    Kokkos::deep_copy(y_new_h, y_new);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    const double dt = (chem_model.tend - chem_model.tstart) / params.num_steps;
    std::cout << "\nChem model 1" << std::endl;
    std::cout << "  t0=" << chem_model.tstart << ", tn=" << chem_model.tend << std::endl;
    std::cout << "  T0=" << chem_model.T0 << ", Tn=" << chem_model.T1 << std::endl;
    std::cout << "  dt=" << dt << std::endl;
    std::cout << "  y(t0)={" << y_old_h(0) << ", " << y_old_h(1) << "}" << std::endl;
    std::cout << "  y(tn)={" << y_new_h(0) << ", " << y_new_h(1) << "}" << std::endl;
#endif
  }

  {
    chem_model_2 chem_model;
    const int neqs      = chem_model.neqs;
    const int num_steps = 1500;

    KokkosODE::Experimental::ODE_params params(num_steps);
    vec_type tmp("tmp vector", neqs);
    mv_type kstack("k stack", solver_type::num_stages(), neqs);

    // Set initial conditions
    vec_type y_new("solution", neqs);
    vec_type y_old("initial conditions", neqs);
    auto y_old_h = Kokkos::create_mirror(y_old);
    y_old_h(0)   = 0.25;
    y_old_h(1)   = 0.25;
    y_old_h(2)   = 0;
    y_old_h(3)   = 0;
    y_old_h(4)   = 0;
    y_old_h(5)   = 0;
    y_old_h(6)   = 0;
    Kokkos::deep_copy(y_old, y_old_h);
    Kokkos::deep_copy(y_new, y_old_h);

    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    RKSolve_wrapper<chem_model_2, RK_type::RKCK, vec_type, mv_type, double> solve_wrapper(
        chem_model, params, chem_model.tstart, chem_model.tend, y_old, y_new, tmp, kstack);
    Kokkos::parallel_for(my_policy, solve_wrapper);

    auto y_new_h = Kokkos::create_mirror(y_new);
    Kokkos::deep_copy(y_new_h, y_new);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
    const double dt = (chem_model.tend - chem_model.tstart) / params.num_steps;
    std::cout << "\nChem model 2" << std::endl;
    std::cout << "  t0=" << chem_model.tstart << ", tn=" << chem_model.tend << std::endl;
    std::cout << "  T0=" << chem_model.T0 << ", Tn=" << chem_model.T1 << std::endl;
    std::cout << "  dt=" << dt << std::endl;
    std::cout << "  y(t0)={" << y_old_h(0) << ", " << y_old_h(1) << ", " << y_old_h(2) << ", " << y_old_h(3) << ", "
              << y_old_h(4) << ", " << y_old_h(5) << ", " << y_old_h(6) << "}" << std::endl;
    std::cout << "  y(tn)={" << y_new_h(0) << ", " << y_new_h(1) << ", " << y_new_h(2) << ", " << y_new_h(3) << ", "
              << y_new_h(4) << ", " << y_new_h(5) << ", " << y_new_h(6) << "}" << std::endl;
#endif
  }
}  // test_chem
}  // namespace Test

int test_chem_models() {
  Test::test_chem<TestDevice>();

  return 1;
}

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, RK_chem_models) { test_chem_models(); }
#endif
