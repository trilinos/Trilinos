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

#include "KokkosODE_BDF.hpp"

namespace Test {

// Logistic equation
// Used to model population growth
// it is a simple nonlinear ODE with
// a lot of literature.
//
// Equation: y'(t) = r*y(t)*(1-y(t)/K)
// Jacobian: df/dy = r - 2*r*y/K
// Solution: y = K / (1 + ((K - y0) / y0)*exp(-rt))
struct Logistic {
  static constexpr int neqs = 1;

  const double r, K;

  Logistic(double r_, double K_) : r(r_), K(K_){};

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    f(0) = r * y(0) * (1.0 - y(0) / K);
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/, const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
    jac(0, 0) = r - 2 * r * y(0) / K;
  }

  template <class vec_type>
  KOKKOS_FUNCTION void solution(const double t, const vec_type& y0, const vec_type& y) const {
    y(0) = K / (1 + (K - y0) / y0 * Kokkos::exp(-r * t));
  }

};  // Logistic

// Lotka-Volterra equation
// A predator-prey model that describe
// population dynamics when two species
// interact.
//
// Equations: y0'(t) = alpha*y0(t) - beta*y0(t)*y1(t)
//            y1'(t) = delta*y0(t)*y1(t) - gamma*y1(t)
// Jacobian: df0/dy = [alpha-beta*y1(t);     beta*y0(t)]
//           df1/dy = [delta*y1(t);          delta*y0(t)-gamma]
// Solution: y = K / (1 + ((K - y0) / y0)*exp(-rt))
struct LotkaVolterra {
  static constexpr int neqs = 2;

  const double alpha, beta, delta, gamma;

  LotkaVolterra(double alpha_, double beta_, double delta_, double gamma_)
      : alpha(alpha_), beta(beta_), delta(delta_), gamma(gamma_){};

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    f(0) = alpha * y(0) - beta * y(0) * y(1);
    f(1) = delta * y(0) * y(1) - gamma * y(1);
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/, const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
    jac(0, 0) = alpha - beta * y(1);
    jac(0, 1) = -beta * y(0);
    jac(1, 0) = delta * y(1);
    jac(1, 1) = delta * y(0) - gamma;
  }

};  // LotkaVolterra

// Robertson's autocatalytic chemical reaction:
// H. H. Robertson, The solution of a set of reaction rate equations,
// in J. Walsh (Ed.), Numerical Analysis: An Introduction,
// pp. 178–182, Academic Press, London (1966).
//
// Equations: y0' = -0.04*y0 + 10e4*y1*y2
//            y1' =  0.04*y0 - 10e4*y1*y2 - 3e7 * y1**2
//            y2' = 3e7 * y1**2
struct StiffChemistry {
  static constexpr int neqs = 3;

  StiffChemistry() {}

  template <class vec_type1, class vec_type2>
  KOKKOS_FUNCTION void evaluate_function(const double /*t*/, const double /*dt*/, const vec_type1& y,
                                         const vec_type2& f) const {
    f(0) = -0.04 * y(0) + 1.e4 * y(1) * y(2);
    f(1) = 0.04 * y(0) - 1.e4 * y(1) * y(2) - 3.e7 * y(1) * y(1);
    f(2) = 3.e7 * y(1) * y(1);
  }

  template <class vec_type, class mat_type>
  KOKKOS_FUNCTION void evaluate_jacobian(const double /*t*/, const double /*dt*/, const vec_type& y,
                                         const mat_type& jac) const {
    jac(0, 0) = -0.04;
    jac(0, 1) = 1.e4 * y(2);
    jac(0, 2) = 1.e4 * y(1);
    jac(1, 0) = 0.04;
    jac(1, 1) = -1.e4 * y(2) - 3.e7 * 2.0 * y(1);
    jac(1, 2) = -1.e4 * y(1);
    jac(2, 0) = 0.0;
    jac(2, 1) = 3.e7 * 2.0 * y(1);
    jac(2, 2) = 0.0;
  }
};

template <class ode_type, KokkosODE::Experimental::BDF_type bdf_type, class vec_type, class mv_type, class mat_type,
          class scalar_type>
struct BDFSolve_wrapper {
  ode_type my_ode;
  scalar_type tstart, tend;
  int num_steps;
  vec_type y_old, y_new, rhs, update, scale;
  mv_type y_vecs, kstack;
  mat_type temp, jac;

  BDFSolve_wrapper(const ode_type& my_ode_, const scalar_type tstart_, const scalar_type tend_, const int num_steps_,
                   const vec_type& y_old_, const vec_type& y_new_, const vec_type& rhs_, const vec_type& update_,
                   const vec_type& scale_, const mv_type& y_vecs_, const mv_type& kstack_, const mat_type& temp_,
                   const mat_type& jac_)
      : my_ode(my_ode_),
        tstart(tstart_),
        tend(tend_),
        num_steps(num_steps_),
        y_old(y_old_),
        y_new(y_new_),
        rhs(rhs_),
        update(update_),
        scale(scale_),
        y_vecs(y_vecs_),
        kstack(kstack_),
        temp(temp_),
        jac(jac_) {}

  KOKKOS_FUNCTION
  void operator()(const int /*idx*/) const {
    KokkosODE::Experimental::BDF<bdf_type>::Solve(my_ode, tstart, tend, num_steps, y_old, y_new, rhs, update, scale,
                                                  y_vecs, kstack, temp, jac);
  }
};

template <class ode_type, class mat_type, class vec_type, class scalar_type>
struct BDF_Solve_wrapper {
  const ode_type my_ode;
  const scalar_type t_start, t_end, dt, max_step;
  const vec_type y0, y_new;
  const mat_type temp, temp2;

  BDF_Solve_wrapper(const ode_type& my_ode_, const scalar_type& t_start_, const scalar_type& t_end_,
                    const scalar_type& dt_, const scalar_type& max_step_, const vec_type& y0_, const vec_type& y_new_,
                    const mat_type& temp_, const mat_type& temp2_)
      : my_ode(my_ode_),
        t_start(t_start_),
        t_end(t_end_),
        dt(dt_),
        max_step(max_step_),
        y0(y0_),
        y_new(y_new_),
        temp(temp_),
        temp2(temp2_) {}

  KOKKOS_FUNCTION void operator()(const int) const {
    KokkosODE::Experimental::BDFSolve(my_ode, t_start, t_end, dt, max_step, y0, y_new, temp, temp2);
  }
};

template <class device_type, class scalar_type>
void test_BDF_Logistic() {
  using execution_space = typename device_type::execution_space;
  using vec_type        = Kokkos::View<scalar_type*, execution_space>;
  using mv_type         = Kokkos::View<scalar_type**, execution_space>;
  using mat_type        = Kokkos::View<scalar_type**, execution_space>;

  Kokkos::RangePolicy<execution_space> myPolicy(0, 1);
  Logistic mySys(1, 1);

  constexpr int num_tests   = 7;
  int num_steps[num_tests]  = {512, 256, 128, 64, 32, 16, 8};
  double errors[num_tests]  = {0};
  const scalar_type t_start = 0.0, t_end = 6.0;
  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);
  vec_type rhs("rhs", mySys.neqs), update("update", mySys.neqs);
  vec_type scale("scaling factors", mySys.neqs);
  mat_type jac("jacobian", mySys.neqs, mySys.neqs), temp("temp storage", mySys.neqs, mySys.neqs + 4);
  mv_type kstack("Startup RK vectors", 6, mySys.neqs);

  Kokkos::deep_copy(scale, 1);

  scalar_type measured_order;

  // Test BDF1
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\nBDF1 convergence test" << std::endl;
#endif
  for (int idx = 0; idx < num_tests; idx++) {
    mv_type y_vecs("history vectors", mySys.neqs, 1);

    Kokkos::deep_copy(y0, 0.5);
    Kokkos::deep_copy(y_vecs, 0.5);

    BDFSolve_wrapper<Logistic, KokkosODE::Experimental::BDF_type::BDF1, vec_type, mv_type, mat_type, scalar_type>
        solve_wrapper(mySys, t_start, t_end, num_steps[idx], y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
    Kokkos::parallel_for(myPolicy, solve_wrapper);
    Kokkos::fence();

    auto y_new_h = Kokkos::create_mirror_view(y_new);
    Kokkos::deep_copy(y_new_h, y_new);

    errors[idx] = Kokkos::abs(y_new_h(0) - 1 / (1 + Kokkos::exp(-t_end))) / Kokkos::abs(1 / (1 + Kokkos::exp(-t_end)));
  }
  measured_order = Kokkos::pow(errors[num_tests - 1] / errors[0], 1.0 / (num_tests - 1));
  EXPECT_NEAR_KK_REL(measured_order, 2.0, 0.15);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "expected ratio: 2, actual ratio: " << measured_order
            << ", order error=" << Kokkos::abs(measured_order - 2.0) / 2.0 << std::endl;
#endif

  // Test BDF2
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\nBDF2 convergence test" << std::endl;
#endif
  for (int idx = 0; idx < num_tests; idx++) {
    mv_type y_vecs("history vectors", mySys.neqs, 2);
    Kokkos::deep_copy(y0, 0.5);

    BDFSolve_wrapper<Logistic, KokkosODE::Experimental::BDF_type::BDF2, vec_type, mv_type, mat_type, scalar_type>
        solve_wrapper(mySys, t_start, t_end, num_steps[idx], y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
    Kokkos::parallel_for(myPolicy, solve_wrapper);
    Kokkos::fence();

    auto y_new_h = Kokkos::create_mirror_view(y_new);
    Kokkos::deep_copy(y_new_h, y_new);

    errors[idx] = Kokkos::abs(y_new_h(0) - 1 / (1 + Kokkos::exp(-t_end))) / Kokkos::abs(1 / (1 + Kokkos::exp(-t_end)));
  }
  measured_order = Kokkos::pow(errors[num_tests - 1] / errors[0], 1.0 / (num_tests - 1));
  EXPECT_NEAR_KK_REL(measured_order, 4.0, 0.15);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "expected ratio: 4, actual ratio: " << measured_order
            << ", order error=" << Kokkos::abs(measured_order - 4.0) / 4.0 << std::endl;
#endif

  // Test BDF3
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\nBDF3 convergence test" << std::endl;
#endif
  for (int idx = 0; idx < num_tests; idx++) {
    mv_type y_vecs("history vectors", mySys.neqs, 3);
    Kokkos::deep_copy(y0, 0.5);

    BDFSolve_wrapper<Logistic, KokkosODE::Experimental::BDF_type::BDF3, vec_type, mv_type, mat_type, scalar_type>
        solve_wrapper(mySys, t_start, t_end, num_steps[idx], y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
    Kokkos::parallel_for(myPolicy, solve_wrapper);
    Kokkos::fence();

    auto y_new_h = Kokkos::create_mirror_view(y_new);
    Kokkos::deep_copy(y_new_h, y_new);

    errors[idx] = Kokkos::abs(y_new_h(0) - 1 / (1 + Kokkos::exp(-t_end))) / Kokkos::abs(1 / (1 + Kokkos::exp(-t_end)));
  }
  measured_order = Kokkos::pow(errors[num_tests - 1] / errors[0], 1.0 / (num_tests - 1));
  EXPECT_NEAR_KK_REL(measured_order, 8.0, 0.15);
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "expected ratio: 8, actual ratio: " << measured_order
            << ", order error=" << Kokkos::abs(measured_order - 8.0) / 8.0 << std::endl;
#endif

  // Test BDF4
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\nBDF4 convergence test" << std::endl;
#endif
  for (int idx = 0; idx < num_tests; idx++) {
    mv_type y_vecs("history vectors", mySys.neqs, 4);
    Kokkos::deep_copy(y0, 0.5);

    BDFSolve_wrapper<Logistic, KokkosODE::Experimental::BDF_type::BDF4, vec_type, mv_type, mat_type, scalar_type>
        solve_wrapper(mySys, t_start, t_end, num_steps[idx], y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
    Kokkos::parallel_for(myPolicy, solve_wrapper);
    Kokkos::fence();

    auto y_new_h = Kokkos::create_mirror_view(y_new);
    Kokkos::deep_copy(y_new_h, y_new);

    errors[idx] = Kokkos::abs(y_new_h(0) - 1 / (1 + Kokkos::exp(-t_end))) / Kokkos::abs(1 / (1 + Kokkos::exp(-t_end)));
  }
  measured_order = Kokkos::pow(errors[num_tests - 1] / errors[0], 1.0 / (num_tests - 1));
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "expected ratio: 16, actual ratio: " << measured_order
            << ", order error=" << Kokkos::abs(measured_order - 16.0) / 16.0 << std::endl;
#endif

  // Test BDF5
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "\nBDF5 convergence test" << std::endl;
#endif
  for (int idx = 0; idx < num_tests; idx++) {
    mv_type y_vecs("history vectors", mySys.neqs, 5);
    Kokkos::deep_copy(y0, 0.5);

    BDFSolve_wrapper<Logistic, KokkosODE::Experimental::BDF_type::BDF5, vec_type, mv_type, mat_type, scalar_type>
        solve_wrapper(mySys, t_start, t_end, num_steps[idx], y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
    Kokkos::parallel_for(myPolicy, solve_wrapper);
    Kokkos::fence();

    auto y_new_h = Kokkos::create_mirror_view(y_new);
    Kokkos::deep_copy(y_new_h, y_new);

    errors[idx] = Kokkos::abs(y_new_h(0) - 1 / (1 + Kokkos::exp(-t_end))) / Kokkos::abs(1 / (1 + Kokkos::exp(-t_end)));
  }
  measured_order = Kokkos::pow(errors[num_tests - 1] / errors[0], 1.0 / (num_tests - 1));
#if defined(HAVE_KOKKOSKERNELS_DEBUG)
  std::cout << "expected ratio: 32, actual ratio: " << measured_order
            << ", order error=" << Kokkos::abs(measured_order - 32.0) / 32.0 << std::endl;
#endif

}  // test_BDF_Logistic

template <class device_type, class scalar_type>
void test_BDF_LotkaVolterra() {
  using execution_space = typename device_type::execution_space;
  using vec_type        = Kokkos::View<scalar_type*, execution_space>;
  using mv_type         = Kokkos::View<scalar_type**, execution_space>;
  using mat_type        = Kokkos::View<scalar_type**, execution_space>;

  LotkaVolterra mySys(1.1, 0.4, 0.1, 0.4);

  const scalar_type t_start = 0.0, t_end = 100.0;
  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);
  vec_type rhs("rhs", mySys.neqs), update("update", mySys.neqs);
  vec_type scale("scaling factors", mySys.neqs);
  mat_type jac("jacobian", mySys.neqs, mySys.neqs), temp("temp storage", mySys.neqs, mySys.neqs + 4);

  Kokkos::deep_copy(scale, 1);

  // Test BDF5
  mv_type kstack("Startup RK vectors", 6, mySys.neqs);
  mv_type y_vecs("history vectors", mySys.neqs, 5);

  Kokkos::deep_copy(y0, 10.0);
  Kokkos::deep_copy(y_vecs, 10.0);

  Kokkos::RangePolicy<execution_space> myPolicy(0, 1);
  BDFSolve_wrapper<LotkaVolterra, KokkosODE::Experimental::BDF_type::BDF5, vec_type, mv_type, mat_type, scalar_type>
      solve_wrapper(mySys, t_start, t_end, 1000, y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
  Kokkos::parallel_for(myPolicy, solve_wrapper);
}

template <class device_type, class scalar_type>
void test_BDF_StiffChemistry() {
  using execution_space = typename device_type::execution_space;
  using vec_type        = Kokkos::View<scalar_type*, execution_space>;
  using mv_type         = Kokkos::View<scalar_type**, execution_space>;
  using mat_type        = Kokkos::View<scalar_type**, execution_space>;

  StiffChemistry mySys{};

  const scalar_type t_start = 0.0, t_end = 500.0;
  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);
  vec_type rhs("rhs", mySys.neqs), update("update", mySys.neqs);
  vec_type scale("scaling factors", mySys.neqs);
  mat_type jac("jacobian", mySys.neqs, mySys.neqs), temp("temp storage", mySys.neqs, mySys.neqs + 4);

  Kokkos::deep_copy(scale, 1);

  // Test BDF5
  mv_type kstack("Startup RK vectors", 6, mySys.neqs);
  mv_type y_vecs("history vectors", mySys.neqs, 5);

  auto y0_h = Kokkos::create_mirror_view(y0);
  y0_h(0)   = 1.0;
  y0_h(1)   = 0.0;
  y0_h(2)   = 0.0;
  Kokkos::deep_copy(y0, y0_h);
  Kokkos::deep_copy(y_vecs, 0.0);

  Kokkos::RangePolicy<execution_space> myPolicy(0, 1);
  BDFSolve_wrapper<StiffChemistry, KokkosODE::Experimental::BDF_type::BDF5, vec_type, mv_type, mat_type, scalar_type>
      solve_wrapper(mySys, t_start, t_end, 110000, y0, y_new, rhs, update, scale, y_vecs, kstack, temp, jac);
  Kokkos::parallel_for(myPolicy, solve_wrapper);
}

// template <class ode_type, KokkosODE::Experimental::BDF_type bdf_type,
//           class vec_type, class mv_type, class mat_type, class scalar_type>
// struct BDFSolve_parallel {
//   ode_type my_ode;
//   scalar_type tstart, tend;
//   int num_steps;
//   vec_type y_old, y_new, rhs, update, scale;
//   mv_type y_vecs, kstack;
//   mat_type temp, jac;

//   BDFSolve_parallel(const ode_type& my_ode_, const scalar_type tstart_,
//                     const scalar_type tend_, const int num_steps_,
//                     const vec_type& y_old_, const vec_type& y_new_,
//                     const vec_type& rhs_, const vec_type& update_,
// 		    const vec_type& scale_,
//                     const mv_type& y_vecs_, const mv_type& kstack_,
//                     const mat_type& temp_, const mat_type& jac_)
//       : my_ode(my_ode_),
//         tstart(tstart_),
//         tend(tend_),
//         num_steps(num_steps_),
//         y_old(y_old_),
//         y_new(y_new_),
//         rhs(rhs_),
//         update(update_),
// 	scale(scale_),
//         y_vecs(y_vecs_),
//         kstack(kstack_),
//         temp(temp_),
//         jac(jac_) {}

//   KOKKOS_FUNCTION
//   void operator()(const int idx) const {
//     auto local_y_old = Kokkos::subview(
//         y_old,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)));
//     auto local_y_new = Kokkos::subview(
//         y_new,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)));
//     auto local_rhs = Kokkos::subview(
//         rhs,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)));
//     auto local_update = Kokkos::subview(
//         update,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)));

//     auto local_y_vecs = Kokkos::subview(
//         y_vecs,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)),
//         Kokkos::ALL());
//     auto local_kstack = Kokkos::subview(
//         kstack, Kokkos::ALL(),
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)));
//     auto local_temp = Kokkos::subview(
//         temp,
//         Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx + 1)),
//         Kokkos::ALL());
//     auto local_jac = Kokkos::subview(
//         jac, Kokkos::pair<int, int>(my_ode.neqs * idx, my_ode.neqs * (idx +
//         1)), Kokkos::ALL());

//     KokkosODE::Experimental::BDF<bdf_type>::Solve(
//         my_ode, tstart, tend, num_steps, local_y_old, local_y_new, local_rhs,
//         local_update, scale, local_y_vecs, local_kstack, local_temp,
//         local_jac);
//   }
// };

// template <class device_type, class scalar_type>
// void test_BDF_parallel() {
//   using execution_space = typename device_type::execution_space;
//   using vec_type = Kokkos::View<scalar_type*, execution_space>;
//   using mv_type  = Kokkos::View<scalar_type**, execution_space>;
//   using mat_type = Kokkos::View<scalar_type**, execution_space>;

//   LotkaVolterra mySys(1.1, 0.4, 0.1, 0.4);
//   constexpr int num_solves = 1000;

//   vec_type scale("scaling factors", mySys.neqs);
//   Kokkos::deep_copy(scale, 1);

//   const scalar_type t_start = 0.0, t_end = 100.0;
//   vec_type y0("initial conditions", mySys.neqs * num_solves),
//       y_new("solution", mySys.neqs * num_solves);
//   vec_type rhs("rhs", mySys.neqs * num_solves),
//       update("update", mySys.neqs * num_solves);
//   mat_type jac("jacobian", mySys.neqs * num_solves, mySys.neqs),
//       temp("temp storage", mySys.neqs * num_solves, mySys.neqs + 4);

//   // Test BDF5
//   mv_type y_vecs("history vectors", mySys.neqs * num_solves, 5),
//       kstack("Startup RK vectors", 6, mySys.neqs * num_solves);

//   Kokkos::deep_copy(y0, 10.0);
//   Kokkos::deep_copy(y_vecs, 10.0);

//   Kokkos::RangePolicy<execution_space> myPolicy(0, num_solves);
//   BDFSolve_parallel<LotkaVolterra, KokkosODE::Experimental::BDF_type::BDF5,
//                     vec_type, mv_type, mat_type, scalar_type>
//       solve_wrapper(mySys, t_start, t_end, 1000, y0, y_new, rhs, update,
//       scale, y_vecs,
//                     kstack, temp, jac);
//   Kokkos::parallel_for(myPolicy, solve_wrapper);

//   Kokkos::fence();
// }

template <class mat_type, class scalar_type>
void compute_coeffs(const int order, const scalar_type factor, const mat_type& coeffs) {
  std::cout << "compute_coeffs" << std::endl;

  coeffs(0, 0) = 1.0;
  for (int colIdx = 0; colIdx < order; ++colIdx) {
    coeffs(0, colIdx + 1) = 1.0;
    for (int rowIdx = 0; rowIdx < order; ++rowIdx) {
      coeffs(rowIdx + 1, colIdx + 1) =
          ((rowIdx - factor * (colIdx + 1.0)) / (rowIdx + 1.0)) * coeffs(rowIdx, colIdx + 1);
    }
  }
}

template <class mat_type, class scalar_type>
void update_D(const int order, const scalar_type factor, const mat_type& coeffs, const mat_type& tempD,
              const mat_type& D) {
  auto subD     = Kokkos::subview(D, Kokkos::pair<int, int>(0, order + 1), Kokkos::ALL);
  auto subTempD = Kokkos::subview(tempD, Kokkos::pair<int, int>(0, order + 1), Kokkos::ALL);

  compute_coeffs(order, factor, coeffs);
  auto R = Kokkos::subview(coeffs, Kokkos::pair<int, int>(0, order + 1), Kokkos::pair<int, int>(0, order + 1));
  std::cout << "SerialGemm" << std::endl;
  KokkosBatched::SerialGemm<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose,
                            KokkosBatched::Algo::Gemm::Blocked>::invoke(1.0, R, subD, 0.0, subTempD);

  compute_coeffs(order, 1.0, coeffs);
  auto U = Kokkos::subview(coeffs, Kokkos::pair<int, int>(0, order + 1), Kokkos::pair<int, int>(0, order + 1));
  std::cout << "SerialGemm" << std::endl;
  KokkosBatched::SerialGemm<KokkosBatched::Trans::Transpose, KokkosBatched::Trans::NoTranspose,
                            KokkosBatched::Algo::Gemm::Blocked>::invoke(1.0, U, subTempD, 0.0, subD);
}

template <class device_type, class scalar_type>
void test_Nordsieck() {
  using execution_space = Kokkos::HostSpace;
  StiffChemistry mySys{};

  Kokkos::View<scalar_type**, execution_space> R("coeffs", 6, 6), U("coeffs", 6, 6);
  Kokkos::View<scalar_type**, execution_space> D("D", 8, mySys.neqs), tempD("tmp", 8, mySys.neqs);
  int order     = 1;
  double factor = 0.8;

  constexpr double t_start = 0.0, t_end = 500.0;
  int max_steps = 200000;
  double dt     = (t_end - t_start) / max_steps;

  auto y0 = Kokkos::subview(D, 0, Kokkos::ALL());
  auto f  = Kokkos::subview(D, 1, Kokkos::ALL());
  y0(0)   = 1.0;

  mySys.evaluate_function(0, 0, y0, f);
  for (int eqIdx = 0; eqIdx < mySys.neqs; ++eqIdx) {
    f(eqIdx) *= dt;
  }

  compute_coeffs(order, factor, R);
  compute_coeffs(order, 1.0, U);

  {
    std::cout << "R: " << std::endl;
    for (int i = 0; i < order + 1; ++i) {
      std::cout << "{ ";
      for (int j = 0; j < order + 1; ++j) {
        std::cout << R(i, j) << ", ";
      }
      std::cout << "}" << std::endl;
    }
  }

  std::cout << "D before update:" << std::endl;
  std::cout << "  { " << D(0, 0) << ", " << D(0, 1) << ", " << D(0, 2) << " }" << std::endl;
  std::cout << "  { " << D(1, 0) << ", " << D(1, 1) << ", " << D(1, 2) << " }" << std::endl;
  update_D(order, factor, R, tempD, D);

  std::cout << "D after update:" << std::endl;
  std::cout << "  { " << D(0, 0) << ", " << D(0, 1) << ", " << D(0, 2) << " }" << std::endl;
  std::cout << "  { " << D(1, 0) << ", " << D(1, 1) << ", " << D(1, 2) << " }" << std::endl;
}

template <class device_type, class scalar_type>
void test_adaptive_BDF() {
  using execution_space = typename device_type::execution_space;
  using vec_type        = Kokkos::View<scalar_type*, execution_space>;
  using mat_type        = Kokkos::View<scalar_type**, execution_space>;

  Logistic mySys(1, 1);

  constexpr double t_start = 0.0, t_end = 6.0, atol = 1.0e-6, rtol = 1.0e-4;
  constexpr int num_steps = 512, max_newton_iters = 5;
  int order = 1, num_equal_steps = 0;
  double dt = (t_end - t_start) / num_steps;
  double t  = t_start;

  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);
  vec_type rhs("rhs", mySys.neqs), update("update", mySys.neqs);
  mat_type temp("buffer1", mySys.neqs, 23 + 2 * mySys.neqs + 4), temp2("buffer2", 6, 7);

  // Initial condition
  Kokkos::deep_copy(y0, 0.5);

  // Initialize D
  auto D  = Kokkos::subview(temp, Kokkos::ALL(), Kokkos::pair<int, int>(2, 10));
  D(0, 0) = y0(0);
  mySys.evaluate_function(0, 0, y0, rhs);
  D(0, 1) = dt * rhs(0);
  Kokkos::deep_copy(rhs, 0);

  std::cout << "**********************\n"
            << "        Step 1\n"
            << "**********************" << std::endl;

  std::cout << "Initial conditions" << std::endl;
  std::cout << "   y0=" << y0(0) << ", t=" << t << ", dt=" << dt << std::endl;

  std::cout << "Initial D: {" << D(0, 0) << ", " << D(0, 1) << ", " << D(0, 2) << ", " << D(0, 3) << ", " << D(0, 4)
            << ", " << D(0, 5) << ", " << D(0, 6) << ", " << D(0, 7) << "}" << std::endl;

  KokkosODE::Impl::BDFStep(mySys, t, dt, t_end, order, num_equal_steps, max_newton_iters, atol, rtol, 0.2, y0, y_new,
                           rhs, update, temp, temp2);

  for (int eqIdx = 0; eqIdx < mySys.neqs; ++eqIdx) {
    y0(eqIdx) = y_new(eqIdx);
  }

  std::cout << "**********************\n"
            << "        Step 2\n"
            << "**********************" << std::endl;

  std::cout << "   y0=" << y0(0) << ", t=" << t << ", dt: " << dt << std::endl;

  std::cout << "Initial D: {" << D(0, 0) << ", " << D(0, 1) << ", " << D(0, 2) << ", " << D(0, 3) << ", " << D(0, 4)
            << ", " << D(0, 5) << ", " << D(0, 6) << ", " << D(0, 7) << "}" << std::endl;

  KokkosODE::Impl::BDFStep(mySys, t, dt, t_end, order, num_equal_steps, max_newton_iters, atol, rtol, 0.2, y0, y_new,
                           rhs, update, temp, temp2);

  for (int eqIdx = 0; eqIdx < mySys.neqs; ++eqIdx) {
    y0(eqIdx) = y_new(eqIdx);
  }

  std::cout << "**********************\n"
            << "        Step 3\n"
            << "**********************" << std::endl;

  std::cout << "   y0=" << y0(0) << ", t=" << t << ", dt: " << dt << std::endl;

  std::cout << "Initial D: {" << D(0, 0) << ", " << D(0, 1) << ", " << D(0, 2) << ", " << D(0, 3) << ", " << D(0, 4)
            << ", " << D(0, 5) << ", " << D(0, 6) << ", " << D(0, 7) << "}" << std::endl;

  KokkosODE::Impl::BDFStep(mySys, t, dt, t_end, order, num_equal_steps, max_newton_iters, atol, rtol, 0.2, y0, y_new,
                           rhs, update, temp, temp2);

  std::cout << "Final t: " << t << ", y=" << y_new(0) << std::endl;

}  // test_adaptive_BDF()

template <class execution_space, class scalar_type>
void test_adaptive_BDF_v2() {
  using vec_type = Kokkos::View<scalar_type*, execution_space>;
  using mat_type = Kokkos::View<scalar_type**, execution_space>;
  using KAT      = Kokkos::ArithTraits<scalar_type>;

  std::cout << "\n\n\nBDF_v2 test starting\n" << std::endl;

  Logistic mySys(1, 1);

  const scalar_type t_start = KAT::zero(),
                    t_end   = 6 * KAT::one();  //, atol = 1.0e-6, rtol = 1.0e-4;
  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);
  Kokkos::deep_copy(y0, 0.5);

  mat_type temp("buffer1", mySys.neqs, 23 + 2 * mySys.neqs + 4), temp2("buffer2", 6, 7);

  {
    scalar_type dt = KAT::zero();
    vec_type f0("initial value f", mySys.neqs);
    mySys.evaluate_function(t_start, KAT::zero(), y0, f0);
    KokkosODE::Impl::initial_step_size(mySys, 1, t_start, 1e-6, 1e-3, y0, f0, temp, dt);

    std::cout << "Initial Step Size: dt=" << dt << std::endl;
  }

  KokkosODE::Experimental::BDFSolve(mySys, t_start, t_end, 0.0117188, (t_end - t_start) / 10, y0, y_new, temp, temp2);
}

template <class Device, class scalar_type>
void test_BDF_adaptive_stiff() {
  using execution_space = typename Device::execution_space;
  using vec_type        = Kokkos::View<scalar_type*, execution_space>;
  using mat_type        = Kokkos::View<scalar_type**, execution_space>;
  using KAT             = Kokkos::ArithTraits<scalar_type>;

  StiffChemistry mySys{};

  const scalar_type t_start = KAT::zero(), t_end = 350 * KAT::one();
  scalar_type dt = KAT::zero();
  vec_type y0("initial conditions", mySys.neqs), y_new("solution", mySys.neqs);

  // Set initial conditions
  auto y0_h = Kokkos::create_mirror_view(y0);
  y0_h(0)   = KAT::one();
  y0_h(1)   = KAT::zero();
  y0_h(2)   = KAT::zero();
  Kokkos::deep_copy(y0, y0_h);

  mat_type temp("buffer1", mySys.neqs, 23 + 2 * mySys.neqs + 4), temp2("buffer2", 6, 7);

  Kokkos::RangePolicy<execution_space> policy(0, 1);
  BDF_Solve_wrapper bdf_wrapper(mySys, t_start, t_end, dt, (t_end - t_start) / 10, y0, y_new, temp, temp2);

  Kokkos::parallel_for(policy, bdf_wrapper);

  auto y_new_h = Kokkos::create_mirror_view(y_new);
  Kokkos::deep_copy(y_new_h, y_new);
  std::cout << "Stiff Chemistry solution at t=500: {" << y_new_h(0) << ", " << y_new_h(1) << ", " << y_new_h(2) << "}"
            << std::endl;
}

}  // namespace Test

TEST_F(TestCategory, BDF_Logistic_serial) { ::Test::test_BDF_Logistic<TestDevice, double>(); }
TEST_F(TestCategory, BDF_LotkaVolterra_serial) { ::Test::test_BDF_LotkaVolterra<TestDevice, double>(); }
TEST_F(TestCategory, BDF_StiffChemistry_serial) { ::Test::test_BDF_StiffChemistry<TestDevice, double>(); }
// TEST_F(TestCategory, BDF_parallel_serial) {
//   ::Test::test_BDF_parallel<TestDevice, double>();
// }
TEST_F(TestCategory, BDF_Nordsieck) { ::Test::test_Nordsieck<TestDevice, double>(); }
// TEST_F(TestCategory, BDF_adaptive) {
//   ::Test::test_adaptive_BDF<TestDevice, double>();
//   ::Test::test_adaptive_BDF_v2<TestDevice, double>();
// }
TEST_F(TestCategory, BDF_StiffChemistry_adaptive) { ::Test::test_BDF_adaptive_stiff<TestDevice, double>(); }
