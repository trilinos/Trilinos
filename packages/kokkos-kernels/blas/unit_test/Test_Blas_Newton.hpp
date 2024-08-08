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

#include <KokkosBlas_Newton_impl.hpp>
#include <KokkosKernels_TestUtils.hpp>

namespace Test {

// Logistic equation
// dy/dt=y(1-y)
//
// solution y = 1/(1+exp(-t))
// y(0)=0.5
//
// Using BDF1 to integrate:
// y-y_n=dt*y*(1-y)
//
// Residual: r = y - y_n - dt*y*(1-y)
// Jacobian: J = 1 - dt + 2*dt*y
template <typename scalar_type, typename execution_space>
struct LogisticEquation {
  using vec_type = Kokkos::View<scalar_type*, execution_space>;
  using mat_type = Kokkos::View<scalar_type**, execution_space>;

  const int neqs = 1;
  scalar_type dt;
  vec_type state;

  LogisticEquation(const scalar_type dt_, vec_type initial_state) : dt(dt_), state(initial_state) {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& dydt) const {
    dydt(0) = y(0) - state(0) - dt * y(0) * (1 - y(0));
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const { jac(0, 0) = 1 - dt + 2 * dt * y(0); }

  KOKKOS_FUNCTION scalar_type expected_val(const scalar_type t) const {
    using Kokkos::exp;

    return static_cast<scalar_type>(1 / (1 + exp(-t)));
  }

  KOKKOS_FUNCTION int num_equations() const { return neqs; }
};

// Intersection of square and hyperbola
// x^2 + y^2 = 20
// x^2 - y^2 = -2
//
// solution: x = +/- 3
//           y = +/- sqrt(11)
//
// Residual: r = [x^2 + y^2 - 20]
//               [x^2 - y^2 +  2]
// Jacobian: J = [2*x,  2*y]
//               [2*x, -2*y]
template <typename scalar_type, typename execution_space>
struct Intersection {
  using vec_type = Kokkos::View<scalar_type*, execution_space>;
  using mat_type = Kokkos::View<scalar_type**, execution_space>;

  const int neqs = 2;

  Intersection() = default;

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& dydt) const {
    dydt(0) = y(0) * y(0) + y(1) * y(1) - 20;
    dydt(1) = y(0) * y(0) - y(1) * y(1) + 2;
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    jac(0, 0) = 2 * y(0);
    jac(0, 1) = 2 * y(1);
    jac(1, 0) = 2 * y(0);
    jac(1, 1) = -2 * y(1);
  }

  KOKKOS_FUNCTION int num_equations() const { return neqs; }
};

template <class solver>
struct NewtonWrapper {
  solver newton_solver;

  NewtonWrapper(solver newton_solver_) : newton_solver(newton_solver_){};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int /* system_index */) const { newton_solver.solve(); }
};

template <typename execution_space, typename scalar_type>
int test_logistic() {
  using vec_type    = typename Kokkos::View<scalar_type*, execution_space>;
  using mat_type    = typename Kokkos::View<scalar_type**, execution_space>;
  using norm_type   = typename Kokkos::View<scalar_type*, execution_space>;
  using handle_type = KokkosBlas::Impl::NewtonHandle<norm_type>;
  using system_type = LogisticEquation<scalar_type, execution_space>;
  using newton_type = KokkosBlas::Impl::NewtonFunctor<system_type, mat_type, vec_type, vec_type, handle_type>;

  // Create the non-linear system and initialize data
  vec_type state("state", 1);
  Kokkos::deep_copy(state, 0.5);
  system_type ode(0.1, state);

  vec_type x("solution vector", 1), rhs("right hand side vector", 1);
  Kokkos::deep_copy(x, 0.5);

  // Create the solver and wrapper
  handle_type handle;
  handle.debug_mode = false;
  newton_type newton_solver(ode, x, rhs, handle);
  NewtonWrapper<newton_type> wrapper(newton_solver);

  // Launch the problem in a parallel_for
  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  Kokkos::parallel_for(my_policy, wrapper);

  // Get the solution back and test it
  auto x_h = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(x_h, x);
  printf("Non-linear problem solution:\n");
  printf("  [%f]\n", x_h(0));

  return 0;
}

template <typename execution_space, typename scalar_type>
int test_intersection() {
  using vec_type    = typename Kokkos::View<scalar_type*, execution_space>;
  using mat_type    = typename Kokkos::View<scalar_type**, execution_space>;
  using norm_type   = typename Kokkos::View<scalar_type*, execution_space>;
  using handle_type = KokkosBlas::Impl::NewtonHandle<norm_type>;
  using system_type = Intersection<scalar_type, execution_space>;
  using newton_type = KokkosBlas::Impl::NewtonFunctor<system_type, mat_type, vec_type, vec_type, handle_type>;

  // Create the non-linear system and initialize data
  system_type intersection;
  vec_type x("solution vector", 2), rhs("right hand side vector", 2);
  {
    typename vec_type::HostMirror x_h = Kokkos::create_mirror_view(x);
    x_h(0)                            = 2.5;
    x_h(1)                            = 3.0;
    Kokkos::deep_copy(x, x_h);
  }

  // Create the solver and wrapper
  handle_type handle;
  handle.debug_mode = false;
  newton_type newton_solver(intersection, x, rhs, handle);
  NewtonWrapper<newton_type> wrapper(newton_solver);

  // Launch the problem in a parallel_for
  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  Kokkos::parallel_for(my_policy, wrapper);

  // Get the solution back and test it
  auto x_h = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(x_h, x);
  printf("Non-linear problem solution:\n");
  for (int idx = 0; idx < x_h.extent_int(0); ++idx) {
    printf("  [%f]\n", x_h(idx));
  }
  EXPECT_NEAR_KK(x_h(0), 3.0, 3.0e-4);
  EXPECT_NEAR_KK(x_h(1), 3.3166247903553998, 3.3166247903553998 * 1.0e-4);

  return 0;
}

}  // namespace Test

template <class scalar_type>
int test_newton() {
  Test::test_logistic<TestDevice, scalar_type>();
  Test::test_intersection<TestDevice, scalar_type>();

  return 1;
}

#if defined(KOKKOSKERNELS_INST_DOUBLE)
TEST_F(TestCategory, newton_serial) { test_newton<double>(); }
#endif
