// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include "KokkosKernels_TestUtils.hpp"

#include "KokkosODE_Newton.hpp"

namespace Test {

template <class system_type, class mat_type, class vec_type, class status_view, class scale_type>
struct NewtonSolve_wrapper {
  using newton_params = KokkosODE::Experimental::Newton_params;

  system_type my_nls;
  newton_params params;

  vec_type x, rhs, update;
  mat_type J, tmp;
  status_view status;

  scale_type scale;

  NewtonSolve_wrapper(const system_type& my_nls_, const newton_params& params_, const vec_type& x_,
                      const vec_type& rhs_, const vec_type& update_, const mat_type& J_, const mat_type& tmp_,
                      const status_view& status_, const scale_type& scale_)
      : my_nls(my_nls_),
        params(params_),
        x(x_),
        rhs(rhs_),
        update(update_),
        J(J_),
        tmp(tmp_),
        status(status_),
        scale(scale_) {}

  KOKKOS_FUNCTION
  void operator()(const int idx) const {
    // Take subviews to create the local problem
    auto local_x = Kokkos::subview(
        x, Kokkos::pair<int, int>(static_cast<int>(my_nls.neqs * idx), static_cast<int>(my_nls.neqs * (idx + 1))));
    auto local_rhs = Kokkos::subview(
        rhs, Kokkos::pair<int, int>(static_cast<int>(my_nls.neqs * idx), static_cast<int>(my_nls.neqs * (idx + 1))));
    auto local_update = Kokkos::subview(
        update, Kokkos::pair<int, int>(static_cast<int>(my_nls.neqs * idx), static_cast<int>(my_nls.neqs * (idx + 1))));
    auto local_J = Kokkos::subview(
        J, Kokkos::pair<int, int>(static_cast<int>(my_nls.neqs * idx), static_cast<int>(my_nls.neqs * (idx + 1))),
        Kokkos::ALL());
    auto local_tmp = Kokkos::subview(
        tmp, Kokkos::pair<int, int>(static_cast<int>(my_nls.neqs * idx), static_cast<int>(my_nls.neqs * (idx + 1))),
        Kokkos::ALL());

    // Run Newton nonlinear solver
    int newton_iterations;
    status(idx) = KokkosODE::Experimental::Newton::Solve(my_nls, params, local_J, local_tmp, local_x, local_rhs,
                                                         local_update, scale, newton_iterations);
  }
};

template <class system_type, class Device, class scalar_type>
void run_newton_test(const system_type& mySys, KokkosODE::Experimental::Newton_params& params,
                     const scalar_type* const initial_val, const scalar_type* const solution) {
  using execution_space      = typename Device::execution_space;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using vec_type             = typename Kokkos::View<scalar_type*, Device>;
  using mat_type             = typename Kokkos::View<scalar_type**, Device>;

  Kokkos::View<newton_solver_status*, Device> status("Newton status", 1);

  vec_type scale("scaling factors", mySys.neqs);
  Kokkos::deep_copy(scale, 1);

  vec_type x("solution vector", mySys.neqs), rhs("right hand side vector", mySys.neqs);
  auto x_h = Kokkos::create_mirror_view(x);
  auto r_h = Kokkos::create_mirror_view(rhs);

  vec_type update("update", mySys.neqs);
  mat_type J("jacobian", mySys.neqs, mySys.neqs), tmp("temp mem", mySys.neqs, mySys.neqs + 4);

  // Initial values
  for (int eqIdx = 0; eqIdx < mySys.neqs; ++eqIdx) {
    x_h(eqIdx) = initial_val[eqIdx];
  }
  Kokkos::deep_copy(x, x_h);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  NewtonSolve_wrapper solve_wrapper(mySys, params, x, rhs, update, J, tmp, status, scale);

  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto status_h = Kokkos::create_mirror_view(status);
  Kokkos::deep_copy(status_h, status);
  EXPECT_TRUE(status_h(0) == newton_solver_status::NLS_SUCCESS);

  Kokkos::deep_copy(x_h, x);
  Kokkos::deep_copy(r_h, rhs);
#ifndef NDEBUG
  std::cout << "Non-linear problem solution and residual:" << std::endl;
  std::cout << "  [(";
  for (int eqIdx = 0; eqIdx < mySys.neqs; ++eqIdx) {
    std::cout << " " << x_h(eqIdx);
    if (Kokkos::abs(solution[eqIdx]) > 0)
      std::cout << ", " << Kokkos::abs(x_h(eqIdx) - solution[eqIdx]) / Kokkos::abs(solution[eqIdx]);
    else
      std::cout << ", " << Kokkos::abs(x_h(eqIdx) - solution[eqIdx]);
  }
  std::cout << " )]" << std::endl;
#else
  (void)solution;
#endif
}

// Quadratic equation
// x^2 - x - 2 = 0
// Solution: x = 2 or x = -1
// Derivative 2*x - 1
template <typename Device, typename scalar_type>
struct QuadraticEquation {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 1;

  QuadraticEquation() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const { f(0) = y(0) * y(0) - y(0) - 2; }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const { jac(0, 0) = 2 * y(0) - 1; }
};

// Trigonometric equation
// f(x) = cos(x) - x = 0
// Solution: 0.739085
// f'(x) = -sin(x) - 1
template <typename Device, typename scalar_type>
struct TrigonometricEquation {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 1;

  TrigonometricEquation() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const { f(0) = Kokkos::cos(y(0)) - y(0); }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const { jac(0, 0) = -Kokkos::sin(y(0)) - 1; }
};

// Logarithmic equation
// f(x) = 7x - log(7x) - 1 = 0
// Solution: 1/7 = 0.14285714285
// f'(x) = 7 - (1 / x)
template <typename Device, typename scalar_type>
struct LogarithmicEquation {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 1;

  LogarithmicEquation() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const {
    f(0) = 7 * y(0) - Kokkos::log(7 * y(0)) - 1;
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const { jac(0, 0) = 7 - 1 / y(0); }
};

template <typename Device, typename scalar_type>
void test_newton_status() {
  using execution_space      = typename Device::execution_space;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using vec_type             = typename Kokkos::View<scalar_type*, Device>;
  using mat_type             = typename Kokkos::View<scalar_type**, Device>;

  vec_type scale("scaling factors", 1);
  Kokkos::deep_copy(scale, 1);

  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);
  Kokkos::View<newton_solver_status*, Device> status("newton solver status", 1);
  auto status_h = Kokkos::create_mirror_view(status);

  // Create the non-linear system and initialize data
  QuadraticEquation<Device, scalar_type> my_system{};

  scalar_type initial_value[3] = {1.0, -0.5, 0.5};
#ifndef NDEBUG
  scalar_type solution[3] = {2.0, -1.0, 0.0};
#endif
  newton_solver_status newton_status[3] = {newton_solver_status::NLS_SUCCESS, newton_solver_status::NLS_DIVERGENCE,
                                           newton_solver_status::LIN_SOLVE_FAIL};
  vec_type x("solution vector", 1), rhs("right hand side vector", 1);
  auto x_h = Kokkos::create_mirror_view(x);
  auto r_h = Kokkos::create_mirror_view(rhs);

  vec_type update("update", 1);
  mat_type J("jacobian", 1, 1), tmp("temp mem", 1, 5);

  for (int idx = 0; idx < 3; ++idx) {
    params.max_iters = (idx == 1) ? 2 : 50;
    Kokkos::deep_copy(x, initial_value[idx]);

    Kokkos::RangePolicy<execution_space> my_policy(0, 1);
    NewtonSolve_wrapper solve_wrapper(my_system, params, x, rhs, update, J, tmp, status, scale);
    Kokkos::parallel_for(my_policy, solve_wrapper);

    Kokkos::deep_copy(status_h, status);
    EXPECT_TRUE(status_h(0) == newton_status[idx]);

#ifndef NDEBUG
    Kokkos::deep_copy(x_h, x);
    Kokkos::deep_copy(r_h, rhs);
    printf("Non-linear problem solution and residual with initial value %f:\n", initial_value[idx]);
    printf("  [%f, %g, %g]\n", x_h(0), r_h(0), Kokkos::abs(x_h(0) - solution[idx]) / Kokkos::abs(solution[idx]));
#endif
  }
}

template <typename Device, typename scalar_type>
void test_simple_problems() {
  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);

  {
    // Test the Newton solver on a quadratci equation
    // with two different initial guess that lead to
    // the two solutions of the equation.
#ifndef NDEBUG
    std::cout << "\nStarting Quadratic Equation problem" << std::endl;
#endif
    using system_type = QuadraticEquation<Device, scalar_type>;
    system_type mySys{};
    scalar_type initial_value[2] = {1.0, -0.5}, solution[2] = {2.0, -1.0};
    for (int idx = 0; idx < 2; ++idx) {
      run_newton_test<system_type, Device, scalar_type>(mySys, params, &(initial_value[idx]), &(solution[idx]));
    }
#ifndef NDEBUG
    std::cout << "Finished Quadratic Equation problem" << std::endl;
#endif
  }

  {
    // Test the Newton solver on a trigonometric equation
#ifndef NDEBUG
    std::cout << "\nStarting Trigonometric Equation problem" << std::endl;
#endif
    using system_type = TrigonometricEquation<Device, scalar_type>;
    system_type mySys{};
    scalar_type initial_value[1] = {0.1}, solution[1] = {0.739085};
    run_newton_test<system_type, Device, scalar_type>(mySys, params, initial_value, solution);
#ifndef NDEBUG
    std::cout << "Finished Trigonometric Equation problem" << std::endl;
#endif
  }

  {
    // Test the Newton solver on a logarithmic equation
#ifndef NDEBUG
    std::cout << "\nStarting Logarithmic Equation problem" << std::endl;
#endif
    using system_type = LogarithmicEquation<Device, scalar_type>;
    system_type mySys{};
    scalar_type initial_value[1] = {static_cast<scalar_type>(0.5)},
                solution[1]      = {static_cast<scalar_type>(1.0) / static_cast<scalar_type>(7.0)};
    run_newton_test<system_type, Device, scalar_type>(mySys, params, initial_value, solution);
#ifndef NDEBUG
    std::cout << "Finished Logarithmic Equation problem" << std::endl;
#endif
  }
}

///////////////////////////////////////
// Now solving systems of equations  //
// To make things more realistic and //
// interesting.                      //
///////////////////////////////////////

// Intersections of two circles
// Equations:  f0 = (x-0)**2 + (y-0)**2 - 4.00 = 0
//             f1 = (x-3)**2 + (y-0)**2 - 2.25 = 0
//
// Jacobian:   J00 = 2*x       J01 = 2*y
//             J10 = 2*(x-3)   J11 = 2*y
//
// Solution:   x = 10.75/6     y = +/- sqrt(2.25 + 7.25/6)
//               ~ 1.7916666     ~ +/- 0.8887803753
template <typename Device, typename scalar_type>
struct CirclesIntersections {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 2;

  CirclesIntersections() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const {
    f(0) = y(0) * y(0) + y(1) * y(1) - 4;
    f(1) = (y(0) - 3) * (y(0) - 3) + y(1) * y(1) - 2.25;
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    jac(0, 0) = 2 * y(0);
    jac(0, 1) = 2 * y(1);
    jac(1, 0) = 2 * (y(0) - 3);
    jac(1, 1) = 2 * y(1);
  }
};

// Intersections of a circle and an hyperbola
// Equations:  f0 = x**2 + y**2 - 4.00 = 0
//             f1 = x*y - 1 = 0             --> also y = 1 / x
//
// Jacobian:   J00 = 2*x       J01 = 2*y
//             J10 = y         J11 = x
//
// Solution:   x = +/- sqrt( (4 +/- sqrt(12)) / 2); y = 1 / x
//             x0~  1.9318516525   y0~  0.5176380902
//             x1~  0.5176380902   y1~  1.9318516525
//             x2~ -0.5176380902   y2~ -1.9318516525
//             x3~ -1.9318516525   y3~ -0.5176380902
template <typename Device, typename scalar_type>
struct CircleHyperbolaIntersection {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 2;

  CircleHyperbolaIntersection() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const {
    f(0) = y(0) * y(0) + y(1) * y(1) - 4;
    f(1) = y(0) * y(1) - 1;
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    jac(0, 0) = 2 * y(0);
    jac(0, 1) = 2 * y(1);
    jac(1, 0) = y(1);
    jac(1, 1) = y(0);
  }
};

template <typename Device, typename scalar_type>
void test_simple_systems() {
  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);

  {
    // First problem: intersection of two circles
#ifndef NDEBUG
    std::cout << "\nStarting Circles Intersetcion problem" << std::endl;
#endif
    using system_type = CirclesIntersections<Device, scalar_type>;
    system_type mySys{};
    scalar_type initial_values[2] = {1.5, 1.5};
    scalar_type solution[2]       = {10.75 / 6, 0.8887803753};
    run_newton_test<system_type, Device, scalar_type>(mySys, params, initial_values, solution);
#ifndef NDEBUG
    std::cout << "Finished Circles Intersetcion problem" << std::endl;
#endif
  }

  {
    // Second problem: circle / hyperbola intersection
#ifndef NDEBUG
    std::cout << "\nStarting Circle/Hyperbola Intersetcion problem" << std::endl;
#endif
    using system_type = CircleHyperbolaIntersection<Device, scalar_type>;
    system_type mySys{};

    scalar_type init_vals[2] = {0.0, 1.0};
    scalar_type solutions[2] = {
        KokkosKernels::ArithTraits<scalar_type>::one() /
            Kokkos::sqrt(static_cast<scalar_type>(4 + Kokkos::sqrt(static_cast<scalar_type>(12.0)) / 2)),
        Kokkos::sqrt(static_cast<scalar_type>((4 + Kokkos::sqrt(static_cast<scalar_type>(12.0))) / 2))};
    run_newton_test<system_type, Device, scalar_type>(mySys, params, init_vals, solutions);
#ifndef NDEBUG
    std::cout << "Finished Circle/Hyperbola Intersetcion problem" << std::endl;
#endif
  }
}

// Equation whose Jacobian becomes singular at the first Newton iterate, used
// to guard against a regression where the update was applied to y0 before
// checking the linear solver status. The residual is linear, f = x - 5, but
// the Jacobian deliberately reports the wrong slope (2 instead of 1) for
// x < 0 so that the first iteration moves the iterate from -1 to 2 without
// converging, and reports 0 for x >= 0 so that the second iteration's linear
// solve fails. The solver must return LIN_SOLVE_FAIL and leave the iterate
// exactly where the last successful iteration put it. Before the fix, the
// stale update from iteration 0 was re-applied to y0 (corrupting it from 2
// back to -1) and the divergence check fired first, masking the linear solve
// failure as NLS_DIVERGENCE.
template <typename Device, typename scalar_type>
struct SingularJacobianEquation {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 1;

  SingularJacobianEquation() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const { f(0) = y(0) - 5; }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const { jac(0, 0) = (y(0) < 0) ? 2 : 0; }
};

template <typename Device, typename scalar_type>
void test_newton_lin_solve_fail_keeps_iterate() {
  using execution_space      = typename Device::execution_space;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using vec_type             = typename Kokkos::View<scalar_type*, Device>;
  using mat_type             = typename Kokkos::View<scalar_type**, Device>;
  using system_type          = SingularJacobianEquation<Device, scalar_type>;

  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);

  system_type my_system{};

  Kokkos::View<newton_solver_status*, Device> status("Newton status", 1);

  vec_type scale("scaling factors", my_system.neqs);
  Kokkos::deep_copy(scale, 1);

  vec_type x("solution vector", my_system.neqs), rhs("right hand side vector", my_system.neqs);
  vec_type update("update", my_system.neqs);
  mat_type J("jacobian", my_system.neqs, my_system.neqs), tmp("temp mem", my_system.neqs, my_system.neqs + 4);

  // Iteration 0: f(-1) = -6, jac = 2, so the iterate moves to -1 + 6/2 = 2.
  // Iteration 1: f(2) = -3 != 0 but jac = 0, so the linear solve must fail.
  Kokkos::deep_copy(x, -1);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  NewtonSolve_wrapper solve_wrapper(my_system, params, x, rhs, update, J, tmp, status, scale);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto status_h = Kokkos::create_mirror_view(status);
  Kokkos::deep_copy(status_h, status);
  auto x_h = Kokkos::create_mirror_view(x);
  Kokkos::deep_copy(x_h, x);

  // The failed linear solve must be reported as such, not misclassified by the
  // divergence check that used to run first.
  EXPECT_TRUE(status_h(0) == newton_solver_status::LIN_SOLVE_FAIL);

  // y0 must be left at the last successful iterate (2), not corrupted by
  // re-applying the stale update from iteration 0 (which brought it back to -1).
  EXPECT_EQ(x_h(0), static_cast<scalar_type>(2));
}

////////////////////////////////////////////
// Finally, solving systems of equations  //
// within a parallel_for loop as it would //
// happen within a FE/FD code.            //
////////////////////////////////////////////

template <class Device, class scalar_type>
void test_newton_on_device() {
  using execution_space      = typename Device::execution_space;
  using vec_type             = Kokkos::View<scalar_type*, Device>;
  using mat_type             = Kokkos::View<scalar_type**, Device>;
  using newton_params        = KokkosODE::Experimental::Newton_params;
  using system_type          = CircleHyperbolaIntersection<Device, scalar_type>;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;

  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }

  constexpr int num_systems = 1000;
  const newton_params params(50, abs_tol, rel_tol);

  system_type mySys{};

  vec_type scale("scaling factors", mySys.neqs);
  Kokkos::deep_copy(scale, 1);

  vec_type x("solution vector", mySys.neqs * num_systems);
  vec_type rhs("right hand side vector", mySys.neqs * num_systems);
  vec_type update("update", mySys.neqs * num_systems);
  mat_type J("jacobian", mySys.neqs * num_systems, mySys.neqs);
  mat_type tmp("temp mem", mySys.neqs * num_systems, mySys.neqs + 4);

  Kokkos::View<newton_solver_status*, execution_space> status("solver status", num_systems);

  auto x_h = Kokkos::create_mirror_view(x);
  auto r_h = Kokkos::create_mirror_view(rhs);

  // Initial values
  scalar_type initial_val[2] = {0.0, 1.0};
  for (int sysIdx = 0; sysIdx < num_systems; ++sysIdx) {
    x_h(2 * sysIdx)     = initial_val[0];
    x_h(2 * sysIdx + 1) = initial_val[1];
  }
  Kokkos::deep_copy(x, x_h);

  Kokkos::RangePolicy<execution_space> my_policy(0, num_systems);
  NewtonSolve_wrapper solve_wrapper(mySys, params, x, rhs, update, J, tmp, status, scale);

  Kokkos::parallel_for(my_policy, solve_wrapper);
  Kokkos::fence();

  auto status_h = Kokkos::create_mirror_view(status);
  Kokkos::deep_copy(status_h, status);
  Kokkos::deep_copy(x_h, x);
  for (int sysIdx = 0; sysIdx < num_systems; ++sysIdx) {
    EXPECT_TRUE(status_h(sysIdx) == newton_solver_status::NLS_SUCCESS)
        << "System " << sysIdx << " did not report a successful convergence!";
  }
}

// Two decoupled equations used to guard against a regression where the norm of
// the scaled update was overwritten instead of accumulated across equations
// (and not reset between iterations). The last equation is linear and is solved
// exactly in a single Newton step, so from the second iteration on its update is
// zero. If only that last equation feeds the scaled-update norm, the norm
// becomes zero on the second iteration and the solver falsely reports
// convergence -- even though the nonlinear first component, started far from its
// root, is still nowhere near converged.
// Equations:  f0 = x**2 - 2 = 0   (solution: sqrt(2) ~ 1.41421356)
//             f1 = y - 5   = 0   (solution: 5, reached in one step)
//
// Jacobian:   J00 = 2*x    J01 = 0
//             J10 = 0      J11 = 1
template <typename Device, typename scalar_type>
struct DecoupledEquations {
  using vec_type = Kokkos::View<scalar_type*, Device>;
  using mat_type = Kokkos::View<scalar_type**, Device>;

  static constexpr int neqs = 2;

  DecoupledEquations() {}

  KOKKOS_FUNCTION void residual(const vec_type& y, const vec_type& f) const {
    f(0) = y(0) * y(0) - 2;
    f(1) = y(1) - 5;
  }

  KOKKOS_FUNCTION void jacobian(const vec_type& y, const mat_type& jac) const {
    jac(0, 0) = 2 * y(0);
    jac(0, 1) = 0;
    jac(1, 0) = 0;
    jac(1, 1) = 1;
  }
};

template <typename Device, typename scalar_type>
void test_newton_norm_accumulation() {
  using execution_space      = typename Device::execution_space;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using vec_type             = typename Kokkos::View<scalar_type*, Device>;
  using mat_type             = typename Kokkos::View<scalar_type**, Device>;
  using system_type          = DecoupledEquations<Device, scalar_type>;

  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);

  system_type my_system{};

  Kokkos::View<newton_solver_status*, Device> status("Newton status", 1);

  vec_type scale("scaling factors", my_system.neqs);
  Kokkos::deep_copy(scale, 1);

  vec_type x("solution vector", my_system.neqs), rhs("right hand side vector", my_system.neqs);
  vec_type update("update", my_system.neqs);
  mat_type J("jacobian", my_system.neqs, my_system.neqs), tmp("temp mem", my_system.neqs, my_system.neqs + 4);

  // Start the nonlinear component far from its root so it needs several Newton
  // iterations, while the linear component converges in one step.
  auto x_h = Kokkos::create_mirror_view(x);
  x_h(0)   = 100.0;
  x_h(1)   = 0.0;
  Kokkos::deep_copy(x, x_h);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  NewtonSolve_wrapper solve_wrapper(my_system, params, x, rhs, update, J, tmp, status, scale);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto status_h = Kokkos::create_mirror_view(status);
  Kokkos::deep_copy(status_h, status);
  Kokkos::deep_copy(x_h, x);

  EXPECT_TRUE(status_h(0) == newton_solver_status::NLS_SUCCESS);

  // The nonlinear first component must actually be converged. With the norm bug
  // the scaled-update norm sees only the already-solved linear equation, hits
  // zero on the second iteration, and the solver reports success with x(0) still
  // far from sqrt(2) (~50 rather than ~1.414).
  const scalar_type solution[2] = {static_cast<scalar_type>(Kokkos::sqrt(static_cast<scalar_type>(2.0))),
                                   static_cast<scalar_type>(5.0)};
  const scalar_type rel_check   = static_cast<scalar_type>(1e-2);
  EXPECT_NEAR(x_h(0), solution[0], rel_check * solution[0]);
  EXPECT_NEAR(x_h(1), solution[1], rel_check * solution[1]);
}

////////////////////////////////////////////
// Check that Newton::Solve reports the   //
// number of iterations it performed      //
////////////////////////////////////////////

template <class system_type, class mat_type, class vec_type, class status_view, class iters_view, class scale_type>
struct NewtonSolve_iters_wrapper {
  using newton_params = KokkosODE::Experimental::Newton_params;

  system_type my_nls;
  newton_params params;

  vec_type x, rhs, update;
  mat_type J, tmp;
  status_view status;
  iters_view iters;

  scale_type scale;

  NewtonSolve_iters_wrapper(const system_type& my_nls_, const newton_params& params_, const vec_type& x_,
                            const vec_type& rhs_, const vec_type& update_, const mat_type& J_, const mat_type& tmp_,
                            const status_view& status_, const iters_view& iters_, const scale_type& scale_)
      : my_nls(my_nls_),
        params(params_),
        x(x_),
        rhs(rhs_),
        update(update_),
        J(J_),
        tmp(tmp_),
        status(status_),
        iters(iters_),
        scale(scale_) {}

  KOKKOS_FUNCTION
  void operator()(const int idx) const {
    int newton_iterations;
    status(idx) =
        KokkosODE::Experimental::Newton::Solve(my_nls, params, J, tmp, x, rhs, update, scale, newton_iterations);
    iters(idx) = newton_iterations;
  }
};

// Solve a problem that requires multiple Newton iterations
// and check that the iteration count performed by the solver
// is reported in `newton_iterations`. If Solve runs on an internal
// copy of params instead, the count seen by the caller stays
// at zero; the BDF integrator relies on this count to compute
// its step size safety factor.
template <class Device, class scalar_type>
void test_newton_iteration_count() {
  using execution_space      = typename Device::execution_space;
  using newton_solver_status = KokkosODE::Experimental::newton_solver_status;
  using vec_type             = typename Kokkos::View<scalar_type*, Device>;
  using mat_type             = typename Kokkos::View<scalar_type**, Device>;
  using system_type          = QuadraticEquation<Device, scalar_type>;

  double abs_tol, rel_tol;
  if (std::is_same_v<scalar_type, float>) {
    rel_tol = 10e-5;
    abs_tol = 10e-7;
  } else if (std::is_same_v<scalar_type, double>) {
    rel_tol = 10e-8;
    abs_tol = 10e-15;
  } else {
    throw std::runtime_error("scalar_type is neither float, nor double!");
  }
  KokkosODE::Experimental::Newton_params params(50, abs_tol, rel_tol);

  system_type mySys{};

  vec_type scale("scaling factors", mySys.neqs);
  Kokkos::deep_copy(scale, 1);

  vec_type x("solution vector", mySys.neqs), rhs("right hand side vector", mySys.neqs);
  vec_type update("update", mySys.neqs);
  mat_type J("jacobian", mySys.neqs, mySys.neqs), tmp("temp mem", mySys.neqs, mySys.neqs + 4);

  Kokkos::View<newton_solver_status*, Device> status("Newton status", 1);
  Kokkos::View<int*, Device> iters("Newton iteration count", 1);

  // Initial guess 1.0 converges to the root x=2
  // after a few Newton iterations.
  Kokkos::deep_copy(x, 1.0);

  Kokkos::RangePolicy<execution_space> my_policy(0, 1);
  NewtonSolve_iters_wrapper solve_wrapper(mySys, params, x, rhs, update, J, tmp, status, iters, scale);
  Kokkos::parallel_for(my_policy, solve_wrapper);

  auto status_h = Kokkos::create_mirror_view(status);
  Kokkos::deep_copy(status_h, status);
  auto iters_h = Kokkos::create_mirror_view(iters);
  Kokkos::deep_copy(iters_h, iters);

  EXPECT_TRUE(status_h(0) == newton_solver_status::NLS_SUCCESS);
  EXPECT_LT(1, iters_h(0)) << "Newton::Solve did not report the number of iterations performed in `newton_iterations`";
  EXPECT_LE(iters_h(0), params.max_iters);
}

}  // namespace Test

// No ETI is performed for these device routines
// Just pick scalar types at will...
TEST_F(TestCategory, Newton_status_float) { ::Test::test_newton_status<TestDevice, float>(); }
TEST_F(TestCategory, Newton_status_double) { ::Test::test_newton_status<TestDevice, double>(); }

TEST_F(TestCategory, Newton_simple_float) { ::Test::test_simple_problems<TestDevice, float>(); }
TEST_F(TestCategory, Newton_simple_double) { ::Test::test_simple_problems<TestDevice, double>(); }

TEST_F(TestCategory, Newton_system_float) { ::Test::test_simple_systems<TestDevice, float>(); }
TEST_F(TestCategory, Newton_system_double) { ::Test::test_simple_systems<TestDevice, double>(); }

TEST_F(TestCategory, Newton_lin_solve_fail_keeps_iterate_float) {
  ::Test::test_newton_lin_solve_fail_keeps_iterate<TestDevice, float>();
}
TEST_F(TestCategory, Newton_lin_solve_fail_keeps_iterate_double) {
  ::Test::test_newton_lin_solve_fail_keeps_iterate<TestDevice, double>();
}

TEST_F(TestCategory, Newton_parallel_float) { ::Test::test_newton_on_device<TestDevice, float>(); }
TEST_F(TestCategory, Newton_parallel_double) { ::Test::test_newton_on_device<TestDevice, double>(); }

TEST_F(TestCategory, Newton_norm_accumulation_float) { ::Test::test_newton_norm_accumulation<TestDevice, float>(); }
TEST_F(TestCategory, Newton_norm_accumulation_double) { ::Test::test_newton_norm_accumulation<TestDevice, double>(); }

TEST_F(TestCategory, Newton_iteration_count_float) { ::Test::test_newton_iteration_count<TestDevice, float>(); }
TEST_F(TestCategory, Newton_iteration_count_double) { ::Test::test_newton_iteration_count<TestDevice, double>(); }
