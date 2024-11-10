// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <gtest/gtest.h>
#include <MiniTensor_FunctionSet.h>
#include "ROL_MiniTensor_BoundConstraint.hpp"
#include "ROL_MiniTensor_EqualityConstraint.hpp"
#include "ROL_MiniTensor_Function.hpp"
#include "ROL_MiniTensor_MiniSolver.hpp"

using Real = double;

int
main(int ac, char * av[])
{
  Kokkos::initialize();

  ::testing::GTEST_FLAG(print_time) = (ac > 1) ? true : false;

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

TEST(MiniTensor_ROL, Rosenbrock_Unconstrained)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  DIM{2};

  Real const
  a = 1.0;

  Real const
  b = 100.0;

  minitensor::Rosenbrock<Real, DIM>
  fn(a, b);

  // Define algorithm.
  std::string const
  algoname{"Line Search"};

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist("Line Search").sublist("Descent Method").
    set("Type", "Newton-Krylov");

  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Set initial guess
  minitensor::Vector<Real, DIM>
  x(minitensor::Filler::RANDOM);

  ROL::MiniTensor_Minimizer<Real, DIM>
  minimizer;

  minimizer.solve(algoname, params, fn, x);

  minimizer.printReport(os);

  Real const
  tol{1.0e-10};

  minitensor::Vector<Real, DIM>
  soln(a, a * a);

  Real const
  error = minitensor::norm(soln - x);

  ASSERT_LE(error, tol);
}

TEST(MiniTensor_ROL, Paraboloid_BoundConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_VAR{2};

  minitensor::Vector<Real, NUM_VAR>
  lo(1.0, -10.0);

  minitensor::Vector<Real, NUM_VAR>
  hi(10.0, 10.0);

  // Function to optimize
  minitensor::Paraboloid<Real, NUM_VAR>
  fn;

  // Constraint that defines the feasible region
  minitensor::Bounds<Real, NUM_VAR>
  bounds(lo, hi);

  // Define algorithm.
  std::string const
  algoname{"Line Search"};

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist("Line Search").sublist("Descent Method").
    set("Type", "Newton-Krylov");

  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Set initial guess
  minitensor::Vector<Real, NUM_VAR>
  x(minitensor::Filler::RANDOM);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, bounds, x);

  minimizer.printReport(os);

  Real const
  tol{0.04};

  minitensor::Vector<Real, NUM_VAR>
  soln(1.0, 0.0);

  Real const
  error = minitensor::norm(soln - x);

  ASSERT_LE(error, tol);
}

TEST(MiniTensor_ROL, Paraboloid_EqualityConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_VAR{2};

  constexpr minitensor::Index
  NUM_CONSTR{1};

  Real const
  a = 2.0;

  Real const
  b = 0.0;

  Real const
  r = 1.0;

  // Function to optimize
  minitensor::Paraboloid<Real, NUM_VAR>
  fn;

  // Constraint that defines the feasible region
  minitensor::Circumference<Real, NUM_CONSTR, NUM_VAR>
  eq_constr(r, a, b);

  // Define algorithm.
  std::string const
  algoname{"Composite Step"};

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist(algoname).
      sublist("Optimality System Solver").set(
      "Nominal Relative Tolerance",
      1.e-8);

  params.sublist("Step").sublist(algoname).
      sublist("Optimality System Solver").set("Fix Tolerance", true);

  params.sublist("Step").sublist(algoname).
      sublist("Tangential Subproblem Solver").set("Iteration Limit", 128);

  params.sublist("Step").sublist(algoname).
      sublist("Tangential Subproblem Solver").set("Relative Tolerance", 1e-6);

  params.sublist("Step").sublist(algoname).set("Output Level", 0);
  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-12);
  params.sublist("Status Test").set("Constraint Tolerance", 1.0e-12);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-18);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Set initial guess
  minitensor::Vector<Real, NUM_VAR>
  x(minitensor::Filler::RANDOM);

  // Set constraint vector
  minitensor::Vector<Real, NUM_CONSTR>
  c(minitensor::Filler::RANDOM);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, eq_constr, x, c);

  minimizer.printReport(os);

  Real const
  tol{0.04};

  minitensor::Vector<Real, NUM_VAR>
  soln(1.0, 0.0);

  Real const
  error = minitensor::norm(soln - x);

  ASSERT_LE(error, tol);
}

TEST(MiniTensor_ROL, Paraboloid_InequalityConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_VAR{2};

  constexpr minitensor::Index
  NUM_CONSTR{1};

  Real const
  a = 2.0;

  Real const
  b = 0.0;

  Real const
  r = 1.0;

  // Function to optimize
  minitensor::Paraboloid<Real, NUM_VAR>
  fn;

  // Constraint that defines the feasible region
  minitensor::Circle<Real, NUM_CONSTR, NUM_VAR>
  ineq_constr(r, a, b);

  // Define algorithm.
  std::string const
  algoname{"Composite Step"};

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist(algoname).
      sublist("Optimality System Solver").set(
      "Nominal Relative Tolerance",
      1.e-8);

  params.sublist("Step").sublist(algoname).
      sublist("Optimality System Solver").set("Fix Tolerance", true);

  params.sublist("Step").sublist(algoname).
      sublist("Tangential Subproblem Solver").set("Iteration Limit", 128);

  params.sublist("Step").sublist(algoname).
      sublist("Tangential Subproblem Solver").set("Relative Tolerance", 1e-6);

  params.sublist("Step").sublist(algoname).set("Output Level", 0);
  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-12);
  params.sublist("Status Test").set("Constraint Tolerance", 1.0e-12);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-18);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Set initial guess
  minitensor::Vector<Real, NUM_VAR>
  x(minitensor::Filler::RANDOM);

  // Set constraint vector
  minitensor::Vector<Real, NUM_CONSTR>
  c(minitensor::Filler::RANDOM);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, ineq_constr, x, c);

  minimizer.printReport(os);

  Real const
  tol{0.04};

  minitensor::Vector<Real, NUM_VAR>
  soln(1.0, 0.0);

  Real const
  error = minitensor::norm(soln - x);

  ASSERT_LE(error, tol);
}

namespace
{

template<typename S, minitensor::Index M = 2>
class HS24 : public minitensor::Function_Base<HS24<S, M>, S, M>
{
public:

  HS24()
  {
  }

  static constexpr
  char const * const
  NAME{"HS24's Function"};

  using Base = minitensor::Function_Base<HS24<S, M>, S, M>;

  // Explicit value.
  template<typename T, minitensor::Index N>
  T
  value(minitensor::Vector<T, N> const & x)
  {
    T const
    sqrt3 = std::sqrt(3.0);

    T const
    a = x(1) * x(1) * x(1);

    return sqrt3 * x(0) * a * (x(0) - 6.0) / 81.0;
  }

  // Default AD gradient.
  template<typename T, minitensor::Index N>
  minitensor::Vector<T, N>
  gradient(minitensor::Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, minitensor::Index N>
  minitensor::Tensor<T, N>
  hessian(minitensor::Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

template<typename S, minitensor::Index M = 2>
class HS4 : public minitensor::Function_Base<HS4<S, M>, S, M>
{
public:

  HS4()
  {
  }

  static constexpr
  char const * const
  NAME{"HS4's Function"};

  using Base = minitensor::Function_Base<HS4<S, M>, S, M>;

  // Explicit value.
  template<typename T, minitensor::Index N>
  T
  value(minitensor::Vector<T, N> const & x)
  {
    T const
    a = x(0) + 1.0;

    return a * a * a / 3.0 + x(1);
  }

  // Default AD gradient.
  template<typename T, minitensor::Index N>
  minitensor::Vector<T, N>
  gradient(minitensor::Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, minitensor::Index N>
  minitensor::Tensor<T, N>
  hessian(minitensor::Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// HS24 feasible region
//
template<typename S, minitensor::Index NC = 3, minitensor::Index NV = 2>
class HS24_Region : public minitensor::Inequality_Constraint<HS24_Region<S, NC, NV>, S, NC, NV>
{
public:

  HS24_Region()
  {
  }

  static constexpr
  char const * const
  NAME{"HS24 feasible region"};

  using Base = minitensor::Inequality_Constraint<HS24_Region<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  template<typename T, minitensor::Index N = 2>
  minitensor::Vector<T, NC>
  value(minitensor::Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    T const
    c = std::sqrt(3.0);

    minitensor::Vector<T, NC>
    f(minitensor::Filler::ZEROS);

    f(0) = x(0) / c- x(1);

    f(1) = x(0) + c * x(1);

    f(2) = -x(0) - c * x(1) + 6.0;

    return f;
  }

  // Default AD gradient.
  template<typename T, minitensor::Index N = 2>
  minitensor::Matrix<T, NC, NV>
  gradient(minitensor::Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

};

} // anonymous namespace

TEST(MiniTensor_ROL, HS24_BoundOnlyMod)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_VAR{2};

  // Function to optimize
  HS24<Real, NUM_VAR>
  fn;

  // Define algorithm.
  std::string const
  algoname{"Line Search"};

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist("Line Search").sublist("Descent Method").
    set("Type", "Newton-Krylov");

  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // These are not the original bounds for these problems.
  // We use them to check the algorithm.
  minitensor::Vector<Real, NUM_VAR>
  lo(-1.0, -1.0);

  minitensor::Vector<Real, NUM_VAR>
  hi(1.0, 1.0);

  // Constraint that defines the feasible region
  minitensor::Bounds<Real, NUM_VAR>
  bounds(lo, hi);

  // Set initial guess
  minitensor::Vector<Real, NUM_VAR>
  x(-0.9, -0.9);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, bounds, x);

  minimizer.printReport(os);

  Real const
  tol{1.0e-10};

  minitensor::Vector<Real, NUM_VAR>
  soln(-1.0, -1.0);

  Real const
  error = minitensor::norm(soln - x);

  ASSERT_LE(error, tol);
}
