// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#include <gtest/gtest.h>
#include <Intrepid2_MiniTensor_FunctionSet.h>
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
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  DIM{2};

  Real const
  a = 1.0;

  Real const
  b = 100.0;

  Intrepid2::Rosenbrock<Real, DIM>
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
  Intrepid2::Vector<Real, DIM>
  x(Intrepid2::RANDOM);

  ROL::MiniTensor_Minimizer<Real, DIM>
  minimizer;

  minimizer.solve(algoname, params, fn, x);

  minimizer.printReport(os);

  Real const
  tol{1.0e-10};

  Intrepid2::Vector<Real, DIM>
  soln(a, a * a);

  Real const
  error = Intrepid2::norm(soln - x);

  ASSERT_LE(error, tol);
}

TEST(MiniTensor_ROL, Paraboloid_EqualityConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  NUM_VAR{2};

  constexpr Intrepid2::Index
  NUM_CONSTR{1};

  Real const
  a = 2.0;

  Real const
  b = 0.0;

  Real const
  r = 1.0;

  // Function to optimize
  Intrepid2::Paraboloid<Real, NUM_VAR>
  fn;

  // Constraint that defines the feasible region
  Intrepid2::Circumference<Real, NUM_CONSTR, NUM_VAR>
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
  Intrepid2::Vector<Real, NUM_VAR>
  x(Intrepid2::RANDOM);

  // Set constraint vector
  Intrepid2::Vector<Real, NUM_CONSTR>
  c(Intrepid2::RANDOM);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, eq_constr, x, c);

  minimizer.printReport(os);

  Real const
  tol{0.04};

  Intrepid2::Vector<Real, NUM_VAR>
  soln(1.0, 0.0);

  Real const
  error = Intrepid2::norm(soln - x);

  ASSERT_LE(error, tol);
}

TEST(MiniTensor_ROL, Paraboloid_BoundConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  NUM_VAR{2};

  Intrepid2::Vector<Real, NUM_VAR>
  lo(1.0, -10.0);

  Intrepid2::Vector<Real, NUM_VAR>
  hi(10.0, 10.0);

  // Function to optimize
  Intrepid2::Paraboloid<Real, NUM_VAR>
  fn;

  // Constraint that defines the feasible region
  Intrepid2::Bounds<Real, NUM_VAR>
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
  Intrepid2::Vector<Real, NUM_VAR>
  x(Intrepid2::RANDOM);

  ROL::MiniTensor_Minimizer<Real, NUM_VAR>
  minimizer;

  minimizer.solve(algoname, params, fn, bounds, x);

  minimizer.printReport(os);

  Real const
  tol{0.04};

  Intrepid2::Vector<Real, NUM_VAR>
  soln(1.0, 0.0);

  Real const
  error = Intrepid2::norm(soln - x);

  ASSERT_LE(error, tol);
}
