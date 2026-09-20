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
#include "ROL_Algorithm.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_MiniTensor_EqualityConstraint.hpp"
#include "ROL_MiniTensor_Function.hpp"
#include "ROL_NonlinearLeastSquaresObjective.hpp"
#include "ROL_StatusTest.hpp"

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

TEST(MiniTensor_ROL, Paraboloid)
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

  using MSFN = minitensor::Paraboloid<Real, DIM>;

  minitensor::Vector<Real, DIM>
  min(0.0, 0.0);

  MSFN
  msfn(0.0, 0.0);

  ROL::MiniTensor_Objective<MSFN, Real, DIM>
  obj(msfn);

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist("Line Search").sublist("Descent Method")
      .set("Type", "Newton-Krylov");

  params.sublist("Status Test").set("Gradient Tolerance", 10.e-12);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-14);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Define algorithm.
  ROL::Ptr<ROL::Step<Real>>
  step = ROL::makePtr<ROL::LineSearchStep<Real>>(params);
  ROL::Ptr<ROL::StatusTest<Real>>
  status = ROL::makePtr<ROL::StatusTest<Real>>(params);
  ROL::Algorithm<Real>
  algo(step,status,false);

  // Set Initial Guess
  minitensor::Vector<Real, DIM>
  xval(minitensor::Filler::RANDOM);

  ROL::MiniTensorVector<Real, DIM> x(xval);

  // Run Algorithm
  algo.run(x, obj, true, os);

  minitensor::Vector<Real, DIM> const
  sol = ROL::MTfromROL<Real, DIM>(x);

  os << "Solution : " << sol << '\n';

  Real const
  epsilon{minitensor::machine_epsilon<Real>()};

  Real const
  error = minitensor::norm(sol - min);

  ASSERT_LE(error, epsilon);
}

TEST(MiniTensor_ROL, Rosenbrock)
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

  using MSFN = minitensor::Rosenbrock<Real, DIM>;

  Real const
  a = 1.0;

  Real const
  b = 100.0;

  MSFN
  msfn(a, b);

  ROL::MiniTensor_Objective<MSFN, Real, DIM>
  obj(msfn);

  // Set parameters.
  Teuchos::ParameterList
  params;

  params.sublist("Step").sublist("Line Search").sublist("Descent Method")
      .set("Type", "Newton-Krylov");

  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-16);
  params.sublist("Status Test").set("Iteration Limit", 128);

  // Define algorithm.
  ROL::Ptr<ROL::Step<Real>>
  step = ROL::makePtr<ROL::LineSearchStep<Real>>(params);
  ROL::Ptr<ROL::StatusTest<Real>>
  status = ROL::makePtr<ROL::StatusTest<Real>>(params);
  ROL::Algorithm<Real>
  algo(step,status,false);

  // Set Initial Guess
  minitensor::Vector<Real, DIM>
  xval(minitensor::Filler::RANDOM);

  ROL::MiniTensorVector<Real, DIM>
  x(xval);

  // Run Algorithm
  algo.run(x, obj, true, os);

  minitensor::Vector<Real, DIM> const
  sol = ROL::MTfromROL<Real, DIM>(x);

  os << "Solution : " << sol << '\n';

  Real const
  epsilon{2.0 * minitensor::machine_epsilon<Real>()};

  xval(0) = a;
  xval(1) = a * a;

  Real const
  error = minitensor::norm(sol - xval);

  ASSERT_LE(error, epsilon);
}

// Disabled for now
#if 0
TEST(MiniTensor_ROL, NLLS01)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_CONSTR{3};

  constexpr minitensor::Index
  NUM_VAR{5};

  using MSEC = minitensor::Nonlinear01<Real, NUM_CONSTR>;

  MSEC
  msec;

  ROL::MiniTensor_EqualityConstraint<MSEC, Real, NUM_CONSTR, NUM_VAR>
  constr(msec);

  minitensor::Vector<Real, NUM_VAR>
  xval(minitensor::Filler::ZEROS);

  minitensor::Vector<Real, NUM_CONSTR>
  cval(minitensor::Filler::ZEROS);

  minitensor::Vector<Real, NUM_VAR>
  solval(minitensor::Filler::ZEROS);

  // Set initial guess.
  xval(0) = -1.8;
  xval(1) =  1.7;
  xval(2) =  1.9;
  xval(3) = -0.8;
  xval(4) = -0.8;

  // Set solution.
  solval(0) = -1.717143570394391e+00;
  solval(1) =  1.595709690183565e+00;
  solval(2) =  1.827245752927178e+00;
  solval(3) = -7.636430781841294e-01;
  solval(4) = -7.636430781841294e-01;

  Real const
  error_full_hess{2.3621708067012991e-02};

  Real const
  error_gn_hess{2.3669791103726853e-02};

  Real const
  tol{1.0e-08};

  ROL::MiniTensorVector<Real, NUM_VAR>
  x(xval);

  ROL::MiniTensorVector<Real, NUM_CONSTR>
  c(cval);

  ROL::MiniTensorVector<Real, NUM_VAR>
  sol(solval);

  ROL::Ptr<ROL::EqualityConstraint<Real>>
  pconstr = &constr, false;

  // Define algorithm.
  Teuchos::ParameterList
  params;

  std::string
  step{"Trust Region"};

  params.sublist("Step").sublist(step).set("Subproblem Solver", "Truncated CG");
  params.sublist("Status Test").set("Gradient Tolerance", 1.0e-10);
  params.sublist("Status Test").set("Constraint Tolerance", 1.0e-10);
  params.sublist("Status Test").set("Step Tolerance", 1.0e-18);
  params.sublist("Status Test").set("Iteration Limit", 128);

  ROL::Ptr<ROL::Step<Real>>
  step = ROL::makePtr<ROL::LineSearchStep<Real>>(params);
  ROL::Ptr<ROL::StatusTest<Real>>
  status = ROL::makePtr<ROL::StatusTest<Real>>(params);
  ROL::Algorithm<Real>
  algo(step,status,false);

  ROL::NonlinearLeastSquaresObjective<Real>
  nlls(pconstr, x, c, false);

  os << "\nSOLVE USING FULL HESSIAN\n";

  algo.run(x, nlls, true, os);

  minitensor::Vector<Real, NUM_VAR>
  xfinal = ROL::MTfromROL<Real, NUM_VAR>(x);

  os << "\nfinal x : " << xfinal << "\n";

  Real
  error = std::abs(minitensor::norm(xfinal - solval) - error_full_hess);

  os << "\nerror : " << error << "\n";

  ASSERT_LE(error, tol);

  algo.reset();
  x.set(xval);

  ROL::NonlinearLeastSquaresObjective<Real>
  gnnlls(pconstr, x, c, true);

  os << "\nSOLVE USING GAUSS-NEWTON HESSIAN\n";

  algo.run(x, gnnlls, true, os);

  xfinal = ROL::MTfromROL<Real, NUM_VAR>(x);

  os << "\nfinal x : " << xfinal << "\n";

  error = std::abs(minitensor::norm(xfinal - solval) - error_gn_hess);

  os << "\nerror : " << error << "\n";

  ASSERT_LE(error, tol);
}
#endif
