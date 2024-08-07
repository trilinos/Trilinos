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

TEST(MiniTensor_ROL, MT_Basics)
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

  minitensor::Paraboloid<Real>
  p;

  minitensor::Vector<Real, DIM> const
  x(0.0, 0.0);

  Real const
  f = p.value(x);

  minitensor::Vector<Real, DIM> const
  df = p.gradient(x);

  minitensor::Tensor<Real, DIM> const
  ddf = p.hessian(x);

  os << "Point   : " << x << '\n';
  os << "Value   : " << f << '\n';
  os << "Gradient: " << df << '\n';
  os << "Hessian : " << ddf << '\n';

  minitensor::Tensor<Real, DIM> const
  I = minitensor::identity<Real, DIM>(DIM);

  Real const
  error = std::sqrt(f) + minitensor::norm(df) + minitensor::norm(ddf - 2.0 * I);

  ASSERT_LE(error, minitensor::machine_epsilon<Real>());
}

TEST(MiniTensor_ROL, Objective)
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

  MSFN
  msfn(0.0, 0.0);

  ROL::MiniTensor_Objective<MSFN, Real, DIM>
  obj(msfn);

  minitensor::Vector<Real, DIM>
  xval(minitensor::Filler::RANDOM);

  os << "xval:" << xval << '\n';

  minitensor::Vector<Real, DIM>
  dval(minitensor::Filler::RANDOM);

  os << "dval:" << dval << '\n';

  minitensor::Vector<Real, DIM>
  vval(minitensor::Filler::RANDOM);

  os << "vval:" << vval << '\n';

  ROL::MiniTensorVector<Real, DIM>
  x(xval);

  ROL::MiniTensorVector<Real, DIM>
  d(dval);

  ROL::MiniTensorVector<Real, DIM>
  v(dval);

  std::vector<std::vector<Real>>
  grad_check = obj.checkGradient(x, d, print_output, os);

  std::vector<std::vector<Real>>
  hess_check = obj.checkHessVec(x, v, print_output, os);

  std::vector<Real>
  symm_check = obj.checkHessSym(x, d, v, print_output, os);

  Real
  error1{1.0};

  Real
  error2{1.0};

  Real const
  error3 = symm_check[3];

  for (minitensor::Index i = 0; i < grad_check.size(); ++i) {
    error1 = std::min(error1, grad_check[i][3]);
    error2 = std::min(error2, hess_check[i][3]);
  }

  Real const
  epsilon{minitensor::machine_epsilon<Real>()};

  Real const
  tol{std::sqrt(epsilon)};

  ASSERT_LE(error1, tol);
  ASSERT_LE(error2, epsilon);
  ASSERT_LE(error3, epsilon);
}

TEST(MiniTensor_ROL, EqualityConstraintId)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  NUM_CONSTR{2};

  constexpr minitensor::Index
  NUM_VAR{2};

  using MSEC = minitensor::Identity<Real, NUM_CONSTR, NUM_VAR>;

  MSEC
  msec;

  ROL::MiniTensor_EqualityConstraint<MSEC, Real, NUM_CONSTR, NUM_VAR>
  constr(msec);

  minitensor::Vector<Real, NUM_VAR>
  xval(minitensor::Filler::RANDOM);

  os << "xval:" << xval << '\n';

  minitensor::Vector<Real, NUM_VAR>
  vval(minitensor::Filler::RANDOM);

  os << "vval:" << vval << '\n';

  minitensor::Vector<Real, NUM_CONSTR>
  jvval(minitensor::Filler::RANDOM);

  os << "jvval:" << jvval << '\n';

  minitensor::Vector<Real, NUM_VAR>
  ajvval(minitensor::Filler::RANDOM);

  os << "ajvval:" << ajvval << '\n';

  ROL::MiniTensorVector<Real, NUM_VAR>
  x(xval);

  ROL::MiniTensorVector<Real, NUM_VAR>
  v(vval);

  ROL::MiniTensorVector<Real, NUM_CONSTR>
  jv(jvval);

  ROL::MiniTensorVector<Real, NUM_VAR>
  ajv(ajvval);

  std::vector<std::vector<Real>>
  jac_check = constr.checkApplyJacobian(x, v, jv, print_output, os);

  std::vector<std::vector<Real>>
  ajac_check = constr.checkApplyAdjointJacobian(x, jv, jv, ajv, print_output, os);

  Real
  error1{1.0};

  Real
  error2{1.0};

  for (minitensor::Index i = 0; i < jac_check.size(); ++i) {
    error1 = std::min(error1, jac_check[i][3]);
    error2 = std::min(error2, jac_check[i][3]);
  }

  Real const
  tol{1.0e-07};

  ASSERT_LE(error1, tol);
  ASSERT_LE(error2, tol);
}

TEST(MiniTensor_ROL, EqualityConstraint01)
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

  using MSEC = minitensor::Nonlinear01<Real, NUM_CONSTR, NUM_VAR>;

  MSEC
  msec;

  ROL::MiniTensor_EqualityConstraint<MSEC, Real, NUM_CONSTR, NUM_VAR>
  constr(msec);

  minitensor::Vector<Real, NUM_VAR>
  xval(minitensor::Filler::RANDOM);

  os << "xval:" << xval << '\n';

  minitensor::Vector<Real, NUM_VAR>
  vval(minitensor::Filler::RANDOM);

  os << "vval:" << vval << '\n';

  minitensor::Vector<Real, NUM_CONSTR>
  jvval(minitensor::Filler::RANDOM);

  os << "jvval:" << jvval << '\n';

  minitensor::Vector<Real, NUM_VAR>
  ajvval(minitensor::Filler::RANDOM);

  os << "ajvval:" << ajvval << '\n';

  ROL::MiniTensorVector<Real, NUM_VAR>
  x(xval);

  ROL::MiniTensorVector<Real, NUM_VAR>
  v(vval);

  ROL::MiniTensorVector<Real, NUM_CONSTR>
  jv(jvval);

  ROL::MiniTensorVector<Real, NUM_VAR>
  ajv(ajvval);

  std::vector<std::vector<Real>>
  jac_check = constr.checkApplyJacobian(x, v, jv, print_output, os);

  std::vector<std::vector<Real>>
  ajac_check = constr.checkApplyAdjointJacobian(x, jv, jv, ajv, print_output, os);

  Real
  error1{1.0};

  Real
  error2{1.0};

  for (minitensor::Index i = 0; i < jac_check.size(); ++i) {
    error1 = std::min(error1, jac_check[i][3]);
    error2 = std::min(error2, ajac_check[i][3]);
  }

  Real const
  tol{1.0e-07};

  ASSERT_LE(error1, tol);
  ASSERT_LE(error2, tol);
}

TEST(MiniTensor_ROL, BoundConstraint)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr minitensor::Index
  DIM{16};

  minitensor::Vector<Real, DIM>
  lo_mt(minitensor::Filler::ONES);

  lo_mt *= -0.5;

  os << "Lower bound:" << lo_mt << '\n';

  minitensor::Vector<Real, DIM>
  hi_mt(minitensor::Filler::ONES);

  hi_mt *= 0.5;

  os << "Upper bound:" << hi_mt << '\n';

  minitensor::Vector<Real, DIM>
  x_mt(minitensor::Filler::RANDOM);

  os << "Initial x  :" << x_mt << '\n';

  ROL::MiniTensorVector<Real, DIM>
  lo_rol(lo_mt);

  ROL::MiniTensorVector<Real, DIM>
  hi_rol(hi_mt);

  ROL::MiniTensorVector<Real, DIM>
  x_rol(x_mt);

  ROL::MiniTensor_BoundConstraint<Real, DIM>
  rol_bounds(lo_rol, hi_rol);

  rol_bounds.project(x_rol);

  x_mt = ROL::MTfromROL<Real, DIM>(x_rol);

  os << "Pruned x   :" << x_mt << '\n';

  bool const
  is_feasible = rol_bounds.isFeasible(x_rol);

  ASSERT_EQ(is_feasible, true);
}
