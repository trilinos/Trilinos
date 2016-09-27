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
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  DIM{2};

  Intrepid2::Paraboloid<Real>
  p;

  Intrepid2::Vector<Real, DIM> const
  x(0.0, 0.0);

  Real const
  f = p.value(x);

  Intrepid2::Vector<Real, DIM> const
  df = p.gradient(x);

  Intrepid2::Tensor<Real, DIM> const
  ddf = p.hessian(x);

  os << "Point   : " << x << '\n';
  os << "Value   : " << f << '\n';
  os << "Gradient: " << df << '\n';
  os << "Hessian : " << ddf << '\n';

  Intrepid2::Tensor<Real, DIM> const
  I = Intrepid2::identity<Real, DIM>(DIM);

  Real const
  error = std::sqrt(f) + Intrepid2::norm(df) + Intrepid2::norm(ddf - 2.0 * I);

  ASSERT_LE(error, Intrepid2::machine_epsilon<Real>());
}

TEST(MiniTensor_ROL, Objective)
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

  using MSFN = Intrepid2::Paraboloid<Real, DIM>;

  MSFN
  msfn(0.0, 0.0);

  ROL::MiniTensor_Objective<MSFN, Real, DIM>
  obj(msfn);

  Intrepid2::Vector<Real, DIM>
  xval(Intrepid2::RANDOM);

  os << "xval:" << xval << '\n';

  Intrepid2::Vector<Real, DIM>
  dval(Intrepid2::RANDOM);

  os << "dval:" << dval << '\n';

  Intrepid2::Vector<Real, DIM>
  vval(Intrepid2::RANDOM);

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

  for (Intrepid2::Index i = 0; i < grad_check.size(); ++i) {
    error1 = std::min(error1, grad_check[i][3]);
    error2 = std::min(error2, hess_check[i][3]);
  }

  Real const
  epsilon{Intrepid2::machine_epsilon<Real>()};

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
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  ROWS{2};

  constexpr Intrepid2::Index
  COLS{2};

  using MSEC = Intrepid2::Identity<Real, ROWS>;

  MSEC
  msec;

  ROL::MiniTensor_EqualityConstraint<MSEC, Real, ROWS, COLS>
  constr(msec);

  Intrepid2::Vector<Real, COLS>
  xval(Intrepid2::RANDOM);

  os << "xval:" << xval << '\n';

  Intrepid2::Vector<Real, COLS>
  vval(Intrepid2::RANDOM);

  os << "vval:" << vval << '\n';

  Intrepid2::Vector<Real, ROWS>
  jvval(Intrepid2::RANDOM);

  os << "jvval:" << jvval << '\n';

  Intrepid2::Vector<Real, COLS>
  ajvval(Intrepid2::RANDOM);

  os << "ajvval:" << ajvval << '\n';

  ROL::MiniTensorVector<Real, COLS>
  x(xval);

  ROL::MiniTensorVector<Real, COLS>
  v(vval);

  ROL::MiniTensorVector<Real, ROWS>
  jv(jvval);

  ROL::MiniTensorVector<Real, COLS>
  ajv(ajvval);

  std::vector<std::vector<Real>>
  jac_check = constr.checkApplyJacobian(x, v, jv, print_output, os);

  std::vector<std::vector<Real>>
  ajac_check = constr.checkApplyAdjointJacobian(x, jv, jv, ajv, print_output, os);

  Real
  error1{1.0};

  Real
  error2{1.0};

  for (Intrepid2::Index i = 0; i < jac_check.size(); ++i) {
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
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  ROWS{3};

  constexpr Intrepid2::Index
  COLS{5};

  using MSEC = Intrepid2::Nonlinear01<Real, ROWS>;

  MSEC
  msec;

  ROL::MiniTensor_EqualityConstraint<MSEC, Real, ROWS, COLS>
  constr(msec);

  Intrepid2::Vector<Real, COLS>
  xval(Intrepid2::RANDOM);

  os << "xval:" << xval << '\n';

  Intrepid2::Vector<Real, COLS>
  vval(Intrepid2::RANDOM);

  os << "vval:" << vval << '\n';

  Intrepid2::Vector<Real, ROWS>
  jvval(Intrepid2::RANDOM);

  os << "jvval:" << jvval << '\n';

  Intrepid2::Vector<Real, COLS>
  ajvval(Intrepid2::RANDOM);

  os << "ajvval:" << ajvval << '\n';

  ROL::MiniTensorVector<Real, COLS>
  x(xval);

  ROL::MiniTensorVector<Real, COLS>
  v(vval);

  ROL::MiniTensorVector<Real, ROWS>
  jv(jvval);

  ROL::MiniTensorVector<Real, COLS>
  ajv(ajvval);

  std::vector<std::vector<Real>>
  jac_check = constr.checkApplyJacobian(x, v, jv, print_output, os);

  std::vector<std::vector<Real>>
  ajac_check = constr.checkApplyAdjointJacobian(x, jv, jv, ajv, print_output, os);

  Real
  error1{1.0};

  Real
  error2{1.0};

  for (Intrepid2::Index i = 0; i < jac_check.size(); ++i) {
    error1 = std::min(error1, jac_check[i][3]);
    error2 = std::min(error2, jac_check[i][3]);
  }

  Real const
  tol{1.0e-07};

  ASSERT_LE(error1, tol);
  ASSERT_LE(error2, tol);
}
