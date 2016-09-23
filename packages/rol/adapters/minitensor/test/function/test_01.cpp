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

TEST(MiniTensor_ROL, ROL_Gradient)
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

  ROL::MiniTensorVector<Real, DIM>
  x(xval);

  ROL::MiniTensorVector<Real, DIM>
  d(dval);

  std::vector<std::vector<Real>>
  grad_check = obj.checkGradient(x, d);

  Real
  error{1.0};

  for (Intrepid2::Index i = 0; i < grad_check.size(); ++i) {
    if (i == 0) {
      os << "\n";
      os << std::right;
      os << std::setw(20) << "Step size";
      os << std::setw(20) << "grad'*dir";
      os << std::setw(20) << "FD approx";
      os << std::setw(20) << "abs error";
      os << "\n";
    }
    os << std::scientific << std::setprecision(8) << std::right;
    os << std::setw(20) << grad_check[i][0];
    os << std::setw(20) << grad_check[i][1];
    os << std::setw(20) << grad_check[i][2];
    os << std::setw(20) << grad_check[i][3];
    os << "\n";
    error = std::min(error, grad_check[i][3] * grad_check[i][3]);
  }

  error = std::sqrt(error);

  Real const
  tol{1.0e-6};

  ASSERT_LE(error, tol);
}
