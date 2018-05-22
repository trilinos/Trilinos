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
#include <ROL_MiniTensor_Vector.hpp>

using T = double;

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

TEST(MiniTensor_ROL, VectorAdaptor)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  T const
  epsilon{minitensor::machine_epsilon<T>()};

  T const
  error_tol{16.0 * epsilon};

  constexpr minitensor::Index
  N{16};

  minitensor::Vector<T, N>
  vx(minitensor::Filler::RANDOM);

  minitensor::Vector<T, N>
  vy(minitensor::Filler::RANDOM);

  minitensor::Vector<T, N>
  vz(minitensor::Filler::RANDOM);

  ROL::MiniTensorVector<T, N>
  x(vx);

  ROL::MiniTensorVector<T, N>
  y(vy);

  ROL::MiniTensorVector<T, N>
  z(vz);

  // Standard tests.
  std::vector<T>
  consistency = x.checkVector(y, z, true, os);

  minitensor::Index const
  num_tests = consistency.size();

  minitensor::Vector<T>
  checkvec(num_tests);

  checkvec.fill(&consistency[0]);

  ASSERT_LE(minitensor::norm(checkvec), std::sqrt(epsilon));

  // Basis tests.
  // set x to first basis vector
  ROL::Ptr<ROL::Vector<T>>
  w = x.clone();

  w = x.basis(0);

  T
  wnorm = w->norm();

  os << "Norm of ROL::Vector w (first basis vector): ";
  os << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // set x to middle basis vector
  w = x.basis(N / 2);

  wnorm = w->norm();

  os << "\nNorm of ROL::Vector w ('middle' basis vector): ";
  os << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // set x to last basis vector
  w = x.basis(N - 1);

  wnorm = w->norm();

  os << "\nNorm of ROL::Vector w (last basis vector): ";
  os << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // Repeat the checkVector tests with a zero vector.
  x.scale(0.0);

  consistency = x.checkVector(x, x, true, os);

  checkvec.fill(&consistency[0]);

  ASSERT_EQ(minitensor::norm(checkvec), 0.0);
}

TEST(MiniTensor_ROL, VectorValue)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  ROL::nullstream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  T const
  epsilon{minitensor::machine_epsilon<T>()};

  constexpr
  minitensor::Index
  N{4};

  minitensor::Vector<T, N> const
  x(minitensor::Filler::ONES);

  os << "minitensor::Vector x     : " << x << '\n';

  ROL::MiniTensorVector<T, N>
  y(x);

  os << "ROL::MiniTensorVector y : " << y << '\n';

  T const
  error0 = std::abs(y.norm() - minitensor::norm(x));

  ASSERT_LE(error0, epsilon);

  ROL::Vector<T> &
  z = y;

  minitensor::Vector<T, N> const
  a = ROL::MTfromROL<T, N>(z);

  os << "minitensor::Vector a     : " << a << '\n';

  T const
  error1 = minitensor::norm(x - a);

  ASSERT_LE(error1, epsilon);

  minitensor::Vector<T, N> const
  b(minitensor::Filler::SEQUENCE);

  os << "minitensor::Vector b     : " << b << '\n';

  ROL::MTtoROL(b, z);

  minitensor::Vector<T, N> const
  c = ROL::MTfromROL<T, N>(z);

  os << "minitensor::Vector c     : " << c << '\n';

  T const
  error2 = minitensor::norm(b - c);

  ASSERT_LE(error2, epsilon);
}
