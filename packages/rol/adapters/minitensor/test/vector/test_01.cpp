// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
