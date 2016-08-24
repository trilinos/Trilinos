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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#include <gtest/gtest.h>
#include <Intrepid2_MiniTensor_ROL_Vector.h>

typedef double T;

int
main(int ac, char * av[])
{
  Kokkos::initialize();

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

TEST(MiniTensor_ROL, VectorAdaptor)
{
  int
  error_flag{0};

  T const
  epsilon{Intrepid2::machine_epsilon<T>()};

  T const
  error_tol{16.0 * epsilon};

  constexpr Intrepid2::Index
  N{16};

  Intrepid2::Vector<T, N>
  vx(Intrepid2::RANDOM);

  Intrepid2::Vector<T, N>
  vy(Intrepid2::RANDOM);

  Intrepid2::Vector<T, N>
  vz(Intrepid2::RANDOM);

  ROL::MiniTensorVector<T, N>
  x(vx);

  ROL::MiniTensorVector<T, N>
  y(vy);

  ROL::MiniTensorVector<T, N>
  z(vz);

  // Standard tests.
  std::vector<T>
  consistency = x.checkVector(y, z, true, std::cout);

  Intrepid2::Index const
  num_tests = consistency.size();

  Intrepid2::Vector<T>
  checkvec(num_tests);

  checkvec.fill(&consistency[0]);

  ASSERT_LE(Intrepid2::norm(checkvec), std::sqrt(epsilon));

  // Basis tests.
  // set x to first basis vector
  Teuchos::RCP<ROL::Vector<T>>
  w = x.clone();

  w = x.basis(0);

  T
  wnorm = w->norm();

  std::cout << "Norm of ROL::Vector w (first basis vector): ";
  std::cout << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // set x to middle basis vector
  w = x.basis(N / 2);

  wnorm = w->norm();

  std::cout << "\nNorm of ROL::Vector w ('middle' basis vector): ";
  std::cout << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // set x to last basis vector
  w = x.basis(N - 1);

  wnorm = w->norm();

  std::cout << "\nNorm of ROL::Vector w (last basis vector): ";
  std::cout << wnorm << "\n";

  ASSERT_LE(std::abs(wnorm - 1.0), error_tol);

  // Repeat the checkVector tests with a zero vector.
  x.scale(0.0);

  consistency = x.checkVector(x, x, true, std::cout);

  checkvec.fill(&consistency[0]);

  ASSERT_EQ(Intrepid2::norm(checkvec), 0.0);
}

