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
#include <ROL_MiniTensor_Function.hpp>

using Real = double;

int
main(int ac, char * av[])
{
  Kokkos::initialize();

  // Disables elapsed time and output by default.
  ::testing::GTEST_FLAG(print_time) = false;

  ::testing::InitGoogleTest(&ac, av);

  auto const
  retval = RUN_ALL_TESTS();

  Kokkos::finalize();

  return retval;
}

//
// Paraboloid of revolution
//
template<typename S>
class Paraboloid : public Intrepid2::Function_Base<Paraboloid<S>, S>
{
public:

  Paraboloid() {}

  static constexpr
  Intrepid2::Index
  DIMENSION{2};

  static constexpr
  char const * const
  NAME{"Paraboloid"};

  // Explicit value.
  template<typename T, Intrepid2::Index N>
  T
  value(Intrepid2::Vector<T, N> const & x)
  {
    Intrepid2::Index const
    dimension = x.get_dimension();

    assert(dimension == DIMENSION);

    T const
    f = (x(0) * x(0) + x(1) * x(1));

    return f;
  }

  // Default AD gradient.
  template<typename T, Intrepid2::Index N>
  Intrepid2::Vector<T, N>
  gradient(Intrepid2::Vector<T, N> const & x)
  {
    return Intrepid2::Function_Base<Paraboloid<S>, S>::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Intrepid2::Index N>
  Intrepid2::Tensor<T, N>
  hessian(Intrepid2::Vector<T, N> const & x)
  {
    return Intrepid2::Function_Base<Paraboloid<S>, S>::hessian(*this, x);
  }

};

TEST(MiniTensor_ROL, FunctionAdaptor)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Intrepid2::Index
  dimension{2};

  Paraboloid<Real>
  p;

  Intrepid2::Vector<Real, dimension> const
  x(0.0, 0.0);

  Real const
  f = p.value(x);

  Intrepid2::Vector<Real, dimension> const
  df = p.gradient(x);

  Intrepid2::Tensor<Real, dimension> const
  ddf = p.hessian(x);

  os << "Point   : " << x << '\n';
  os << "Value   : " << f << '\n';
  os << "Gradient: " << df << '\n';
  os << "Hessian : " << ddf << '\n';

  Intrepid2::Tensor<Real, dimension> const
  I = Intrepid2::identity<Real, dimension>(dimension);

  Real const
  error = std::sqrt(f) + Intrepid2::norm(df) + Intrepid2::norm(ddf - 2.0 * I);

  ASSERT_LE(error, Intrepid2::machine_epsilon<Real>());
}

