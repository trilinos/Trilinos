// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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

#include "gtest/gtest.h"
#include "Intrepid2_MiniTensor_FunctionSet.h"

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

namespace Intrepid2 {

namespace {

//
// Test the solution methods by themselves.
//

// Test one function with one method.
template <typename STEP, typename FN, typename T, Index N>
bool
solveFNwithSTEP(STEP & step_method, FN & function, Vector<T, N> & x)
{
  Minimizer<T, N>
  minimizer;

  minimizer.solve(step_method, function, x);

  minimizer.printReport(std::cout);

  return minimizer.converged;
}

// Test one system with various methods.
template <typename FN, typename T, Index N>
bool
solveFN(FN & function, Vector<T, N> const & x)
{
  bool
  all_ok = true;

  Vector<T, N>
  y;

  NewtonStep<FN, T, N>
  newton_step;

  y = x;

  bool const
  newton_ok = solveFNwithSTEP(newton_step, function, y);

  all_ok = all_ok && newton_ok;

  TrustRegionStep<FN, T, N>
  trust_region_step;

  y = x;

  bool const
  trust_region_ok = solveFNwithSTEP(trust_region_step, function, y);

  all_ok = all_ok && trust_region_ok;

  ConjugateGradientStep<FN, T, N>
  pcg_step;

  y = x;

  bool const
  pcg_ok = solveFNwithSTEP(pcg_step, function, y);

  all_ok = all_ok && pcg_ok;

  LineSearchRegularizedStep<FN, T, N>
  line_search_step;

  y = x;

  bool const
  line_search_ok = solveFNwithSTEP(line_search_step, function, y);

  all_ok = all_ok && line_search_ok;
  
  NewtonWithLineSearchStep<FN, T, N>
  newton_line_search_step;

  y = x;

  bool const
  newton_line_search_ok = solveFNwithSTEP(newton_line_search_step, function, y);

  all_ok = all_ok && newton_line_search_ok;
  
  return all_ok;
}

// Test various systems with various methods.
bool testSystemsAndMethods()
{
  constexpr Index
  max_dimension{2};

  bool
  all_ok = true;

  Vector<Real, max_dimension>
  x;

  SquareRoot<Real>
  square_root(2.0);

  x.set_dimension(SquareRoot<Real>::DIMENSION);

  x(0) = 10.0;

  bool const
  square_root_ok = solveFN(square_root, x);

  all_ok = all_ok && square_root_ok;

  Quadratic<Real>
  quadratic(10.0, 15.0, 1.0);

  x.set_dimension(Quadratic<Real>::DIMENSION);

  x(0) = -15.0;
  x(1) = -10.0;

  bool const
  quadratic_ok = solveFN(quadratic, x);

  all_ok = all_ok && quadratic_ok;

  Gaussian<Real>
  gaussian(1.0, 2.0, 0.125);

  x.set_dimension(Gaussian<Real>::DIMENSION);

  x(0) = 0.0;
  x(1) = 0.0;

  bool const
  gaussian_ok = solveFN(gaussian, x);

  all_ok = all_ok && gaussian_ok;

  Banana<Real>
  banana;

  x.set_dimension(Banana<Real>::DIMENSION);

  x(0) = 0.0;
  x(1) = 3.0;

  bool const
  banana_ok = solveFN(banana, x);

  all_ok = all_ok && banana_ok;

  Matyas<Real>
  matyas;

  x.set_dimension(Matyas<Real>::DIMENSION);

  x(0) = 10.0;
  x(1) =  0.0;

  bool const
  matyas_ok = solveFN(matyas, x);

  all_ok = all_ok && matyas_ok;

  McCormick<Real>
  mccormick;

  x.set_dimension(McCormick<Real>::DIMENSION);

  x(0) = -0.5;
  x(1) = -1.5;

  bool const
  mccormick_ok = solveFN(mccormick, x);

  all_ok = all_ok && mccormick_ok;

  StyblinskiTang<Real>
  styblinski_tang;

  x.set_dimension(StyblinskiTang<Real>::DIMENSION);

  x(0) = -4.0;
  x(1) = -4.0;

  bool const
  styblinski_tang_ok = solveFN(styblinski_tang, x);

  all_ok = all_ok && styblinski_tang_ok;

  Paraboloid<Real>
  paraboloid;

  x.set_dimension(Paraboloid<Real>::DIMENSION);

  x(0) = 128.0;
  x(1) = 256.0;;

  bool const
  paraboloid_ok = solveFN(paraboloid, x);

  all_ok = all_ok && paraboloid_ok;

  Beale<Real>
  beale;

  x.set_dimension(Beale<Real>::DIMENSION);

  x(0) = -4.5;
  x(1) = -4.5;

  bool const
  beale_ok = solveFN(beale, x);

  all_ok = all_ok && beale_ok;

  Booth<Real>
  booth;

  x.set_dimension(Booth<Real>::DIMENSION);

  x(0) = -10.0;
  x(1) = -10.0;

  bool const
  booth_ok = solveFN(booth, x);

  all_ok = all_ok && booth_ok;

  GoldsteinPrice<Real>
  goldstein_price;

  x.set_dimension(GoldsteinPrice<Real>::DIMENSION);

  x(0) = 2.0;
  x(1) = 2.0;

  bool const
  goldstein_price_ok = solveFN(goldstein_price, x);

  all_ok = all_ok && goldstein_price_ok;

  return all_ok;
}

} // anonymous namespace

TEST(NonlinearSystems, NonlinearMethods)
{
  bool const
  passed = testSystemsAndMethods();

  ASSERT_EQ(passed, true);
}

TEST(Testing, OptimizationMethods)
{
  constexpr Index
  dimension{2};

  using MIN = Minimizer<Real, dimension>;
  using FN = Banana<Real>;
  using STEP = NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  FN
  banana;

  STEP
  step;

  Vector<Real, dimension>
  x;

  x(0) = 0.0;
  x(1) = 3.0;

  minimizer.solve(step, banana, x);

  minimizer.printReport(std::cout);

  ASSERT_EQ(minimizer.converged, true);
}

TEST(Testing, ValueGradientHessian)
{
  constexpr Index
  dimension{2};

  Paraboloid<Real>
  p;

  Vector<Real, dimension> const
  x(0.0, 0.0);

  Real const
  f = p.value(x);

  Vector<Real, dimension> const
  df = p.gradient(x);

  Tensor<Real, dimension> const
  ddf = p.hessian(x);

  std::cout << "Point   : " << x << '\n';
  std::cout << "Value   : " << f << '\n';
  std::cout << "Gradient: " << df << '\n';
  std::cout << "Hessian : " << ddf << '\n';

  Tensor<Real, dimension> const
  I = identity<Real, dimension>(dimension);

  Real const
  error = std::sqrt(f) + norm(df) + norm(ddf - 2.0 * I);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(Testing, MixedStorage)
{
  Index const
  dimension{2};

  std::cout << '\n';

  Vector<Real, 3>
  v(1.0, 2.0, 3.0);

  v.set_dimension(dimension);

  std::cout << "Vector   : " << v << '\n';

  Tensor<Real, 3>
  A(1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 7.0, 8.0, 9.0);

  A.set_dimension(dimension);

  std::cout << "Tensor   : " << A << '\n';

  Matrix<Real, 3, 4>
  B(ONES);

  B.set_dimensions(4, 2);

  std::cout << "Matrix   : " << B << '\n';

  bool const
  passed = v.get_dimension() == dimension && A.get_dimension() == dimension &&
    B.get_num_rows() == 4 && B.get_num_cols() == 2;

  ASSERT_EQ(passed, true);
}

TEST(Testing, FailedFlag)
{
  constexpr Index
  dimension{1};

  using MIN = Minimizer<Real, dimension>;
  using FN = Failure<Real>;
  using STEP = NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  FN
  fn;

  STEP
  step;

  Vector<Real, dimension>
  x;

  x(0) = 0.0;

  minimizer.solve(step, fn, x);

  ASSERT_EQ(minimizer.failed, true);
}

TEST(Testing, Monotonicity)
{
  constexpr Index
  dimension{1};

  using MIN = Minimizer<Real, dimension>;
  using FN = Mesa<Real>;
  using STEP = NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  minimizer.enforce_monotonicity = true;

  FN
  fn;

  STEP
  step;

  Vector<Real, dimension>
  x;

  x(0) = 2.0;

  minimizer.solve(step, fn, x);

  ASSERT_EQ(minimizer.monotonic, false);
  ASSERT_EQ(minimizer.failed, true);
}

TEST(Testing, Boundedness)
{
  constexpr Index
  dimension{1};

  using MIN = Minimizer<Real, dimension>;
  using FN = Sigmoid<Real>;
  using STEP = NewtonStep<FN, Real, dimension>;

  MIN
  minimizer;

  minimizer.enforce_boundedness = true;

  FN
  fn;

  STEP
  step;

  Vector<Real, dimension>
  x;

  x(0) = 0.5;

  minimizer.solve(step, fn, x);

  minimizer.printReport(std::cout);

  ASSERT_EQ(minimizer.bounded, true);
  ASSERT_EQ(minimizer.failed, false);
}

TEST(Testing, ConstraintIdentity)
{
  constexpr Index
  num_rows{2};

  constexpr Index
  num_cols{2};

  Identity<Real, num_rows>
  id;

  Vector<Real, num_cols> const
  x(ZEROS);

  Vector<Real, num_rows> const
  f = id.value(x);

  Matrix<Real, num_rows, num_cols> const
  df = id.gradient(x);

  std::cout << "Point   : " << x << '\n';
  std::cout << "Value   : " << f << '\n';
  std::cout << "Gradient: " << df << '\n';

  Real
  error{0.0};

  for (Index i = 0; i < min(num_rows, num_cols); ++i) {
    error += (df(i, i) - 1.0) * (df(i, i) - 1.0);
  }

  error = std::sqrt(error);

  ASSERT_LE(error, machine_epsilon<Real>());
}

} // namespace Intrepid2
