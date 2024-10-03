// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "gtest/gtest.h"
#include "MiniTensor_FunctionSet.h"
#include "Teuchos_oblackholestream.hpp"

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

namespace minitensor {

namespace {

//
// Test the solution methods by themselves.
//

// Test one function with one method.
template <typename STEP, typename FN, typename T, Index N>
bool
solveFNwithSTEP(STEP & step_method, FN & function, Vector<T, N> & x)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  Minimizer<T, N>
  minimizer;

  minimizer.solve(step_method, function, x);

  minimizer.printReport(os);

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

TEST(MiniTensor, LinearSolver)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Index
  DIM{11};

  Tensor<Real, DIM> const
  A = 2.0 * eye<Real, DIM>() + Tensor<Real, DIM>(Filler::RANDOM_UNIFORM);

  os << "\n\nMatrix A:" << A;

  Vector<Real, DIM> const
  x(Filler::RANDOM_UNIFORM);

  os << "\n\nVector x:" << x;

  Vector<Real, DIM> const
  b = A * x;

  os << "\n\nVector b = A * x:" << b;

  Vector<Real, DIM> const
  y = solve(A, b);

  os << "\n\nVector y = solve(A, b):" << y;

  Real const
  error = norm(y - x);

  os << "\n\nerror = norm(y - x):" << error << '\n';

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (DIM - 1) * machine_epsilon<Real>();

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, Preconditioners)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  constexpr Index
  DIM{2};

  Tensor<Real, DIM> const
  A(2.0e32, 1.0e32, 1.0, 2.0);

  os << "\n\nMatrix A:" << A;

  Vector<Real, DIM> const
  b(1.0e32, 2.0);

  os << "\n\nVector b:" << b;

  Vector<Real, DIM> const
  x(0.0, 1.0);

  os << "\n\nVector x:" << x;

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (DIM - 1) * machine_epsilon<Real>();

  Vector<Real, DIM>
  y = solve(A, b, PreconditionerType::IDENTITY);

  os << "\n\nVector y = solve(A, b, IDENTITY):" << y;

  Real
  error = norm(y - x);

  os << "\n\nerror = norm(y - x):" << error << '\n';

  ASSERT_LE(error, tolerance);

  y = solve(A, b, PreconditionerType::DIAGONAL);

  os << "\n\nVector y = solve(A, b, DIAGONAL):" << y;

  error = norm(y - x);

  os << "\n\nerror = norm(y - x):" << error << '\n';

  ASSERT_LE(error, tolerance);

  y = solve(A, b, PreconditionerType::MAX_ABS_ROW);

  os << "\n\nVector y = solve(A, b, MAX_ABS_ROW):" << y;

  error = norm(y - x);

  os << "\n\nerror = norm(y - x):" << error << '\n';

  ASSERT_LE(error, tolerance);

  y = solve_full_pivot(A, b);

  os << "\n\nVector y = solve_full_pivot(A, b):" << y;

  error = norm(y - x);

  os << "\n\nerror = norm(y - x):" << error << '\n';

  ASSERT_LE(error, tolerance);
}

TEST(MiniTensor, NonlinearMethods)
{
  bool const
  passed = testSystemsAndMethods();

  ASSERT_EQ(passed, true);
}

TEST(MiniTensor, OptimizationMethods)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

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

  minimizer.printReport(os);

  ASSERT_EQ(minimizer.converged, true);
}

TEST(MiniTensor, ValueGradientHessian)
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

  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  os << "Point   : " << x << '\n';
  os << "Value   : " << f << '\n';
  os << "Gradient: " << df << '\n';
  os << "Hessian : " << ddf << '\n';

  Tensor<Real, dimension> const
  I = identity<Real, dimension>(dimension);

  Real const
  error = std::sqrt(f) + norm(df) + norm(ddf - 2.0 * I);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, MixedStorage)
{
  Index const
  dimension{2};

  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  os << '\n';

  Vector<Real, 3>
  v(1.0, 2.0, 3.0);

  v.set_dimension(dimension);

  os << "Vector   : " << v << '\n';

  Tensor<Real, 3>
  A(1.0, 2.0, 3.0, 4.0, 5.0, 5.0, 7.0, 8.0, 9.0);

  A.set_dimension(dimension);

  os << "Tensor   : " << A << '\n';

  Matrix<Real, 3, 4>
  B(Filler::ONES);

  B.set_dimensions(4, 2);

  os << "Matrix   : " << B << '\n';

  bool const
  passed = v.get_dimension() == dimension && A.get_dimension() == dimension &&
    B.get_num_rows() == 4 && B.get_num_cols() == 2;

  ASSERT_EQ(passed, true);
}

TEST(MiniTensor, FailedFlag)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

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

  os << minimizer.failure_message << '\n';

  ASSERT_EQ(minimizer.failed, true);
}

TEST(MiniTensor, Monotonicity)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

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

  os << minimizer.failure_message << '\n';

  ASSERT_EQ(minimizer.monotonic, false);
  ASSERT_EQ(minimizer.failed, true);
}

TEST(MiniTensor, Boundedness)
{
  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

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

  os << minimizer.failure_message << '\n';

  minimizer.printReport(os);

  ASSERT_EQ(minimizer.bounded, true);
  ASSERT_EQ(minimizer.failed, false);
}

TEST(MiniTensor, ConstraintIdentity)
{
  constexpr Index
  NUM_CONSTR{2};

  constexpr Index
  NUM_VAR{2};

  Identity<Real, NUM_CONSTR, NUM_VAR>
  id;

  Vector<Real, NUM_VAR> const
  x(Filler::ZEROS);

  Vector<Real, NUM_CONSTR> const
  f = id.value(x);

  Matrix<Real, NUM_CONSTR, NUM_VAR> const
  df = id.gradient(x);

  bool const
  print_output = ::testing::GTEST_FLAG(print_time);

  // outputs nothing
  Teuchos::oblackholestream
  bhs;

  std::ostream &
  os = (print_output == true) ? std::cout : bhs;

  os << "Point   : " << x << '\n';
  os << "Value   : " << f << '\n';
  os << "Gradient: " << df << '\n';

  Real
  error{0.0};

  for (Index i = 0; i < min(NUM_CONSTR, NUM_VAR); ++i) {
    error += (df(i, i) - 1.0) * (df(i, i) - 1.0);
  }

  error = std::sqrt(error);

  ASSERT_LE(error, machine_epsilon<Real>());
}

} // namespace minitensor
