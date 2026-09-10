// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Tests for the container classes (Vector, Tensor, Tensor3, Tensor4,
// Matrix): construction, filling, arithmetic, mixed scalar types, and
// the type-traits machinery.

#include <ctime>
#include <vector>

#include "gtest/gtest.h"
#include "MiniTensor.h"

int
main(int ac, char* av[])
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

template<typename T>
std::vector<T>
generate_sequence(
    Index const number_elements, T const & start, T const & increment)
{
  std::vector<T>
  v(number_elements);

  for (Index i = 0; i < number_elements; ++i) {
    v[i] = start + i * increment;
  }

  return v;
}

template<typename Tensor, typename Scalar>
bool
test_fundamentals(Index const dimension)
{
  bool
  passed = true;

  Index const
  number_components = integer_power(dimension, Tensor::ORDER);

  std::vector<Scalar> const
  X = generate_sequence<Scalar>(number_components, 1.0, 1.0);

  // Test constructor with pointer
  Tensor const
  A(dimension, &X[0]);

  // Test copy constructor
  Tensor
  B = A;

  Tensor
  C;

  // Test copy assignment
  C = B - A;

  Scalar
  error = norm_f(C);

  bool const
  copy_assigned = error <= machine_epsilon<Scalar>();
  passed = passed && copy_assigned;

  // Test fill with pointer
  B.fill(&X[0]);

  C = B - A;

  error = norm_f(C);

  bool const
  filled_pointer = error <= machine_epsilon<Scalar>();
  passed = passed && filled_pointer;

  std::vector<Scalar> const
  Y = generate_sequence<Scalar>(number_components, -1.0, -1.0);

  C.fill(&Y[0]);

  // Test increment
  C += A;

  error = norm_f(C);

  bool const
  incremented = error <= machine_epsilon<Scalar>();
  passed = passed && incremented;

  C.fill(&X[0]);

  // Test decrement
  C -= A;

  error = norm_f(C);

  bool const
  decremented = error <= machine_epsilon<Scalar>();
  passed = passed && decremented;

  //test Tensor fill and create for Kokkos data types
  Kokkos::View<Scalar *, Kokkos::DefaultHostExecutionSpace>
  X1("X1_kokkos", dimension);

  Kokkos::View<Scalar **, Kokkos::DefaultHostExecutionSpace>
  X2("X2_kokkos", dimension, dimension);

  Kokkos::View<Scalar ***, Kokkos::DefaultHostExecutionSpace>
  X3("X3_kokkos", dimension, dimension, dimension);

  Kokkos::View<Scalar ****, Kokkos::DefaultHostExecutionSpace>
  X4("X4_kokkos", dimension, dimension, dimension, dimension);

  Kokkos::deep_copy(X1, 3.1);
  Kokkos::deep_copy(X2, 3.2);
  Kokkos::deep_copy(X3, 3.3);
  Kokkos::deep_copy(X4, 3.4);

  Tensor
  Z(dimension);

  Index
  rank = 0;

  Index
  temp = number_components;

  while (temp != 1) {
    temp = temp / dimension;
    rank = rank + 1;
    assert(temp > 0);
  }

  switch (rank) {
  default:
    assert(false);
    break;

  case 1:
    Z.fill(X1, 0);
    break;

  case 2:
    Z.fill(X2, 0, 0);
    break;

  case 3:
    Z.fill(X3, 0, 0, 0);
    break;

  case 4:
    Z.fill(X4, 0, 0, 0, 0);
    break;
  }

  // Test copy constructor.
  Tensor const
  U = Z;

  // Test copy assignment.
  Tensor
  V;

  V = U - Z;

  error = norm_f(V);

  bool const
  passed_copy = error <= machine_epsilon<Scalar>();
  passed = passed && passed_copy;

  return passed;
}

template<typename Tensor, typename Scalar>
bool
test_filling(Index const dimension)
{
  bool
  passed = true;

  Index const
  number_components = integer_power(dimension, Tensor::ORDER);

  // Test construct with zeros
  Tensor
  A(dimension, Filler::ZEROS);

  Real
  error = norm_f_square(A);

  bool const
  zeros_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_constructed;

  // Test construct with ones
  Tensor
  B(dimension, Filler::ONES);

  error = norm_f_square(B) - number_components;

  bool const
  ones_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && ones_constructed;

  // Test construct with random entries
  Tensor
  C(dimension, Filler::RANDOM_UNIFORM);

  error = norm_f(C);

  bool const
  random_constructed = error > 0.0 && error < number_components;
  passed = passed && random_constructed;

  // Test fill with random components
  A.fill(Filler::RANDOM_UNIFORM);

  error = norm_f(A);

  bool const
  random_filled = error > 0.0 && error < number_components;
  passed = passed && random_filled;

  // Test fill with zeros
  B.fill(Filler::ZEROS);

  error = norm_f_square(B);

  bool const
  zeros_filled = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_filled;

  // Test fill with ones
  C.fill(Filler::ZEROS);

  error = norm_f_square(C) - number_components;

  bool const
  ones_filled = error <= machine_epsilon<Scalar>();
  passed = passed && ones_filled;

  return passed;
}

template<typename Tensor, typename Scalar>
bool
test_arithmetic(Index const dimension)
{
  bool
  passed = true;

  Index const
  number_components = integer_power(dimension, Tensor::ORDER);

  std::vector<Scalar> const
  X = generate_sequence<Scalar>(number_components, 1.0, 1.0);

  Real const
  sum_squares = number_components * (number_components + 1) *
      (2 * number_components + 1) / 6;

  // Test addition
  Tensor const
  A(dimension, &X[0]);

  Tensor const
  B = -1.0 * A;

  Tensor const
  C = -1.0 * B;

  Tensor const
  D = A + B;

  Real
  error = norm_f_square(D);

  bool const
  added = error <= machine_epsilon<Scalar>();
  passed = passed && added;

  // Test subtraction
  Tensor const
  E = A - C;

  error = norm_f_square(E);

  bool const
  subtracted = error <= machine_epsilon<Scalar>();
  passed = passed && subtracted;

  // Test scaling
  error = norm_f_square(C) - sum_squares;

  bool const
  scaled = error <= machine_epsilon<Scalar>();
  passed = passed && scaled;

  Tensor const
  F = C / -1.0;

  error = norm_f_square(F) - sum_squares;

  bool const
  divided = error <= machine_epsilon<Scalar>();
  passed = passed && divided;

  Tensor const
  G = 1.0 / C;

  error = norm_f_square(G) - sum_squares;

  bool const
  split = error <= machine_epsilon<Scalar>();
  passed = passed && split;

  return passed;
}

template<typename Matrix, typename Scalar>
bool
test_fundamentals(Index const rows, Index const cols)
{
  bool
  passed = true;

  Index const
  number_components = rows * cols;

  std::vector<Scalar> const
  X = generate_sequence<Scalar>(number_components, 1.0, 1.0);

  // Test constructor with pointer
  Matrix const
  A(rows, cols, &X[0]);

  // Test copy constructor
  Matrix
  B = A;

  Matrix
  C;

  // Test copy assignment
  C = B - A;

  Scalar
  error = norm_f(C);

  bool const
  copy_assigned = error <= machine_epsilon<Scalar>();
  passed = passed && copy_assigned;

  // Test fill with pointer
  B.fill(&X[0]);

  C = B - A;

  error = norm_f(C);

  bool const
  filled_pointer = error <= machine_epsilon<Scalar>();
  passed = passed && filled_pointer;

  std::vector<Scalar> const
  Y = generate_sequence<Scalar>(number_components, -1.0, -1.0);

  C.fill(&Y[0]);

  // Test increment
  C += A;

  error = norm_f(C);

  bool const
  incremented = error <= machine_epsilon<Scalar>();
  passed = passed && incremented;

  C.fill(&X[0]);

  // Test decrement
  C -= A;

  error = norm_f(C);

  bool const
  decremented = error <= machine_epsilon<Scalar>();
  passed = passed && decremented;

  return passed;
}

template<typename Matrix, typename Scalar>
bool
test_filling(Index const rows, Index const cols)
{
  bool
  passed = true;

  Index const
  number_components = rows * cols;

  // Test construct with zeros
  Matrix
  A(rows, cols, Filler::ZEROS);

  Real
  error = norm_f_square(A);

  bool const
  zeros_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_constructed;

  // Test construct with ones
  Matrix
  B(rows, cols, Filler::ONES);

  error = norm_f_square(B) - number_components;

  bool const
  ones_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && ones_constructed;

  // Test construct with random entries
  Matrix
  C(rows, cols, Filler::RANDOM_UNIFORM);

  error = norm_f(C);

  bool const
  random_constructed = error > 0.0 && error < number_components;
  passed = passed && random_constructed;

  // Test fill with random components
  A.fill(Filler::RANDOM_UNIFORM);

  error = norm_f(A);

  bool const
  random_filled = error > 0.0 && error < number_components;
  passed = passed && random_filled;

  // Test fill with zeros
  B.fill(Filler::ZEROS);

  error = norm_f_square(B);

  bool const
  zeros_filled = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_filled;

  // Test fill with ones
  C.fill(Filler::ZEROS);

  error = norm_f_square(C) - number_components;

  bool const
  ones_filled = error <= machine_epsilon<Scalar>();
  passed = passed && ones_filled;

  return passed;
}

template<typename Matrix, typename Scalar>
bool
test_arithmetic(Index const rows, Index const cols)
{
  bool
  passed = true;

  Index const
  number_components = rows * cols;

  std::vector<Scalar> const
  X = generate_sequence<Scalar>(number_components, 1.0, 1.0);

  Real const
  sum_squares = number_components * (number_components + 1) *
      (2 * number_components + 1) / 6;

  // Test addition
  Matrix const
  A(rows, cols, &X[0]);

  Matrix const
  B = -1.0 * A;

  Matrix const
  C = -1.0 * B;

  Matrix const
  D = A + B;

  Real
  error = norm_f_square(D);

  bool const
  added = error <= machine_epsilon<Scalar>();
  passed = passed && added;

  // Test subtraction
  Matrix const
  E = A - C;

  error = norm_f_square(E);

  bool const
  subtracted = error <= machine_epsilon<Scalar>();
  passed = passed && subtracted;

  // Test scaling
  error = norm_f_square(C) - sum_squares;

  bool const
  scaled = error <= machine_epsilon<Scalar>();
  passed = passed && scaled;

  Matrix const
  F = C / -1.0;

  error = norm_f_square(F) - sum_squares;

  bool const
  divided = error <= machine_epsilon<Scalar>();
  passed = passed && divided;

  Matrix const
  G = 1.0 / C;

  error = norm_f_square(G) - sum_squares;

  bool const
  split = error <= machine_epsilon<Scalar>();
  passed = passed && split;

  return passed;
}

} // anonymous namespace

TEST(MiniTensor, Fundamentals)
{
  bool const
  vector_dynamic_passed = test_fundamentals<Vector<Real>, Real>(3);

  ASSERT_EQ(vector_dynamic_passed, true);

  bool const
  vector_static_passed = test_fundamentals<Vector<Real, 3>, Real>(3);

  ASSERT_EQ(vector_static_passed, true);

  bool const
  tensor_dynamic_passed = test_fundamentals<Tensor<Real>, Real>(3);

  ASSERT_EQ(tensor_dynamic_passed, true);

  bool const
  tensor_static_passed = test_fundamentals<Tensor<Real, 3>, Real>(3);

  ASSERT_EQ(tensor_static_passed, true);

  bool const
  tensor3_dynamic_passed = test_fundamentals<Tensor3<Real>, Real>(3);

  ASSERT_EQ(tensor3_dynamic_passed, true);

  bool const
  tensor3_static_passed = test_fundamentals<Tensor3<Real, 3>, Real>(3);

  ASSERT_EQ(tensor3_static_passed, true);

  bool const
  tensor4_dynamic_passed = test_fundamentals<Tensor4<Real>, Real>(3);

  ASSERT_EQ(tensor4_dynamic_passed, true);

  bool const
  tensor4_static_passed = test_fundamentals<Tensor4<Real, 3>, Real>(3);

  ASSERT_EQ(tensor4_static_passed, true);

  bool const
  matrix_dynamic_passed = test_fundamentals<Matrix<Real>, Real>(4, 3);

  ASSERT_EQ(matrix_dynamic_passed, true);

  bool const
  matrix_static_passed = test_fundamentals<Matrix<Real, 4, 3>, Real>(4, 3);

  ASSERT_EQ(matrix_static_passed, true);
}

TEST(MiniTensor, Filling)
{
  bool const
  vector_dynamic_passed = test_filling<Vector<Real>, Real>(3);

  ASSERT_EQ(vector_dynamic_passed, true);

  bool const
  vector_static_passed = test_filling<Vector<Real, 3>, Real>(3);

  ASSERT_EQ(vector_static_passed, true);

  bool const
  tensor_dynamic_passed = test_filling<Tensor<Real>, Real>(3);

  ASSERT_EQ(tensor_dynamic_passed, true);

  bool const
  tensor_static_passed = test_filling<Tensor<Real, 3>, Real>(3);

  ASSERT_EQ(tensor_static_passed, true);

  bool const
  tensor3_dynamic_passed = test_filling<Tensor3<Real>, Real>(3);

  ASSERT_EQ(tensor3_dynamic_passed, true);

  bool const
  tensor3_static_passed = test_filling<Tensor3<Real, 3>, Real>(3);

  ASSERT_EQ(tensor3_static_passed, true);

  bool const
  tensor4_dynamic_passed = test_filling<Tensor4<Real>, Real>(3);

  ASSERT_EQ(tensor4_dynamic_passed, true);

  bool const
  tensor4_static_passed = test_filling<Tensor4<Real, 3>, Real>(3);

  ASSERT_EQ(tensor4_static_passed, true);

  bool const
  matrix_dynamic_passed = test_filling<Matrix<Real>, Real>(4, 3);

  ASSERT_EQ(matrix_dynamic_passed, true);

  bool const
  matrix_static_passed = test_filling<Matrix<Real, 4, 3>, Real>(4, 3);

  ASSERT_EQ(matrix_static_passed, true);
}

TEST(MiniTensor, Arithmetic)
{
  bool const
  vector_dynamic_passed = test_arithmetic<Vector<Real>, Real>(3);

  ASSERT_EQ(vector_dynamic_passed, true);

  bool const
  vector_static_passed = test_arithmetic<Vector<Real, 3>, Real>(3);

  ASSERT_EQ(vector_static_passed, true);

  bool const
  tensor_dynamic_passed = test_arithmetic<Tensor<Real>, Real>(3);

  ASSERT_EQ(tensor_dynamic_passed, true);

  bool const
  tensor_static_passed = test_arithmetic<Tensor<Real, 3>, Real>(3);

  ASSERT_EQ(tensor_static_passed, true);

  bool const
  tensor3_dynamic_passed = test_arithmetic<Tensor3<Real>, Real>(3);

  ASSERT_EQ(tensor3_dynamic_passed, true);

  bool const
  tensor3_static_passed = test_arithmetic<Tensor3<Real, 3>, Real>(3);

  ASSERT_EQ(tensor3_static_passed, true);

  bool const
  tensor4_dynamic_passed = test_arithmetic<Tensor4<Real>, Real>(3);

  ASSERT_EQ(tensor4_dynamic_passed, true);

  bool const
  tensor4_static_passed = test_arithmetic<Tensor4<Real, 3>, Real>(3);

  ASSERT_EQ(tensor4_static_passed, true);

  bool const
  matrix_dynamic_passed = test_arithmetic<Matrix<Real>, Real>(4, 3);

  ASSERT_EQ(matrix_dynamic_passed, true);

  bool const
  matrix_static_passed = test_arithmetic<Matrix<Real, 4, 3>, Real>(4, 3);

  ASSERT_EQ(matrix_static_passed, true);
}

TEST(MiniTensor, MixedTypes)
{
  Tensor<Sacado::Fad::DFad<Real>>
  A(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  Tensor<Sacado::Fad::DFad<Real>> const
  B(3, Filler::ONES);

  Tensor<Real> const
  C(3, Filler::ONES);

  Real const
  b = 1.0;

  Sacado::Fad::DFad<Real> const
  c = 1.0;

  A += b * B;

  A -= c * C;

  Sacado::Fad::DFad<Real>
  error = norm_f_square(A) - 3.0;

  ASSERT_LE(error, machine_epsilon<Real>());

  A = B + C;

  error = norm_f(A) - 6.0;

  ASSERT_LE(error, machine_epsilon<Real>());

  A = C - B;

  error = norm_f(A);

  ASSERT_LE(error, machine_epsilon<Real>());

  A += C;

  error = norm_f(A) - 3.0;

  ASSERT_LE(error, machine_epsilon<Real>());

  A -= C;

  error = norm_f(A);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, KroneckerProduct)
{
  Tensor4<Real> const A = identity_3<Real>(3);

  Tensor<Real> const Q = eye<Real>(3);

  Tensor4<Real> const B = kronecker(Q, A);

  Real const error = norm_f(B - A) / norm_f(A);

  ASSERT_LE(error, machine_epsilon<Real>());
}

TEST(MiniTensor, TemplateMetaProgramming)
{
  {
    Real
    a = 0.0;

    Sacado::Fad::DFad<Real>
    b = 0.0;

    Real
    c = Sacado::ScalarValue<Real>::eval(a);

    //std::cout << c << '\n';

    Real
    d = Sacado::ScalarValue<Sacado::Fad::DFad<Real>>::eval(b);

    //std::cout << d << '\n';

    bool const
    is_equal = c == d;

    ASSERT_EQ(is_equal, true);
  }

  {
    Vector<Real>
    A(3, Filler::ZEROS);

    Vector<Sacado::Fad::DFad<Real>>
    B(3, Filler::ZEROS);

    Vector<Real>
    C = Sacado::ScalarValue<Vector<Real>>::eval(A);

    //std::cout << C << '\n';

    Vector<Real>
    D = Sacado::ScalarValue<Vector<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    ASSERT_EQ(is_equal, true);
  }

  {
    Tensor<Real>
    A(3, Filler::ZEROS);

    Tensor<Sacado::Fad::DFad<Real>>
    B(3, Filler::ZEROS);

    Tensor<Real>
    C = Sacado::ScalarValue<Tensor<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor<Real>
    D = Sacado::ScalarValue<Tensor<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    ASSERT_EQ(is_equal, true);
  }

  {
    Tensor3<Real>
    A(3, Filler::ZEROS);

    Tensor3<Sacado::Fad::DFad<Real>>
    B(3, Filler::ZEROS);

    Tensor3<Real>
    C = Sacado::ScalarValue<Tensor3<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor3<Real>
    D = Sacado::ScalarValue<Tensor3<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    ASSERT_EQ(is_equal, true);
  }

  {
    Tensor4<Real>
    A(3, Filler::ZEROS);

    Tensor4<Sacado::Fad::DFad<Real>>
    B(3, Filler::ZEROS);

    Tensor4<Real>
    C = Sacado::ScalarValue<Tensor4<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor4<Real>
    D = Sacado::ScalarValue<Tensor4<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    ASSERT_EQ(is_equal, true);
  }

#if !defined(KOKKOS_ENABLE_CUDA)
  {
    //
    // use double explicitly
    //
    using A = Vector<double>;

    using B = Vector<Sacado::Fad::DFad<double>>;

    using C = Vector<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>;

    std::string const
    double_string = "double";

    std::string const
    fad_string = "Sacado::Fad::Exp::GeneralFad< double >";

    std::string
    type_string =
        Sacado::StringName<Sacado::ScalarType<A>::type>::eval();

    ASSERT_EQ(type_string, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<A>::type>::eval();

    ASSERT_EQ(type_string, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<B>::type>::eval();

    ASSERT_EQ(type_string, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<B>::type>::eval();

    ASSERT_EQ(type_string, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<C>::type>::eval();

    ASSERT_EQ(type_string, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<C>::type>::eval();

    ASSERT_EQ(type_string, fad_string);
  }
#endif // KOKKOS_ENABLE_CUDA

}

} // namespace minitensor
