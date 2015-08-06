// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

#include <ctime>
#include <sstream>
#include <vector>

#include "Sacado.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_MiniTensor.h"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Sacado.hpp"

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::RCP<std::stringstream>
  oss = Teuchos::rcp(new std::stringstream);

  Teuchos::FancyOStream
  fos(oss);

  const bool
  success = Teuchos::UnitTestRepository::runUnitTests(fos);

  if (success == true) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  else {
    std::cout << "\nEnd Result: TEST FAILED" << std::endl;
  }

  return success == true ? 0 : 1;
}

namespace Intrepid {

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
  A(dimension, ZEROS);

  Real
  error = norm_f_square(A);

  bool const
  zeros_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_constructed;

  // Test construct with ones
  Tensor
  B(dimension, ONES);

  error = norm_f_square(B) - number_components;

  bool const
  ones_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && ones_constructed;

  // Test construct with random entries
  Tensor
  C(dimension, RANDOM_UNIFORM);

  error = norm_f(C);

  bool const
  random_constructed = error > 0.0 && error < number_components;
  passed = passed && random_constructed;

  // Test fill with random components
  A.fill(RANDOM_UNIFORM);

  error = norm_f(A);

  bool const
  random_filled = error > 0.0 && error < number_components;
  passed = passed && random_filled;

  // Test fill with zeros
  B.fill(ZEROS);

  error = norm_f_square(B);

  bool const
  zeros_filled = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_filled;

  // Test fill with ones
  C.fill(ZEROS);

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
  A(rows, cols, ZEROS);

  Real
  error = norm_f_square(A);

  bool const
  zeros_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_constructed;

  // Test construct with ones
  Matrix
  B(rows, cols, ONES);

  error = norm_f_square(B) - number_components;

  bool const
  ones_constructed = error <= machine_epsilon<Scalar>();
  passed = passed && ones_constructed;

  // Test construct with random entries
  Matrix
  C(rows, cols, RANDOM_UNIFORM);

  error = norm_f(C);

  bool const
  random_constructed = error > 0.0 && error < number_components;
  passed = passed && random_constructed;

  // Test fill with random components
  A.fill(RANDOM_UNIFORM);

  error = norm_f(A);

  bool const
  random_filled = error > 0.0 && error < number_components;
  passed = passed && random_filled;

  // Test fill with zeros
  B.fill(ZEROS);

  error = norm_f_square(B);

  bool const
  zeros_filled = error <= machine_epsilon<Scalar>();
  passed = passed && zeros_filled;

  // Test fill with ones
  C.fill(ZEROS);

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

TEUCHOS_UNIT_TEST(MiniTensor, Fundamentals)
{
  bool const
  vector_dynamic_passed = test_fundamentals<Vector<Real>, Real>(3);

  TEST_COMPARE(vector_dynamic_passed, ==, true);

  bool const
  vector_static_passed = test_fundamentals<Vector<Real, 3>, Real>(3);

  TEST_COMPARE(vector_static_passed, ==, true);

  bool const
  tensor_dynamic_passed = test_fundamentals<Tensor<Real>, Real>(3);

  TEST_COMPARE(tensor_dynamic_passed, ==, true);

  bool const
  tensor_static_passed = test_fundamentals<Tensor<Real, 3>, Real>(3);

  TEST_COMPARE(tensor_static_passed, ==, true);

  bool const
  tensor3_dynamic_passed = test_fundamentals<Tensor3<Real>, Real>(3);

  TEST_COMPARE(tensor3_dynamic_passed, ==, true);

  bool const
  tensor3_static_passed = test_fundamentals<Tensor3<Real, 3>, Real>(3);

  TEST_COMPARE(tensor3_static_passed, ==, true);

  bool const
  tensor4_dynamic_passed = test_fundamentals<Tensor4<Real>, Real>(3);

  TEST_COMPARE(tensor4_dynamic_passed, ==, true);

  bool const
  tensor4_static_passed = test_fundamentals<Tensor4<Real, 3>, Real>(3);

  TEST_COMPARE(tensor4_static_passed, ==, true);

  bool const
  matrix_dynamic_passed = test_fundamentals<Matrix<Real>, Real>(4, 3);

  TEST_COMPARE(matrix_dynamic_passed, ==, true);

  bool const
  matrix_static_passed = test_fundamentals<Matrix<Real, 4, 3>, Real>(4, 3);

  TEST_COMPARE(matrix_static_passed, ==, true);
}

TEUCHOS_UNIT_TEST(MiniTensor, Filling)
{
  bool const
  vector_dynamic_passed = test_filling<Vector<Real>, Real>(3);

  TEST_COMPARE(vector_dynamic_passed, ==, true);

  bool const
  vector_static_passed = test_filling<Vector<Real, 3>, Real>(3);

  TEST_COMPARE(vector_static_passed, ==, true);

  bool const
  tensor_dynamic_passed = test_filling<Tensor<Real>, Real>(3);

  TEST_COMPARE(tensor_dynamic_passed, ==, true);

  bool const
  tensor_static_passed = test_filling<Tensor<Real, 3>, Real>(3);

  TEST_COMPARE(tensor_static_passed, ==, true);

  bool const
  tensor3_dynamic_passed = test_filling<Tensor3<Real>, Real>(3);

  TEST_COMPARE(tensor3_dynamic_passed, ==, true);

  bool const
  tensor3_static_passed = test_filling<Tensor3<Real, 3>, Real>(3);

  TEST_COMPARE(tensor3_static_passed, ==, true);

  bool const
  tensor4_dynamic_passed = test_filling<Tensor4<Real>, Real>(3);

  TEST_COMPARE(tensor4_dynamic_passed, ==, true);

  bool const
  tensor4_static_passed = test_filling<Tensor4<Real, 3>, Real>(3);

  TEST_COMPARE(tensor4_static_passed, ==, true);

  bool const
  matrix_dynamic_passed = test_filling<Matrix<Real>, Real>(4, 3);

  TEST_COMPARE(matrix_dynamic_passed, ==, true);

  bool const
  matrix_static_passed = test_filling<Matrix<Real, 4, 3>, Real>(4, 3);

  TEST_COMPARE(matrix_static_passed, ==, true);
}

TEUCHOS_UNIT_TEST(MiniTensor, Arithmetic)
{
  bool const
  vector_dynamic_passed = test_arithmetic<Vector<Real>, Real>(3);

  TEST_COMPARE(vector_dynamic_passed, ==, true);

  bool const
  vector_static_passed = test_arithmetic<Vector<Real, 3>, Real>(3);

  TEST_COMPARE(vector_static_passed, ==, true);

  bool const
  tensor_dynamic_passed = test_arithmetic<Tensor<Real>, Real>(3);

  TEST_COMPARE(tensor_dynamic_passed, ==, true);

  bool const
  tensor_static_passed = test_arithmetic<Tensor<Real, 3>, Real>(3);

  TEST_COMPARE(tensor_static_passed, ==, true);

  bool const
  tensor3_dynamic_passed = test_arithmetic<Tensor3<Real>, Real>(3);

  TEST_COMPARE(tensor3_dynamic_passed, ==, true);

  bool const
  tensor3_static_passed = test_arithmetic<Tensor3<Real, 3>, Real>(3);

  TEST_COMPARE(tensor3_static_passed, ==, true);

  bool const
  tensor4_dynamic_passed = test_arithmetic<Tensor4<Real>, Real>(3);

  TEST_COMPARE(tensor4_dynamic_passed, ==, true);

  bool const
  tensor4_static_passed = test_arithmetic<Tensor4<Real, 3>, Real>(3);

  TEST_COMPARE(tensor4_static_passed, ==, true);

  bool const
  matrix_dynamic_passed = test_arithmetic<Matrix<Real>, Real>(4, 3);

  TEST_COMPARE(matrix_dynamic_passed, ==, true);

  bool const
  matrix_static_passed = test_arithmetic<Matrix<Real, 4, 3>, Real>(4, 3);

  TEST_COMPARE(matrix_static_passed, ==, true);
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse2x2)
{
  Index const
  N = 2;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  TEST_COMPARE(error, <=, tolerance);
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse3x3)
{
  Index const
  N = 3;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  TEST_COMPARE(error, <=, tolerance);
}

TEUCHOS_UNIT_TEST(MiniTensor, InverseNxN)
{
  Index const
  N = 11;

  Tensor<Real, N> const
  A = 2.0 * eye<Real, N>() + Tensor<Real, N>(RANDOM_UNIFORM);

  Tensor<Real, N> const
  B = inverse(A);

  Tensor<Real, N> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, N>()) / norm(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (N - 1) * machine_epsilon<Real>();

  TEST_COMPARE(error, <=, tolerance);
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse_4th_NxN)
{
  Index const
  N = 4;

  Tensor4<Real, N> const
  A = 2.0 * identity_1<Real, N>() + Tensor4<Real, N>(RANDOM_UNIFORM);

  Tensor4<Real, N> const
  B = inverse(A);

  Tensor4<Real, N> const
  C = dotdot(A, B);

  Real const
  error = norm_f(C - identity_1<Real, N>()) / norm_f(A);

  // See Golub & Van Loan, Matrix Computations 4th Ed., pp 122-123
  Real const
  tolerance = 2 * (2 * N - 1) * machine_epsilon<Real>();

  TEST_COMPARE(error, <=, tolerance);
}

TEUCHOS_UNIT_TEST(MiniTensor, TensorManipulation)
{
  Tensor<Real> A = eye<Real>(3);
  Tensor<Real> B(3);
  Tensor<Real> C(3);
  Vector<Real> u(3);

  A = 2.0 * A;
  A(1, 0) = A(0, 1) = 1.0;
  A(2, 1) = A(1, 2) = 1.0;

  B = inverse(A);

  C = A * B;

  TEST_COMPARE(norm(C - eye<Real>(3)), <=, machine_epsilon<Real>());

  Real I1_A = I1(A);
  Real I2_A = I2(A);
  Real I3_A = I3(A);

  u(0) = I1_A - 6;
  u(1) = I2_A - 10;
  u(2) = I3_A - 4;

  Real const error = norm(u);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Exponential)
{
  Tensor<Real> const A = eye<Real>(3) + Tensor<Real>(3, ONES);

  Tensor<Real> const B = exp_pade(A);

  Tensor<Real> const C = exp_taylor(A);

  Tensor<Real> const D = B - C;

  Real const error = norm(D) / norm(B);

  TEST_COMPARE(error, <=, 2.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SymmetricEigen)
{
  Tensor<Real> A = eye<Real>(3);
  A(0, 1) = 0.1;
  A(1, 0) = 0.1;

  Tensor<Real> V(3);
  Tensor<Real> D(3);

  boost::tie(V, D) = eig_sym(A);

  TEST_COMPARE(std::abs(D(0, 0) - 1.1), <=, machine_epsilon<Real>());
  TEST_COMPARE(std::abs(D(1, 1) - 1.0), <=, machine_epsilon<Real>());
  TEST_COMPARE(std::abs(D(2, 2) - 0.9), <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, LeftPolarDecomposition)
{
  Tensor<Real> const X(1.1, 0.2, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0, 1.2);

  Real const
  c = sqrt(2.0) / 2.0;

  Tensor<Real> const Y(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const F = X * Y;
  Tensor<Real> V(3);
  Tensor<Real> R(3);

  boost::tie(V, R) = polar_left(F);

  Real const
  error_x = norm(V - X) / norm(X);

  Real const
  error_y = norm(R - Y) / norm(Y);

  TEST_COMPARE(error_x, <=, machine_epsilon<Real>());
  TEST_COMPARE(error_y, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, LogRotation)
{
  // Identity rotation
  Tensor<Real>
  I = identity<Real>(3);

  Tensor<Real>
  i = log_rotation(I);

  Real const
  error_I = norm(i) / norm(I);

  TEST_COMPARE(error_I, <=, machine_epsilon<Real>());

  // Pi / 4 rotation about Z.
  Real const
  c = sqrt(2.0) / 2.0;

  Tensor<Real> const R(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const r = log_rotation(R);

  Real const
  Pi = std::acos(-1.0);

  Real const
  error_R = std::abs(r(0,1) + Pi / 4.0);

  TEST_COMPARE(error_R, <=, machine_epsilon<Real>());

  Real const
  error_r = std::abs(r(0,1) + r(1,0));

  TEST_COMPARE(error_r, <=, machine_epsilon<Real>());

  // Pi rotation about Z
  I(0, 0) = -1.0;
  I(1, 1) = -1.0;

  i = log_rotation(I);

  Real const
  error_pi = 0.5 * (std::abs(i(0, 1) + Pi) + std::abs(i(1, 0) - Pi));

  TEST_COMPARE(error_pi, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, BakerCampbellHausdorff)
{
  Real const
  Pi = std::acos(-1.0);

  Real const
  gamma = 0.1;

  Tensor<Real> const u(0, gamma, 0, gamma, 0, 0, 0, 0, 0);

  Tensor<Real> const r(0, -Pi/4, 0, Pi/4, 0, 0, 0, 0, 0);

  Tensor<Real> const f_bch = bch(r, u);

  Tensor<Real> const U = exp(u);

  Tensor<Real> const R = exp(r);

  Tensor<Real> const F = R * U;

  Tensor<Real> const f = log(F);

  Real const
  error = norm(f_bch - f) / norm(F);

  // Our implementation of the Baker-Campbell-Hausdorff
  // formula uses only 4 terms, so we expect some error here.
  Real const
  tolerance = 1.0e-3;

  TEST_COMPARE(error, <=, tolerance);
}

TEUCHOS_UNIT_TEST(MiniTensor, PolarLeftLog)
{
  Real const
  gamma = 0.1;

  Tensor<Real> const x(0, gamma, 0, gamma, 0, 0, 0, 0, 0);

  Tensor<Real> const X = exp(x);

  Real const
  c = sqrt(2.0) / 2.0;

  Tensor<Real> const Y(c, -c, 0.0, c, c, 0.0, 0.0, 0.0, 1.0);

  Tensor<Real> const F = X * Y;

  Tensor<Real> V(3), R(3), v(3);

  boost::tie(V, R, v) = polar_left_logV(F);

  Real const error = norm(v - x) / norm(x);

  TEST_COMPARE(error, <=, 16 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, VolumetricDeviatoric)
{
  Tensor<Real> A = 3.0 * eye<Real>(3);

  TEST_COMPARE(norm(A - vol(A)), <=, machine_epsilon<Real>());

  Tensor<Real> const B = dev(A);

  A(0, 0) = 0.0;
  A(1, 1) = 0.0;
  A(2, 2) = 0.0;

  TEST_COMPARE(norm(A - B), <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SVD2x2)
{
  Real const phi = 1.0;

  Real const psi = 2.0;

  Real const s0 = sqrt(3.0);

  Real const s1 = sqrt(2.0);

  Real const cl = cos(phi);

  Real const sl = sin(phi);

  Real const cr = cos(psi);

  Real const sr = sin(psi);

  Tensor<Real> const X(cl, -sl, sl, cl);

  Tensor<Real> const Y(cr, -sr, sr, cr);

  Tensor<Real> const D(s0, 0.0, 0.0, s1);

  Tensor<Real> const A = X * D * transpose(Y);

  Tensor<Real> U(2), S(2), V(2);

  boost::tie(U, S, V) = svd(A);

  Tensor<Real> B = U * S * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, 2.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SVD3x3)
{
  Tensor<Real> const A(1, 2, 3, 4, 5, 6, 7, 8, 9);

  Tensor<Real> U(3), S(3), V(3);

  boost::tie(U, S, V) = svd(A);

  Tensor<Real> const B = U * S * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SVD3x3Fad)
{
  Tensor<Sacado::Fad::DFad<Real>> const A(1, 2, 3, 4, 5, 6, 7, 8, 9);

  Tensor<Sacado::Fad::DFad<Real>> U(3), S(3), V(3);

  boost::tie(U, S, V) = svd(A);

  Tensor<Sacado::Fad::DFad<Real>> const B = U * S * transpose(V);

  Sacado::Fad::DFad<Real> const error = norm(B - A) / norm(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, MixedTypes)
{
  Tensor<Sacado::Fad::DFad<Real>> A(1, 0, 0, 0, 1, 0, 0, 0, 1);

  Tensor<Sacado::Fad::DFad<Real>> const B(3, ONES);

  Tensor<Real> const C(3, ONES);

  Real const
  b = 1.0;

  Sacado::Fad::DFad<Real> const
  c = 1.0;

  A += b * B;

  A -= c * C;

  Sacado::Fad::DFad<Real>
  error = norm_f_square(A) - 3.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  A = B + C;

  error = norm_f(A) - 6.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  A = C - B;

  error = norm_f(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  A += C;

  error = norm_f(A) - 3.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  A -= C;

  error = norm_f(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SymmetricEigen2x2)
{
  Tensor<Real> const A(2.0, 1.0, 1.0, 2.0);

  Tensor<Real> V(2), D(2);

  boost::tie(V, D) = eig_sym(A);

  Tensor<Real> const B = V * D * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SymmetricEigen3x3)
{
  Tensor<Real> const A(2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0);

  Tensor<Real> V(3), D(3);

  boost::tie(V, D) = eig_sym(A);

  Tensor<Real> const B = V * D * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, 4.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Polar3x3)
{
  Tensor<Real> const A(2.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0);

  Tensor<Real> R(3), U(3);

  boost::tie(R, U) = polar_right(A);

  Tensor<Real> X(3), D(3), Y(3);

  boost::tie(X, D, Y) = svd(A);

  Tensor<Real> const B = R - X * transpose(Y) + U - Y * D * transpose(Y);

  Real const error = norm(B) / norm(A);

  TEST_COMPARE(error, <=, 8.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Cholesky)
{
  Tensor<Real> const A(1.0, 1.0, 1.0, 1.0, 5.0, 3.0, 1.0, 3.0, 3.0);

  Tensor<Real> G(3);

  bool is_spd;

  boost::tie(G, is_spd) = cholesky(A);

  Tensor<Real> const B(1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 1.0, 1.0);

  Real const error = norm(G - B) / norm(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, MechanicsTransforms)
{
  Tensor<Real> const F(0.0, -6.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0);

  Tensor<Real> sigma(0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0);

  Tensor<Real> const P = piola(F, sigma);

  Real error = std::abs(P(1, 0) - 100.0) / 100.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  sigma = piola_inverse(F, P);

  error = std::abs(sigma(1, 1) - 50.0) / 50.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  Tensor<Real> const E = 0.5 * (t_dot(F, F) - eye<Real>(3));

  Tensor<Real> const e = 0.5 * (eye<Real>(3) - inverse(dot_t(F, F)));

  Tensor<Real> const g = push_forward_covariant(F, E);

  error = norm(g - e) / norm(e);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  Tensor<Real> const G = pull_back_covariant(F, e);

  error = norm(G - E) / norm(E);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, KroneckerProduct)
{
  Tensor4<Real> const A = identity_3<Real>(3);

  Tensor<Real> const Q = eye<Real>(3);

  Tensor4<Real> const B = kronecker(Q, A);

  Real const error = norm_f(B - A) / norm_f(A);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, TemplateMetaProgramming)
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

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Vector<Real>
    A(3, ZEROS);

    Vector<Sacado::Fad::DFad<Real>>
    B(3, ZEROS);

    Vector<Real>
    C = Sacado::ScalarValue<Vector<Real>>::eval(A);

    //std::cout << C << '\n';

    Vector<Real>
    D = Sacado::ScalarValue<Vector<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor<Real>
    A(3, ZEROS);

    Tensor<Sacado::Fad::DFad<Real>>
    B(3, ZEROS);

    Tensor<Real>
    C = Sacado::ScalarValue<Tensor<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor<Real>
    D = Sacado::ScalarValue<Tensor<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor3<Real>
    A(3, ZEROS);

    Tensor3<Sacado::Fad::DFad<Real>>
    B(3, ZEROS);

    Tensor3<Real>
    C = Sacado::ScalarValue<Tensor3<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor3<Real>
    D = Sacado::ScalarValue<Tensor3<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor4<Real>
    A(3, ZEROS);

    Tensor4<Sacado::Fad::DFad<Real>>
    B(3, ZEROS);

    Tensor4<Real>
    C = Sacado::ScalarValue<Tensor4<Real>>::eval(A);

    //std::cout << C << '\n';

    Tensor4<Real>
    D = Sacado::ScalarValue<Tensor4<Sacado::Fad::DFad<Real>>>::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

#ifndef KOKKOS_HAVE_CUDA
  {
    //
    // use double explicitly
    //
    typedef Vector<double> A;

    typedef Vector<Sacado::Fad::DFad<double>> B;

    typedef Vector<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>> C;

    std::string const
    double_string = "double";

    std::string const
    fad_string = "Sacado::Fad::DFad< double >";

    std::string
    type_string =
        Sacado::StringName<Sacado::ScalarType<A>::type>::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<A>::type>::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<B>::type>::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<B>::type>::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<C>::type>::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<C>::type>::eval();

    TEST_COMPARE(type_string, ==, fad_string);
  }
#endif

}

} // namespace Intrepid
