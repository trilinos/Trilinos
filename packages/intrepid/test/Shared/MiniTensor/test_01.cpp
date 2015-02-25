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
#include <vector>

#include "Intrepid_FieldContainer.hpp"
#ifdef HAVE_INTREPID_KOKKOSCORE
#include "Sacado.hpp"
#else
#include "Sacado_No_Kokkos.hpp"
#endif
#include "Intrepid_MiniTensor.h"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  std::cout << "End Result: TEST PASSED";
  std::cout << std::endl;
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
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

template <typename Tensor, typename Scalar>
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

 
#ifdef HAVE_INTREPID_KOKKOSCORE
  //test Tensor fill and create for Kokkos data types
  Kokkos::View<Scalar*, Kokkos::Serial> X1_k ("X1_kokkos",dimension);
  Kokkos::View<Scalar**, Kokkos::Serial> X2_k ("X2_kokkos",dimension,dimension);
  Kokkos::View<Scalar***, Kokkos::Serial> X3_k ("X3_kokkos",dimension,dimension,dimension);
  Kokkos::View<Scalar****, Kokkos::Serial> X4_k ("X4_kokkos",dimension,dimension, dimension, dimension);

  Kokkos::deep_copy(X1_k, 3.1);
  Kokkos::deep_copy(X2_k, 3.2);
  Kokkos::deep_copy(X3_k, 3.3);
  Kokkos::deep_copy(X4_k, 3.4);

  Tensor  A1_k(dimension); //(X1_k,0);

 // Index const  number_components = A1_k.get_number_components();
  Index rank=0;
  Index temp=number_components;

  while (temp!=1){
   temp =temp/dimension;
   rank=rank +1;
   if (temp<1) TEUCHOS_TEST_FOR_EXCEPTION( ( (temp<1)  ), std::invalid_argument,
                                  ">>> ERROR (MiniTensor:: fill): rank calculation is not correct");
  }

  if (rank==1)
   A1_k.fill(X1_k,0);

  if (rank==2)
   A1_k.fill(X2_k,0,0);
 
  if (rank==3)
   A1_k.fill(X3_k,0,0,0);

  if (rank==4)
   A1_k.fill(X4_k,0,0,0,0);

  Tensor B1_k=A1_k;

  Tensor C1_k;
  C1_k = B1_k - A1_k;

  error = norm_f(C1_k);

  bool const
  tensor_create_from_1d_kokkos = error <= machine_epsilon<Scalar>();
  passed = passed && tensor_create_from_1d_kokkos;
#endif 


  return passed;
}

template <typename Tensor, typename Scalar>
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

template <typename Tensor, typename Scalar>
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
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse2x2)
{
  Tensor<Real, 2> const
  A = 2.0 * eye<Real, 2>() + Tensor<Real, 2>(RANDOM_UNIFORM);

  Tensor<Real, 2> const
  B = inverse(A);

  Tensor<Real, 2> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, 2>()) / norm(A);

  TEST_COMPARE(error, <=, 100.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse3x3)
{
  Tensor<Real, 3> const
  A = 2.0 * eye<Real, 3>() + Tensor<Real, 3>(RANDOM_UNIFORM);

  Tensor<Real, 3> const
  B = inverse(A);

  Tensor<Real, 3> const
  C = A * B;

  Real const
  error = norm(C - eye<Real, 3>()) / norm(A);

  TEST_COMPARE(error, <=, 100.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, InverseNxN)
{
  Index const
  N = 11;

  Tensor<Real> const
  A = 2.0 * eye<Real>(N) + Tensor<Real>(N, RANDOM_UNIFORM);

  Tensor<Real> const
  B = inverse(A);

  Tensor<Real> const
  C = A * B;

  Real const
  error = norm(C - eye<Real>(N)) / norm(A);

  TEST_COMPARE(error, <=, 100.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Inverse_4th_NxN)
{
  Index const
  N = 4;

  Tensor4<Real> const
  A = 2.0 * identity_1<Real>(N) + Tensor4<Real>(N, RANDOM_UNIFORM);

  Tensor4<Real> const
  B = inverse(A);

  Tensor4<Real> const
  C = dotdot(A, B);

  Real const
  error = norm_f(C - identity_1<Real>(N)) / norm_f(A);

  TEST_COMPARE(error, <=, 100.0 * machine_epsilon<Real>());
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
  Tensor<Real> const A(1, 2, 3, 4, 5, 6, 7, 8, 9);

  Tensor<Real> const B = exp_pade(A);

  Tensor<Real> const C = exp_taylor(A);

  Tensor<Real> const D = B - C;

  Real const error = norm(D) / norm(B);

  TEST_COMPARE( error, <=, 100.0 * machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SymmetricEigen)
{
  Tensor<Real> A = eye<Real>(3);
  A(0, 1) = 0.1;
  A(1, 0) = 0.1;

  Tensor<Real> V(3);
  Tensor<Real> D(3);

  boost::tie(V, D) = eig_sym(A);

  TEST_COMPARE(std::abs(D(0,0) - 1.1), <=, machine_epsilon<Real>());
  TEST_COMPARE(std::abs(D(1,1) - 1.0), <=, machine_epsilon<Real>());
  TEST_COMPARE(std::abs(D(2,2) - 0.9), <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, LeftPolarDecomposition)
{
  Tensor<Real> V0(1.1, 0.2, 0.0, 0.2, 1.0, 0.0, 0.0, 0.0, 1.2);

  Tensor<Real> R0(sqrt(2) / 2, -sqrt(2) / 2, 0.0, sqrt(2) / 2, sqrt(2) / 2,
      0.0, 0.0, 0.0, 1.0);

  Tensor<Real> F = V0 * R0;
  Tensor<Real> V(3);
  Tensor<Real> R(3);
  boost::tie(V, R) = polar_left(F);

  TEST_COMPARE(norm(V-V0), <=, 10.0*machine_epsilon<Real>());
  TEST_COMPARE(norm(R-R0), <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, LogRotation)
{
  Tensor<Real> R = identity<Real>(3);
  Tensor<Real> R0(sqrt(2) / 2, -sqrt(2) / 2, 0.0, sqrt(2) / 2, sqrt(2) / 2,
      0.0, 0.0, 0.0, 1.0);

  Tensor<Real> r = log_rotation(R);
  Tensor<Real> r0 = log_rotation(R0);

  TEST_COMPARE(norm(r), <=, machine_epsilon<Real>());

  TEST_COMPARE( std::abs(r0(0,1) + 0.785398163397448), <=,
      10.0*machine_epsilon<Real>());

  TEST_COMPARE( std::abs(r0(0,1) + r0(1,0)), <=, machine_epsilon<Real>());

  Real theta = std::acos(-1.0) + 10 * machine_epsilon<Real>();

  R(0, 0) = cos(theta);
  R(1, 1) = cos(theta);
  R(0, 1) = sin(theta);
  R(1, 0) = -sin(theta);
  R(2, 2) = 1.0;

  Tensor<Real> logR = log_rotation(R);

  Tensor<Real> Rref(3, ZEROS);
  Rref(0, 1) = -theta;
  Rref(1, 0) = theta;

  TEST_COMPARE(norm(logR - Rref), <=, 100*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, BakerCampbellHausdorff)
{
  Tensor<Real> F = 3.0 * identity<Real>(3);
  Tensor<Real> V(3), R(3), logV(3), logR(3);

  boost::tie(V, R, logV) = polar_left_logV(F);
  logR = log_rotation(R);

  Tensor<Real> f = bch(logV, logR);

  TEST_COMPARE( std::abs(f(0,0) - std::log(3.0)), <=,
      machine_epsilon<Real>());

  Vector<Real> u(3);
  u(0) = std::acos(-1.0) / std::sqrt(2.0);
  u(1) = u(0);
  u(2) = 0.0;

  Tensor<Real> R1(3, ZEROS);
  Tensor<Real> logR2(3, ZEROS);
  logR2(0, 2) = u(1);
  logR2(1, 2) = -u(0);
  logR2(2, 0) = -u(1);
  logR2(2, 1) = u(0);
  logR2(0, 1) = -u(2);
  logR2(1, 0) = u(2);

  R1 = exp_skew_symmetric(logR2);
  Tensor<Real> Rref = zero<Real>(3);
  Rref(0, 1) = 1.0;
  Rref(1, 0) = 1.0;
  Rref(2, 2) = -1.0;

  TEST_COMPARE( norm(Rref-R1), <=, 100.0*machine_epsilon<Real>());
  TEST_COMPARE( norm(exp_skew_symmetric(logR) - R), <=,
      100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, PolarLeftLog)
{
  Tensor<Real> const F(3.60070151614402, 0.00545892068653966,
      0.144580850331452, -5.73345529510674, 0.176660910549112,
      1.39627497290058, 2.51510445213514, 0.453212159218359,
      -1.44616077859513);

  Tensor<Real> const L(0.265620603957487, -1.066921781600734,
      -0.089540974250415, -1.066921781600734, 0.927394431410918,
      -0.942214085118614, -0.089540974250415, -0.942214085118613,
      0.105672693695746);

  Tensor<Real> V(3), R(3), v(3), r(3);

  boost::tie(V, R, v) = polar_left_logV(F);

  Real const error = norm(v - L) / norm(L);

  TEST_COMPARE( error, <=, 100*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, VolumetricDeviatoric)
{
  Tensor<Real> A = 3.0 * eye<Real>(3);

  TEST_COMPARE( norm(A - vol(A)), <=, 100.0*machine_epsilon<Real>());

  Tensor<Real> B = dev(A);

  A(0, 0) = 0.0;
  A(1, 1) = 0.0;
  A(2, 2) = 0.0;

  TEST_COMPARE( norm(A - B), <=, 100.0*machine_epsilon<Real>());
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

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SVD3x3)
{
  Tensor<Real> const A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);

  Tensor<Real> U(3), S(3), V(3);

  boost::tie(U, S, V) = svd(A);

  Tensor<Real> const B = U * S * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SVD3x3Fad)
{
  Tensor<Sacado::Fad::DFad<Real> > A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
      9.0);

  Tensor<Sacado::Fad::DFad<Real> > U(3), S(3), V(3);

  boost::tie(U, S, V) = svd(A);

  Tensor<Sacado::Fad::DFad<Real> > B = U * S * transpose(V);

  Sacado::Fad::DFad<Real> const error = norm(B - A) / norm(A);

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, MixedTypes)
{
  Tensor<Sacado::Fad::DFad<Real> >
  A(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  Tensor<Sacado::Fad::DFad<Real> > B(3, ONES);

  Tensor<Real> C(3, ONES);

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

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, SymmetricEigen3x3)
{
  Tensor<Real> const A(2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0);

  Tensor<Real> V(3), D(3);

  boost::tie(V, D) = eig_sym(A);

  Tensor<Real> const B = V * D * transpose(V);

  Real const error = norm(A - B) / norm(A);

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Polar3x3)
{
  Tensor<Real> A(2.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0);

  Tensor<Real> R(3), U(3);

  boost::tie(R, U) = polar_right(A);

  Tensor<Real> X(3), D(3), Y(3);

  boost::tie(X, D, Y) = svd(A);

  Tensor<Real> B = R - X * transpose(Y) + U - Y * D * transpose(Y);

  Real const error = norm(B) / norm(A);

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, Cholesky)
{
  Tensor<Real> A(1.0, 1.0, 1.0, 1.0, 5.0, 3.0, 1.0, 3.0, 3.0);

  Tensor<Real> G(3);

  bool is_spd;

  boost::tie(G, is_spd) = cholesky(A);

  Tensor<Real> B(1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 1.0, 1.0);

  Real const error = norm(G - B) / norm(A);

  TEST_COMPARE(error, <=, 100.0*machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, MechanicsTransforms)
{
  Tensor<Real> F(0.0, -6.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0 / 3.0);

  Tensor<Real> sigma(0.0, 0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0);

  Tensor<Real> P = piola(F, sigma);

  Real error = std::abs(P(1, 0) - 100.0) / 100.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  sigma = piola_inverse(F, P);

  error = std::abs(sigma(1, 1) - 50.0) / 50.0;

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  Tensor<Real> E = 0.5 * (t_dot(F, F) - eye<Real>(3));

  Tensor<Real> e = 0.5 * (eye<Real>(3) - inverse(dot_t(F, F)));

  Tensor<Real> g = push_forward_covariant(F, E);

  error = norm(g - e) / norm(e);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());

  Tensor<Real> G = pull_back_covariant(F, e);

  error = norm(G - E) / norm(E);

  TEST_COMPARE(error, <=, machine_epsilon<Real>());
}

TEUCHOS_UNIT_TEST(MiniTensor, KroneckerProduct)
{
  Tensor4<Real> A = identity_3<Real>(3);

  Tensor<Real> Q = eye<Real>(3);

  Tensor4<Real> B = kronecker(Q, A);

  Real const error = norm_f(B-A) / norm_f(A);

  TEST_COMPARE(error, <=, 100.0 * machine_epsilon<Real>());
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
    d = Sacado::ScalarValue<Sacado::Fad::DFad<Real> >::eval(b);

    //std::cout << d << '\n';

    bool const
    is_equal = c == d;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Vector<Real>
    A(3, ZEROS);

    Vector<Sacado::Fad::DFad<Real> >
    B(3, ZEROS);

    Vector<Real>
    C = Sacado::ScalarValue<Vector<Real> >::eval(A);

    //std::cout << C << '\n';

    Vector<Real>
    D = Sacado::ScalarValue<Vector<Sacado::Fad::DFad<Real> > >::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor<Real>
    A(3, ZEROS);

    Tensor<Sacado::Fad::DFad<Real> >
    B(3, ZEROS);

    Tensor<Real>
    C = Sacado::ScalarValue<Tensor<Real> >::eval(A);

    //std::cout << C << '\n';

    Tensor<Real>
    D = Sacado::ScalarValue<Tensor<Sacado::Fad::DFad<Real> > >::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor3<Real>
    A(3, ZEROS);

    Tensor3<Sacado::Fad::DFad<Real> >
    B(3, ZEROS);

    Tensor3<Real>
    C = Sacado::ScalarValue<Tensor3<Real> >::eval(A);

    //std::cout << C << '\n';

    Tensor3<Real>
    D = Sacado::ScalarValue<Tensor3<Sacado::Fad::DFad<Real> > >::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }

  {
    Tensor4<Real>
    A(3, ZEROS);

    Tensor4<Sacado::Fad::DFad<Real> >
    B(3, ZEROS);

    Tensor4<Real>
    C = Sacado::ScalarValue<Tensor4<Real> >::eval(A);

    //std::cout << C << '\n';

    Tensor4<Real>
    D = Sacado::ScalarValue<Tensor4<Sacado::Fad::DFad<Real> > >::eval(B);

    //std::cout << D << '\n';

    bool const
    is_equal = C == D;

    TEST_COMPARE(is_equal, ==, true);
  }


  {
    //
    // use double explicitly
    //
    typedef Vector<double> A;

    typedef Vector<Sacado::Fad::DFad<double> > B;

    typedef Vector<Sacado::Fad::DFad<Sacado::Fad::DFad<double> > > C;

    std::string const
    double_string = "double";

    std::string const
    fad_string = "Sacado::Fad::DFad< double >";

    std::string
    type_string =
        Sacado::StringName<Sacado::ScalarType<A>::type >::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<A>::type >::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<B>::type >::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<B>::type >::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ScalarType<C>::type >::eval();

    TEST_COMPARE(type_string, ==, double_string);

    type_string =
        Sacado::StringName<Sacado::ValueType<C>::type >::eval();

    TEST_COMPARE(type_string, ==, fad_string);
  }
}

} // namespace Intrepid
