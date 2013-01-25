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

#include<ctime>

#include "Intrepid_FieldContainer.hpp"
#include "Sacado.hpp"
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

namespace Intrepid
{

  TEUCHOS_UNIT_TEST(MiniTensor, Initialization)
  {
    FieldContainer<Real> FC(3, 3);
    FC(0, 0) = 1.0;
    FC(0, 1) = 2.0;
    FC(0, 2) = 3.0;
    FC(1, 0) = 4.0;
    FC(1, 1) = 5.0;
    FC(1, 2) = 6.0;
    FC(2, 0) = 7.0;
    FC(2, 1) = 8.0;
    FC(2, 2) = 9.0;

    Real const * dataPtr0 = &FC(0, 0);

    Index const N = 3;
    Vector<Real> u(N, dataPtr0);

    TEST_COMPARE( u(0), ==, 1.0);
    TEST_COMPARE( u(1), ==, 2.0);
    TEST_COMPARE( u(2), ==, 3.0);

    Real const * dataPtr1 = &FC(1, 0);

    u = Vector<Real>(N, dataPtr1);

    TEST_COMPARE( u(0), ==, 4.0);
    TEST_COMPARE( u(1), ==, 5.0);
    TEST_COMPARE( u(2), ==, 6.0);

    Real const * dataPtr2 = &FC(2, 0);

    u = Vector<Real>(N, dataPtr2);

    TEST_COMPARE( u(0), ==, 7.0);
    TEST_COMPARE( u(1), ==, 8.0);
    TEST_COMPARE( u(2), ==, 9.0);
  }

  TEUCHOS_UNIT_TEST(MiniTensor, VectorAddition)
  {
    Vector<Real> const u(1.0, 0.0, 0.0);
    Vector<Real> const v(0.0, 1.0, 0.0);
    Vector<Real> const w(1.0, 1.0, 0.0);

    TEST_COMPARE( u + v == w, !=, 0);
  }

  TEUCHOS_UNIT_TEST(MiniTensor, VectorSubtraction)
  {
    Vector<Real> u(3);
    Vector<Real> v(3);
    u(0) = 1.0;
    u(1) = 2.0;
    u(2) = 3.0;

    v = u - u;

    TEST_COMPARE(norm(v), <=, machine_epsilon<Real>());
  }

  TEUCHOS_UNIT_TEST(MiniTensor, VectorScalarMultipliaction)
  {
    Vector<Real> u(3);
    Vector<Real> v(3);
    Vector<Real> w(3);
    u(0) = 1.0;
    u(1) = 2.0;
    u(2) = 3.0;

    v(0) = -2.0;
    v(1) = -4.0;
    v(2) = -6.0;

    w = 4.0 * u + 2.0 * v;

    TEST_COMPARE( norm(w), <=, machine_epsilon<Real>());
  }

  TEUCHOS_UNIT_TEST(MiniTensor, TensorInstantiation)
  {
    FieldContainer<Real> FC(2, 3, 3);
    FC(0, 0, 0) = 1.0;
    FC(0, 0, 1) = 2.0;
    FC(0, 0, 2) = 3.0;
    FC(0, 1, 0) = 4.0;
    FC(0, 1, 1) = 5.0;
    FC(0, 1, 2) = 6.0;
    FC(0, 2, 0) = 7.0;
    FC(0, 2, 1) = 8.0;
    FC(0, 2, 2) = 9.0;
    FC(1, 0, 0) = 10.0;
    FC(1, 0, 1) = 11.0;
    FC(1, 0, 2) = 12.0;
    FC(1, 1, 0) = 13.0;
    FC(1, 1, 1) = 14.0;
    FC(1, 1, 2) = 15.0;
    FC(1, 2, 0) = 16.0;
    FC(1, 2, 1) = 17.0;
    FC(1, 2, 2) = 18.0;

    Real const * dataPtr0 = &FC(0, 0, 0);

    Tensor<Real> const A(3, dataPtr0);

    TEST_COMPARE( A(0,0), ==, 1.0);
    TEST_COMPARE( A(0,1), ==, 2.0);
    TEST_COMPARE( A(0,2), ==, 3.0);
    TEST_COMPARE( A(1,0), ==, 4.0);
    TEST_COMPARE( A(1,1), ==, 5.0);
    TEST_COMPARE( A(1,2), ==, 6.0);
    TEST_COMPARE( A(2,0), ==, 7.0);
    TEST_COMPARE( A(2,1), ==, 8.0);
    TEST_COMPARE( A(2,2), ==, 9.0);

    Real const * dataPtr1 = &FC(1, 0, 0);

    Tensor<Real> const B(3, dataPtr1);

    TEST_COMPARE( B(0,0), ==, 10.0);
    TEST_COMPARE( B(0,1), ==, 11.0);
    TEST_COMPARE( B(0,2), ==, 12.0);
    TEST_COMPARE( B(1,0), ==, 13.0);
    TEST_COMPARE( B(1,1), ==, 14.0);
    TEST_COMPARE( B(1,2), ==, 15.0);
    TEST_COMPARE( B(2,0), ==, 16.0);
    TEST_COMPARE( B(2,1), ==, 17.0);
    TEST_COMPARE( B(2,2), ==, 18.0);
  }

  TEUCHOS_UNIT_TEST(MiniTensor, TensorAddition)
  {
    Tensor<Real> const A(3, 1.0);
    Tensor<Real> const B(3, 2.0);
    Tensor<Real> const C(3, 3.0);

    TEST_COMPARE( C == A + B, !=, 0);
  }

  TEUCHOS_UNIT_TEST(MiniTensor, Inverse)
  {
    std::srand(std::time(NULL));
    Index const N = double(std::rand()) / double(RAND_MAX) * 7.0 + 3.0;
    Tensor<Real> A(N);
    Tensor<Real> B(N);
    Tensor<Real> C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        A(i, j) = double(std::rand()) / double(RAND_MAX) * 20.0 - 10.0;
      }
    }

    B = inverse(A);

    C = A * B;

    Real const error = norm(C - eye<Real>(N)) / norm(A);

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

    Tensor<Real> Rref(3, 0.0);
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

    Tensor<Real> R1(3, 0.0);
    Tensor<Real> logR2(3, 0.0);
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

  TEUCHOS_UNIT_TEST(MiniTensor, Inverse4x4)
  {
    Tensor<Real> A = 2.0 * identity<Real>(4);

    A(0, 1) = 1.0;
    A(1, 0) = 1.0;

    A(1, 2) = 1.0;
    A(2, 1) = 1.0;

    A(2, 3) = 1.0;
    A(3, 2) = 1.0;

    Tensor<Real> const B = inverse(A);

    Tensor<Real> const C = A * B;

    Tensor<Real> const I = eye<Real>(4);

    Real const error = norm(C - I) / norm(A);

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

  TEUCHOS_UNIT_TEST(MiniTensor, SVD3x3Fad)
  {
    Tensor<Sacado::Fad::DFad<double> > A(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
        9.0);

    Tensor<Sacado::Fad::DFad<double> > U(3), S(3), V(3);

    boost::tie(U, S, V) = svd(A);

    Tensor<Sacado::Fad::DFad<double> > B = U * S * transpose(V);

    Sacado::Fad::DFad<double> const error = norm(B - A) / norm(A);

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

    Real error = abs(P(1, 0) - 100.0) / 100.0;

    TEST_COMPARE(error, <=, machine_epsilon<Real>());

    sigma = piola_inverse(F, P);

    error = abs(sigma(1, 1) - 50.0) / 50.0;

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

} // namespace Intrepid
