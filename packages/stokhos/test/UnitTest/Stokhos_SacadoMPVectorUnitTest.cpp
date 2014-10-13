// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include <Kokkos_Core.hpp>

//
// Currently this doesn't test:
//   * the device
//   * threaded storage (needs the device)
//   * strided storage with non-trivial stride
//

// Common setup for unit tests
template <typename VectorType>
struct UnitTestSetup {

  typedef VectorType vec_type;

  double rtol, atol;
  double crtol, catol;
  int sz;
  vec_type x, y, cx;
  double a;

  UnitTestSetup() {
    rtol = 1e-4;
    atol = 1e-5;
    crtol = 1e-12;
    catol = 1e-12;
    a = 3.1;
    sz = 8;

    // Create vector
    x.reset(sz);
    y.reset(sz);
    cx.reset(1);
    cx = a;
    for (int i=0; i<sz; i++) {
      x.fastAccessCoeff(i) = 0.1*i;
      y.fastAccessCoeff(i) = 0.25*i;
    }
  }
};

#define UNARY_UNIT_TEST(VEC, OP, OPNAME)                                \
  TEUCHOS_UNIT_TEST( VEC, OPNAME) {                                     \
    UTS::vec_type u = OP(setup.x);                                      \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = OP(setup.x.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_const) {                             \
    UTS::vec_type u = OP(setup.cx);                                     \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
      v.fastAccessCoeff(i) = OP(setup.cx.fastAccessCoeff(0));           \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_resize) {                            \
    UTS::vec_type u;                                                    \
    u = OP(setup.x);                                                    \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = OP(setup.x.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define BINARY_UNIT_TEST(VEC, OP, OPNAME)                               \
  TEUCHOS_UNIT_TEST( VEC, OPNAME) {                                     \
    UTS::vec_type u = setup.x OP setup.y;                               \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const) {                        \
    UTS::vec_type u = setup.a OP setup.y;                               \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.a OP setup.y.fastAccessCoeff(i);     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const) {                       \
    UTS::vec_type u = setup.x OP setup.a ;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.a;                                                        \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_both_const) {                        \
    UTS::vec_type u = setup.cx OP setup.cx;                             \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) OP             \
        setup.cx.fastAccessCoeff(0);                                    \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const2) {                       \
    UTS::vec_type u = setup.cx OP setup.x;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) OP             \
        setup.x.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const2) {                      \
    UTS::vec_type u = setup.x OP setup.cx;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.cx.fastAccessCoeff(0);                                    \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_resize) {                            \
    UTS::vec_type u;                                                    \
    u = setup.x OP setup.y;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const_resize) {                 \
    UTS::vec_type u;                                                    \
    u = setup.a OP setup.y;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.a OP                                 \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const_resize) {                \
    UTS::vec_type u;                                                    \
    u = setup.x OP setup.a;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.a;                                                        \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define BINARYFUNC_UNIT_TEST(VEC, OP, SOP, OPNAME)                      \
  TEUCHOS_UNIT_TEST( VEC, OPNAME) {                                     \
    UTS::vec_type u = OP(setup.x,setup.y);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.y.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const) {                        \
    UTS::vec_type u = OP(setup.a,setup.y);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.a,                               \
                                setup.y.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const) {                       \
    UTS::vec_type u = OP(setup.x,setup.a);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.a);                               \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_both_const) {                        \
    UTS::vec_type u = OP(setup.cx,setup.cx);                            \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.cx.fastAccessCoeff(0),           \
                                 setup.cx.fastAccessCoeff(0));          \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const2) {                       \
    UTS::vec_type u = OP(setup.cx,setup.x);                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.cx.fastAccessCoeff(0),           \
                                setup.x.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const2) {                      \
    UTS::vec_type u = OP(setup.x,setup.cx);                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.cx.fastAccessCoeff(0));           \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_resize) {                            \
    UTS::vec_type u;                                                    \
    u = OP(setup.x,setup.y);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.y.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_left_const_resize) {                 \
    UTS::vec_type u;                                                    \
    u = OP(setup.a,setup.y);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.a,                               \
                                setup.y.fastAccessCoeff(i));            \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_right_const_resize) {                \
    UTS::vec_type u;                                                    \
    u = OP(setup.x,setup.a);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.a);                               \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define OPASSIGN_UNIT_TEST(VEC, OP, OPNAME)                             \
  TEUCHOS_UNIT_TEST( VEC, OPNAME) {                                     \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.x;                                                       \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.x.fastAccessCoeff(i);               \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_const) {                             \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.a;                                                       \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.a;                                  \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_const2) {                            \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.cx;                                                      \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.cx.fastAccessCoeff(0);              \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, OPNAME##_resize) {                            \
    UTS::vec_type u = setup.a;                                          \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    u OP setup.x;                                                       \
    for (int i=0; i<setup.sz; i++) {                                    \
      v.fastAccessCoeff(i) = setup.a;                                   \
      v.fastAccessCoeff(i) OP setup.x.fastAccessCoeff(i);               \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define SAXPY_UNIT_TEST(VEC)                                            \
  TEUCHOS_UNIT_TEST( VEC, saxpy) {                                      \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.x*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.x.fastAccessCoeff(i)*setup.y.fastAccessCoeff(i);          \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, saxpy_resize) {                               \
    UTS::vec_type u = setup.cx;                                         \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    u += setup.x*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) +              \
        setup.x.fastAccessCoeff(i)*setup.y.fastAccessCoeff(i);          \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, saxpy_const) {                                \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.a*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.a*setup.y.fastAccessCoeff(i);                             \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC, saxpy_const2) {                               \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.cx*setup.y;                                              \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.cx.fastAccessCoeff(0)*setup.y.fastAccessCoeff(i);         \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define VECTOR_UNIT_TESTS(VEC)                                  \
  UNARY_UNIT_TEST(VEC, +, UnaryPlus)                            \
  UNARY_UNIT_TEST(VEC, -, UnaryMinus)                           \
  UNARY_UNIT_TEST(VEC, std::exp, Exp)                           \
  UNARY_UNIT_TEST(VEC, std::log, Log)                           \
  UNARY_UNIT_TEST(VEC, std::log10, Log10)                       \
  UNARY_UNIT_TEST(VEC, std::sqrt, Sqrt)                         \
  UNARY_UNIT_TEST(VEC, std::sin, Sin)                           \
  UNARY_UNIT_TEST(VEC, std::cos, Cos)                           \
  UNARY_UNIT_TEST(VEC, std::tan, Tan)                           \
  UNARY_UNIT_TEST(VEC, std::sinh, Sinh)                         \
  UNARY_UNIT_TEST(VEC, std::cosh, Cosh)                         \
  UNARY_UNIT_TEST(VEC, std::tanh, Tanh)                         \
  UNARY_UNIT_TEST(VEC, std::asin, ASin)                         \
  UNARY_UNIT_TEST(VEC, std::acos, ACos)                         \
  UNARY_UNIT_TEST(VEC, std::atan, ATan)                         \
  UNARY_UNIT_TEST(VEC, std::asinh, ASinh)                       \
  UNARY_UNIT_TEST(VEC, std::acosh, ACosh)                       \
  UNARY_UNIT_TEST(VEC, std::atanh, ATanh)                       \
                                                                \
  BINARY_UNIT_TEST(VEC, +, Plus)                                \
  BINARY_UNIT_TEST(VEC, -, Minus)                               \
  BINARY_UNIT_TEST(VEC, *, Times)                               \
  BINARY_UNIT_TEST(VEC, /, Divide)                              \
                                                                \
  BINARYFUNC_UNIT_TEST(VEC, atan2, std::atan2, ATan2)           \
  BINARYFUNC_UNIT_TEST(VEC, pow, std::pow, Pow)                 \
  BINARYFUNC_UNIT_TEST(VEC, max, std::max, Max)                 \
  BINARYFUNC_UNIT_TEST(VEC, min, std::min, Min)                 \
                                                                \
  OPASSIGN_UNIT_TEST(VEC, +=, PlusEqual)                        \
  OPASSIGN_UNIT_TEST(VEC, -=, MinusEqual)                       \
  OPASSIGN_UNIT_TEST(VEC, *=, TimesEqual)                       \
  OPASSIGN_UNIT_TEST(VEC, /=, DivideEqual)                      \
                                                                \
  SAXPY_UNIT_TEST(VEC)

namespace DynamicVecTest {
  typedef Kokkos::Threads device_type;
  typedef Stokhos::DynamicStorage<int,double,device_type> storage_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  typedef UnitTestSetup<vec_type> UTS;
  UTS setup;
  VECTOR_UNIT_TESTS(DynamicVector)
}

namespace DynamicStridedVecTest {
  typedef Kokkos::Threads device_type;
  typedef Stokhos::DynamicStridedStorage<int,double,device_type> storage_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  typedef UnitTestSetup<vec_type> UTS;
  UTS setup;
  VECTOR_UNIT_TESTS(DynamicStridedVector)
}

namespace StaticVecTest {
  typedef Kokkos::Threads device_type;
  typedef Stokhos::StaticStorage<int,double,8,device_type> storage_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  typedef UnitTestSetup<vec_type> UTS;
  UTS setup;
  VECTOR_UNIT_TESTS(StaticVector)
}

namespace StaticFixedVecTest {
  typedef Kokkos::Threads device_type;
  typedef Stokhos::StaticFixedStorage<int,double,8,device_type> storage_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  typedef UnitTestSetup<vec_type> UTS;
  UTS setup;
  VECTOR_UNIT_TESTS(StaticFixedVector)
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
