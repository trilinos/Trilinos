// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_Complex.hpp"

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
  typedef typename vec_type::value_type value_type;

  double rtol, atol;
  double crtol, catol;
  int sz;
  vec_type x, y, cx;
  value_type a;

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

/**
 * @def UNARY_UNIT_TEST
 * Common series of test for any unary operator.
 */
#define UNARY_UNIT_TEST(VEC, SCALAR_T, OP, OPNAME, USING_OP)            \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME) {                        \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.x);                                      \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_OP                                                          \
      v.fastAccessCoeff(i) = OP(setup.x.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_const) {                \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.cx);                                     \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
    {                                                                   \
    USING_OP                                                            \
      v.fastAccessCoeff(i) = OP(setup.cx.fastAccessCoeff(0));           \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_resize) {               \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = OP(setup.x);                                                    \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
    USING_OP                                                            \
      v.fastAccessCoeff(i) = OP(setup.x.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define BINARY_UNIT_TEST(VEC, SCALAR_T, OP, OPNAME)                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME) {                        \
    UTS setup;                                                          \
    UTS::vec_type u = setup.x OP setup.y;                               \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const) {           \
    UTS setup;                                                          \
    UTS::vec_type u = setup.a OP setup.y;                               \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.a OP setup.y.fastAccessCoeff(i);     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const) {          \
    UTS setup;                                                          \
    UTS::vec_type u = setup.x OP setup.a ;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.a;                                                        \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_both_const) {           \
    UTS setup;                                                          \
    UTS::vec_type u = setup.cx OP setup.cx;                             \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) OP             \
        setup.cx.fastAccessCoeff(0);                                    \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const2) {          \
    UTS setup;                                                          \
    UTS::vec_type u = setup.cx OP setup.x;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) OP             \
        setup.x.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const2) {         \
    UTS setup;                                                          \
    UTS::vec_type u = setup.x OP setup.cx;                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.cx.fastAccessCoeff(0);                                    \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_resize) {               \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = setup.x OP setup.y;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const_resize) {    \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = setup.a OP setup.y;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.a OP                                 \
        setup.y.fastAccessCoeff(i);                                     \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const_resize) {   \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = setup.x OP setup.a;                                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.x.fastAccessCoeff(i) OP              \
        setup.a;                                                        \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, OP, SOP, USING_SOP, OPNAME) \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME) {                        \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.x,setup.y);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.y.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const) {           \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.a,setup.y);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.a,                               \
                                setup.y.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const) {          \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.x,setup.a);                              \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.a);                               \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_both_const) {           \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.cx,setup.cx);                            \
    UTS::vec_type v(1, 0.0);                                            \
    for (int i=0; i<v.size(); i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.cx.fastAccessCoeff(0),           \
                                 setup.cx.fastAccessCoeff(0));          \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const2) {          \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.cx,setup.x);                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.cx.fastAccessCoeff(0),           \
                                setup.x.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const2) {         \
    UTS setup;                                                          \
    UTS::vec_type u = OP(setup.x,setup.cx);                             \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.cx.fastAccessCoeff(0));           \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_resize) {               \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = OP(setup.x,setup.y);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.y.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_left_const_resize) {    \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = OP(setup.a,setup.y);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.a,                               \
                                setup.y.fastAccessCoeff(i));            \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_right_const_resize) {   \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = OP(setup.x,setup.a);                                            \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
    {                                                                   \
      USING_SOP                                                         \
      v.fastAccessCoeff(i) = SOP(setup.x.fastAccessCoeff(i),            \
                                setup.a);                               \
    }                                                                   \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define OPASSIGN_UNIT_TEST(VEC, SCALAR_T, OP, OPNAME)                   \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME) {                        \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.x;                                                       \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.x.fastAccessCoeff(i);               \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_const) {                \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.a;                                                       \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.a;                                  \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_const2) {               \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u OP setup.cx;                                                      \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) OP setup.cx.fastAccessCoeff(0);              \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, OPNAME##_resize) {               \
    UTS setup;                                                          \
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

#define SAXPY_UNIT_TEST(VEC, SCALAR_T)                                  \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, saxpy) {                         \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.x*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.x.fastAccessCoeff(i)*setup.y.fastAccessCoeff(i);          \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, saxpy_resize) {                  \
    UTS setup;                                                          \
    UTS::vec_type u = setup.cx;                                         \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    u += setup.x*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = setup.cx.fastAccessCoeff(0) +              \
        setup.x.fastAccessCoeff(i)*setup.y.fastAccessCoeff(i);          \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, saxpy_const) {                   \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.a*setup.y;                                               \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.a*setup.y.fastAccessCoeff(i);                             \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, saxpy_const2) {                  \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = std::sin(setup.x);                                \
    u += setup.cx*setup.y;                                              \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) +=                                           \
        setup.cx.fastAccessCoeff(0)*setup.y.fastAccessCoeff(i);         \
    success = compareVecs(u, "u",v, "v",                                \
                          setup.rtol, setup.atol, out);                 \
  }

#define TERNARY_UNIT_TEST(VEC, SCALAR_T)                                \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, ternay) {                        \
    UTS setup;                                                          \
    UTS::vec_type u = std::sin(setup.x);                                \
    UTS::vec_type v = -std::sin(setup.x);                               \
    u = u >= 0 ? -u : u;                                                 \
    success = compareVecs(u, "u", v, "v",                               \
                          setup.rtol, setup.atol, out);                 \
  }

/**
 * Default series of tests to perform for any type.
 * The set of operators should make sense for both real and complex types.
 */
#define VECTOR_UNIT_TESTS_ANY_TYPE(VEC,SCALAR_T)                        \
  UNARY_UNIT_TEST(VEC, SCALAR_T, +    , UnaryPlus ,)                    \
  UNARY_UNIT_TEST(VEC, SCALAR_T, -    , UnaryMinus,)                    \
  UNARY_UNIT_TEST(VEC, SCALAR_T, exp  , Exp       , using std::exp;  )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, log  , Log       , using std::log;  )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, log10, Log10     , using std::log10;)  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, sqrt , Sqrt      , using std::sqrt; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, sin  , Sin       , using std::sin;  )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, cos  , Cos       , using std::cos;  )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, tan  , Tan       , using std::tan;  )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, sinh , Sinh      , using std::sinh; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, cosh , Cosh      , using std::cosh; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, tanh , Tanh      , using std::tanh; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, asin , ASin      , using std::asin; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, acos , ACos      , using std::acos; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, atan , ATan      , using std::atan; )  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, asinh, ASinh     , using std::asinh;)  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, acosh, ACosh     , using std::acosh;)  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, atanh, ATanh     , using std::atanh;)  \
                                                                \
  BINARY_UNIT_TEST(VEC, SCALAR_T, +, Plus)                              \
  BINARY_UNIT_TEST(VEC, SCALAR_T, -, Minus)                             \
  BINARY_UNIT_TEST(VEC, SCALAR_T, *, Times)                             \
  BINARY_UNIT_TEST(VEC, SCALAR_T, /, Divide)                            \
                                                                \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, pow, pow, using std::pow;, Pow)   \
                                                                \
  OPASSIGN_UNIT_TEST(VEC, SCALAR_T, +=, PlusEqual)                      \
  OPASSIGN_UNIT_TEST(VEC, SCALAR_T, -=, MinusEqual)                     \
  OPASSIGN_UNIT_TEST(VEC, SCALAR_T, *=, TimesEqual)                     \
  OPASSIGN_UNIT_TEST(VEC, SCALAR_T, /=, DivideEqual)                    \
                                                                \
  SAXPY_UNIT_TEST(VEC, SCALAR_T)                                        \
                                                                \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, initializer_list_constructor ) { \
    UTS setup;                                                          \
    UTS::vec_type u{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };          \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = i+1;                                       \
    success = compareVecs(u, "u", v, "v",                               \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, initializer_list_copy ) {        \
    UTS setup;                                                          \
    UTS::vec_type u = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };       \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = i+1;                                       \
    success = compareVecs(u, "u", v, "v",                               \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, initializer_list_assign ) {      \
    UTS setup;                                                          \
    UTS::vec_type u;                                                    \
    u = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };                     \
    UTS::vec_type v(setup.sz, 0.0);                                     \
    for (int i=0; i<setup.sz; i++)                                      \
      v.fastAccessCoeff(i) = i+1;                                       \
    success = compareVecs(u, "u", v, "v",                               \
                          setup.rtol, setup.atol, out);                 \
  }                                                                     \
  TEUCHOS_UNIT_TEST( VEC##_##SCALAR_T, range_based_for ) {              \
    UTS setup;                                                          \
    UTS::vec_type u(setup.sz, 0.0);                                     \
    for (auto& z : u) { z = 3.0; }                                      \
    UTS::vec_type v(setup.sz, 3.0);                                     \
    success = compareVecs(u, "u", v, "v",                               \
                          setup.rtol, setup.atol, out);                 \
  }

  /**
   * Default series of tests to perform for SFS for any value type.
   */
  #define VECTOR_UNIT_TESTS_SFS_ANY_VALUE_TYPE(SCALAR_T)                                     \
  TEUCHOS_UNIT_TEST( StaticFixedVector##_##SCALAR_T, initializer_list_constructor_partial ) {\
    UTS setup;                                                                               \
    UTS::vec_type u{ 1.0};                                                                   \
    UTS::vec_type v(setup.sz, 1.0);                                                          \
    success = compareVecs(u, "u", v, "v",                                                    \
                          setup.rtol, setup.atol, out);                                      \
  }                                                                                          \
  TEUCHOS_UNIT_TEST( StaticFixedVector##_##SCALAR_T, initializer_list_constructor_empty  ) { \
    UTS setup;                                                                               \
    UTS::vec_type u{ std::initializer_list<typename UTS::value_type>()};                     \
    UTS::vec_type v(setup.sz, 0.0);                                                          \
    success = compareVecs(u, "u", v, "v",                                                    \
                          setup.rtol, setup.atol, out);                                      \
  }
/**
 * Series of tests to run for complex type.
 * It will run the series of tests for any type.
 */
#define VECTOR_UNIT_TESTS_COMPLEX_TYPE(VEC,SCALAR_T)                               \
  /* Run the series of tests for any type. */                                      \
  VECTOR_UNIT_TESTS_ANY_TYPE(VEC,SCALAR_T)

/**
 * Series of tests to run for real type.
 * It will run the series of tests for any type, as well as tests that don't work or
 * don't make sense for complex type.
 */
#define VECTOR_UNIT_TESTS_REAL_TYPE(VEC,SCALAR_T)                                  \
  /* Run the series of tests for any type. */                                      \
  VECTOR_UNIT_TESTS_ANY_TYPE(VEC,SCALAR_T)                                         \
  /* Operator cbrt not supported for complex type but supported for real type. */  \
  UNARY_UNIT_TEST(VEC, SCALAR_T, cbrt , Cbrt, using std::cbrt;)                    \
  /* Operator atan2 not supported for complex type but supported for real type. */ \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, atan2, atan2, using std::atan2;, ATan2)      \
  /* Operators min and max are not supported for complex type. */                  \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, max  , max  , using std::max  ;, Max)        \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, min  , min  , using std::min  ;, Min)        \
  /* Operators fmin and fmax are not supported for complex type. */                \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, fmax , fmax , using std::fmax ;, FMax)       \
  BINARYFUNC_UNIT_TEST(VEC, SCALAR_T, fmin , fmin , using std::fmin ;, FMin)       \
  /* Ternary test uses 'operator<' that is not defined for complex type. */        \
  TERNARY_UNIT_TEST(VEC, SCALAR_T)

/**
 * Common test structure for dynamic storage.
 */
#define TEST_DYNAMIC_STORAGE(__storage_type__,__vec_type__,__scalar_type__,__macro_for_tests__)\
  typedef Kokkos::DefaultExecutionSpace execution_space;                                       \
  typedef Stokhos::__storage_type__<int,__scalar_type__,execution_space> storage_type;         \
  typedef Sacado::MP::Vector<storage_type> vec_type;                                           \
  typedef UnitTestSetup<vec_type> UTS;                                                         \
  __macro_for_tests__(__vec_type__,__scalar_type__)

namespace DynamicVecTest
{
  TEST_DYNAMIC_STORAGE(DynamicStorage, DynamicVector, double, VECTOR_UNIT_TESTS_REAL_TYPE)
}

namespace DynamicStridedVecTest
{
  TEST_DYNAMIC_STORAGE(DynamicStridedStorage, DynamicStridedVector, double, VECTOR_UNIT_TESTS_REAL_TYPE)
}

/**
 * Common test structure for static storage.
 */
#define TEST_STATIC_STORAGE(__storage_type__,__vec_type__,__scalar_type__,__scalar_type_name__,__storage_size__,__macro_for_tests__) \
    typedef ::Kokkos::DefaultExecutionSpace execution_space;                                                                         \
    typedef ::Stokhos::__storage_type__<int,__scalar_type__,__storage_size__,execution_space> storage_type;                          \
    typedef ::Sacado::MP::Vector<storage_type> vec_type;                                                                             \
    typedef UnitTestSetup<vec_type> UTS;                                                                                             \
    __macro_for_tests__(__vec_type__,__scalar_type_name__)

namespace StaticVecTest
{
  TEST_STATIC_STORAGE(StaticStorage, StaticVector, double, double, 8, VECTOR_UNIT_TESTS_REAL_TYPE)
}

/**
 * Common test structure for static fixed storage.
 */
#define TEST_STATIC_FIXED_STORAGE(__storage_type__,__vec_type__,__scalar_type__,__scalar_type_name__,__storage_size__,__macro_for_tests__) \
    TEST_STATIC_STORAGE(__storage_type__,__vec_type__,__scalar_type__,__scalar_type_name__,__storage_size__,__macro_for_tests__)           \
    VECTOR_UNIT_TESTS_SFS_ANY_VALUE_TYPE(__scalar_type_name__)


namespace StaticFixedVecTest
{
  namespace Double        {TEST_STATIC_FIXED_STORAGE(StaticFixedStorage, StaticFixedVector,                   double ,                double, 8, VECTOR_UNIT_TESTS_REAL_TYPE   )}

// Skip std::complex when compiling with CUDA, because std::complex isn't supported in that case.
// Note that even though the tests aren't run on the device, nvcc still complains that __device__ code functions are called
// from __host__ code (or vice versa).
#if !defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_HIP)
  namespace Complex_std   {TEST_STATIC_FIXED_STORAGE(StaticFixedStorage, StaticFixedVector,   std   ::complex<double>,    std_complex_double, 8, VECTOR_UNIT_TESTS_COMPLEX_TYPE)}
#endif

  // Always test for Kokkos::complex because it is always shipped as part of Kokkos, whatever the space.
  namespace Complex_Kokkos{TEST_STATIC_FIXED_STORAGE(StaticFixedStorage, StaticFixedVector, ::Kokkos::complex<double>, kokkos_complex_double, 8, VECTOR_UNIT_TESTS_COMPLEX_TYPE)}
}
