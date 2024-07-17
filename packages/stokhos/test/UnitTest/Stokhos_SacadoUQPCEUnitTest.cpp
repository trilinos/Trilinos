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

#include "Stokhos.hpp"
#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

// Tests only run on the host, so use host execution space
typedef Kokkos::DefaultHostExecutionSpace execution_space;
typedef Stokhos::DynamicStorage<int,double,execution_space> storage_type;
typedef Sacado::UQ::PCE<storage_type> pce_type;

namespace Stokhos {

  template<class PCEType, class OrdinalType, class ValueType>
  bool comparePCEs(const PCEType& a1,
                   const std::string& a1_name,
                   const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>&a2,
                   const std::string& a2_name,
                   const ValueType& rel_tol, const ValueType& abs_tol,
                   Teuchos::FancyOStream& out)
  {
    bool success = true;

    out << "Comparing " << a1_name << " == " << a2_name << " ... ";

    const OrdinalType n = a1.size();

    // Compare sizes
    if (a2.size() != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == "
          << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }

    // Compare elements
    for( OrdinalType i = 0; i < n; ++i ) {
      ValueType nrm = std::sqrt(a2.basis()->norm_squared(i));
      ValueType err = std::abs(a1.coeff(i) - a2[i]) / nrm;
      ValueType tol =
        abs_tol + rel_tol*std::max(std::abs(a1.coeff(i)),std::abs(a2[i]))/nrm;
      if (err  > tol) {
        out
          <<"\nError, relErr("<<a1_name<<"["<<i<<"],"
          <<a2_name<<"["<<i<<"]) = relErr("<<a1.coeff(i)<<","<<a2[i]<<") = "
          <<err<<" <= tol = "<<tol<<": failed!\n";
        success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }
    else {
      out << std::endl
          << a1_name << " = " << a1 << std::endl
          << a2_name << " = " << a2 << std::endl;
    }

    return success;
  }

}

namespace SacadoPCEUnitTest {

  // Common setup for unit tests
  template <typename PCEType>
  struct UnitTestSetup {

    typedef PCEType pce_type;
    typedef typename pce_type::ordinal_type ordinal_type;
    typedef typename pce_type::value_type value_type;
    typedef typename pce_type::execution_space execution_space;
    typedef typename pce_type::cijk_type kokkos_cijk_type;
    typedef Stokhos::OrthogPolyApprox<ordinal_type,value_type> opa_type;
    value_type rtol, atol;
    value_type crtol, catol;
    ordinal_type sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > basis;
    Teuchos::RCP<Stokhos::Sparse3Tensor<ordinal_type,value_type> > Cijk;
    kokkos_cijk_type cijk;
    Teuchos::RCP< Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type,value_type> > exp;
    pce_type x, y, sin_x, cos_y, cx, u, u2, cu, cu2, sx, su, su2;
    opa_type x_opa, y_opa, sin_x_opa, cos_y_opa, cx_opa;
    value_type a;

    UnitTestSetup() {
      rtol = 1e-4;
      atol = 1e-5;
      crtol = 1e-12;
      catol = 1e-12;
      a = 3.1;
      const ordinal_type d = 2;
      const ordinal_type p = 7;

      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(d);
      for (ordinal_type i=0; i<d; i++)
        bases[i] =
          Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(p, true));
      basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases));

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();

      // Kokkos triple product tensor
      cijk = Stokhos::create_product_tensor<execution_space>(*basis, *Cijk);

      // Algebraic expansion
      exp = Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<ordinal_type,value_type>(basis, Cijk));

      // Quad expansion for initialization
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad =
        Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
      Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > quad_exp =
        Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

      // Create approximation
      x_opa.reset(basis);
      y_opa.reset(basis);
      sin_x_opa.reset(basis);
      cos_y_opa.reset(basis);
      cx_opa.reset(basis,1);
      x_opa.term(0, 0) = 1.0;
      y_opa.term(0, 0) = 2.0;
      cx_opa.term(0, 0) = a;
      for (int i=0; i<d; i++) {
        x_opa.term(i, 1) = 0.1;
        y_opa.term(i, 1) = 0.25;
      }
      quad_exp->sin(sin_x_opa, x_opa);
      quad_exp->cos(cos_y_opa, y_opa);

      // Create PCEs
      x.reset(cijk);
      y.reset(cijk);
      sin_x.reset(cijk);
      cos_y.reset(cijk);
      cx.reset(cijk, 1);
      x.load(x_opa.coeff());
      y.load(y_opa.coeff());
      sin_x.load(sin_x_opa.coeff());
      cos_y.load(cos_y_opa.coeff());
      cx.load(cx_opa.coeff());

      u.reset(cijk);
      u2.reset(cijk);
      cu.reset(cijk);
      cu2.reset(cijk, 1);
      sx.reset(cijk, d+1);
      su.reset(cijk, d+1);
      su2.reset(cijk, d+1);
      for (ordinal_type i=0; i<d; i++) {
        sx.fastAccessCoeff(i+1) = 0.0;
      }
    }
  };

  typedef UnitTestSetup<pce_type> UTS;

  TEUCHOS_UNIT_TEST( Stokhos_PCE, UMinus) {
    UTS setup;
    UTS::pce_type u = -setup.sin_x;
    UTS::opa_type u_opa(setup.basis);
    setup.exp->unaryMinus(u_opa, setup.sin_x_opa);
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",
                                   setup.rtol, setup.atol, out);
  }


#define UNARY_UNIT_TEST(OP)                                             \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, OP##_const) {                         \
    UTS setup;                                                          \
    UTS::pce_type u = OP(setup.cx);                                     \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->OP(u_opa, setup.cx_opa);                                 \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, OP##_resize) {                        \
    UTS setup;                                                          \
    UTS::pce_type u;                                                    \
    u = OP(setup.cx);                                                   \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->OP(u_opa, setup.cx_opa);                                 \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }

  UNARY_UNIT_TEST(exp)
  UNARY_UNIT_TEST(log)
  UNARY_UNIT_TEST(log10)
  UNARY_UNIT_TEST(sqrt)
  UNARY_UNIT_TEST(cbrt)
  UNARY_UNIT_TEST(sin)
  UNARY_UNIT_TEST(cos)
  UNARY_UNIT_TEST(tan)
  UNARY_UNIT_TEST(sinh)
  UNARY_UNIT_TEST(cosh)
  UNARY_UNIT_TEST(tanh)
  UNARY_UNIT_TEST(asin)
  UNARY_UNIT_TEST(acos)
  UNARY_UNIT_TEST(atan)
  // UNARY_UNIT_TEST(asinh)
  // UNARY_UNIT_TEST(acosh)
  // UNARY_UNIT_TEST(atanh)

#define BINARY_UNIT_TEST(OP, EXPOP)                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP) {                              \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type w = setup.cos_y;                                      \
    UTS::pce_type u = OP(v,w);                                          \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa, setup.cos_y_opa);          \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const) {                 \
    UTS setup;                                                          \
    UTS::pce_type w = setup.sin_x;                                      \
    UTS::pce_type u = OP(setup.a, w);                                   \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.a, setup.sin_x_opa);                  \
    success = Stokhos::comparePCEs(u, "u",  u_opa, "u_opa",             \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const) {                \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type u = OP(v, setup.a);                                   \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa, setup.a);                  \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_both_const) {                 \
    UTS setup;                                                          \
    UTS::pce_type u = OP(setup.cx, setup.cx);                           \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.cx_opa, setup.cx_opa);                \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const2) {                \
    UTS setup;                                                          \
    UTS::pce_type w = setup.sin_x;                                      \
    UTS::pce_type u = OP(setup.cx, w);                                  \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.cx_opa, setup.sin_x_opa);             \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const2) {               \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type u = OP(v, setup.cx);                                  \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa, setup.cx_opa);             \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_resize) {                     \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type w = setup.cos_y;                                      \
    UTS::pce_type u;                                                    \
    u = OP(v, w);                                                       \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa, setup.cos_y_opa);          \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const_resize) {          \
    UTS setup;                                                          \
    UTS::pce_type w = setup.sin_x;                                      \
    UTS::pce_type u;                                                    \
    u = OP(setup.a, w);                                                 \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.a, setup.sin_x_opa);                  \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const_resize) {         \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type u;                                                    \
    u = OP(v, setup.a);                                                 \
    UTS::opa_type u_opa(setup.basis);                                   \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa, setup.a);                  \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }

  BINARY_UNIT_TEST(operator+, plus)
  BINARY_UNIT_TEST(operator-, minus)
  BINARY_UNIT_TEST(operator*, times)
  BINARY_UNIT_TEST(operator/, divide)

#define OPASSIGN_UNIT_TEST(OP, EXPOP)                                   \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP) {                              \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type u = setup.cos_y;                                      \
    u OP v;                                                             \
    UTS::opa_type u_opa = setup.cos_y_opa;                              \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa);                           \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_const) {                      \
    UTS setup;                                                          \
    UTS::pce_type u = setup.sin_x;                                      \
    u OP setup.a;                                                       \
    UTS::opa_type u_opa = setup.sin_x_opa;                              \
    setup.exp->EXPOP(u_opa, setup.a);                                   \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_const2) {                     \
    UTS setup;                                                          \
    UTS::pce_type u = setup.sin_x;                                      \
    u OP setup.cx;                                                      \
    UTS::opa_type u_opa = setup.sin_x_opa;                              \
    setup.exp->EXPOP(u_opa, setup.cx_opa);                              \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_resize) {                     \
    UTS setup;                                                          \
    UTS::pce_type v = setup.sin_x;                                      \
    UTS::pce_type u = setup.a;                                          \
    u OP v;                                                             \
    UTS::opa_type u_opa = setup.cx_opa;                                 \
    setup.exp->EXPOP(u_opa, setup.sin_x_opa);                           \
    success = Stokhos::comparePCEs(u, "u", u_opa, "u_opa",              \
                                   setup.rtol, setup.atol, out);        \
  }

  OPASSIGN_UNIT_TEST(+=, plusEqual)
  OPASSIGN_UNIT_TEST(-=, minusEqual)
  OPASSIGN_UNIT_TEST(*=, timesEqual)
  OPASSIGN_UNIT_TEST(/=, divideEqual)

  TEUCHOS_UNIT_TEST( Stokhos_PCE, initializer_list_copy ) {
    UTS setup;
    UTS::pce_type u = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    UTS::opa_type v(setup.basis, 8);
    for (int i=0; i<8; i++)
      v[i] = i+1;
    success = comparePCEs(u, "u", v, "v",
                          setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, initializer_list_assign ) {
    UTS setup;
    UTS::pce_type u;
    u = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0 };
    UTS::opa_type v(setup.basis, 8);
    for (int i=0; i<8; i++)
      v[i] = i+1;
    success = comparePCEs(u, "u", v, "v",
                          setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, range_based_for ) {
    UTS setup;
    UTS::pce_type u;
    u.reset(UTS::pce_type::cijk_type(), 8);
    for (auto& z : u) { z = 3.0; }
    UTS::opa_type v(setup.basis, 8);
    for (int i=0; i<8; i++)
      v[i] = 3.0;
    success = comparePCEs(u, "u", v, "v",
                          setup.rtol, setup.atol, out);
  }
}
