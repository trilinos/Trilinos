// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
    
namespace SacadoPCEUnitTest {

  // Common setup for unit tests
  template <typename PCEType>
  struct UnitTestSetup {
    
    typedef PCEType pce_type;
    typedef Stokhos::OrthogPolyApprox<int,double> opa_type;
    double rtol, atol;
    double crtol, catol;
    int sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > exp;
    pce_type x, y, u, u2, cx, cu, cu2, sx, su, su2;
    double a;
    
    UnitTestSetup() {
      rtol = 1e-4;
      atol = 1e-5;
      crtol = 1e-12;
      catol = 1e-12;
      a = 3.1;
      const int d = 2;
      const int p = 7;
      
      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
      for (int i=0; i<d; i++)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p));
      basis =
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
      
      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
      
      // Quadrature expansion
      exp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));
      
      // Create approximation
      sz = basis->size();
      x.reset(exp);
      y.reset(exp);
      u.reset(exp); 
      u2.reset(exp);
      cx.reset(exp, 1);
      x.term(0, 0) = 1.0;
      cx.term(0, 0) = a;
      cu.reset(exp);
      cu2.reset(exp, 1);
      sx.reset(exp, d+1);
      su.reset(exp, d+1);
      su2.reset(exp, d+1);
      for (int i=0; i<d; i++) {
	x.term(i, 1) = 0.1;
	sx.term(i, 1) = 0.0;
      }
      y.term(0, 0) = 2.0;
      for (int i=0; i<d; i++)
	y.term(i, 1) = 0.25;
    }
  };

  typedef UnitTestSetup<pce_type> UTS;
  UTS setup;

  TEUCHOS_UNIT_TEST( Stokhos_PCE, UMinus) {
    UTS::pce_type v = std::sin(setup.x);
    UTS::pce_type u = -v;
    UTS::opa_type u_opa(setup.basis);
    setup.exp->unaryMinus(u_opa, v.getOrthogPolyApprox());
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u", 
				   u_opa, "u_opa", 
				   setup.rtol, setup.atol, out);
  }

#define UNARY_UNIT_TEST(OP)						\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, OP) {					\
    UTS::pce_type u = OP(setup.x);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->OP(u_opa, setup.x.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, OP##_const) {				\
    UTS::pce_type u = OP(setup.cx);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->OP(u_opa, setup.cx.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, OP##_resize) {			\
    UTS::pce_type u;							\
    u = OP(setup.x);							\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->OP(u_opa, setup.x.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }
 
  UNARY_UNIT_TEST(exp)
  UNARY_UNIT_TEST(log)
  UNARY_UNIT_TEST(log10)
  UNARY_UNIT_TEST(sqrt)
  UNARY_UNIT_TEST(sin)
  UNARY_UNIT_TEST(cos)
  UNARY_UNIT_TEST(tan)
  UNARY_UNIT_TEST(sinh)
  UNARY_UNIT_TEST(cosh)
  UNARY_UNIT_TEST(tanh)
  UNARY_UNIT_TEST(asin)
  UNARY_UNIT_TEST(acos)
  UNARY_UNIT_TEST(atan)
  UNARY_UNIT_TEST(asinh)
  UNARY_UNIT_TEST(acosh)
  UNARY_UNIT_TEST(atanh)

#define BINARY_UNIT_TEST(OP, EXPOP)					\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP) {				\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u = OP(v,w);						\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const) {			\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u = OP(setup.a, w);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.a,					\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const) {		\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u = OP(v, setup.a);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     setup.a);						\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_both_const) {			\
    UTS::pce_type u = OP(setup.cx, setup.cx);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.cx.getOrthogPolyApprox(),		\
  		     setup.cx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const2) {		\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u = OP(setup.cx, w);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.cx.getOrthogPolyApprox(),		\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const2) {		\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u = OP(v, setup.cx);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     setup.cx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_resize) {			\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u;							\
    u = OP(v, w);							\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_const_resize) {		\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u;							\
    u = OP(setup.a, w);							\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.a,					\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_const_resize) {		\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u;							\
    u = OP(v, setup.a);							\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     setup.a);						\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_short) {			\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u = OP(setup.sx, w);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox(),		\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_short) {		\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u = OP(v, setup.sx);					\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     setup.sx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_both_short) {			\
    UTS::pce_type u = OP(setup.sx, setup.sx);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox(),		\
  		     setup.sx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_short_resize) {		\
    UTS::pce_type w = std::cos(setup.x);				\
    UTS::pce_type u;							\
    u = OP(setup.sx, w);						\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox(),		\
  		     w.getOrthogPolyApprox());				\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_short_resize) {		\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u;							\
    u = OP(v, setup.sx);						\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox(),			\
  		     setup.sx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_short_right_const) {	\
    UTS::pce_type u = OP(setup.sx, setup.a);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox(),		\
  		     setup.a);						\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_short_left_const) {	\
    UTS::pce_type u = OP(setup.a, setup.sx);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.a,					\
  		     setup.sx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_left_short_right_const2) {	\
    UTS::pce_type u = OP(setup.sx, setup.cx);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox(),		\
  		     setup.cx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_right_short_left_const2) {	\
    UTS::pce_type u = OP(setup.cx, setup.sx);				\
    UTS::opa_type u_opa(setup.basis);					\
    setup.exp->EXPOP(u_opa, setup.cx.getOrthogPolyApprox(),		\
  		     setup.sx.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
  				   u_opa, "u_opa",			\
  				   setup.rtol, setup.atol, out);	\
  }

  BINARY_UNIT_TEST(operator+, plus)
  BINARY_UNIT_TEST(operator-, minus)
  BINARY_UNIT_TEST(operator*, times)
  BINARY_UNIT_TEST(operator/, divide)

#define OPASSIGN_UNIT_TEST(OP, EXPOP)					\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP) {				\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u = std::cos(setup.x);				\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP v;								\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_const) {			\
    UTS::pce_type u = std::cos(setup.x);				\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP setup.a;							\
    setup.exp->EXPOP(u_opa, setup.a);					\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_const2) {			\
    UTS::pce_type u = std::cos(setup.x);				\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP setup.cx;							\
    setup.exp->EXPOP(u_opa, setup.cx.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_resize) {			\
    UTS::pce_type v = std::sin(setup.x);				\
    UTS::pce_type u = setup.a;						\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP v;								\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_short) {			\
    UTS::pce_type u = std::cos(setup.x);				\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP setup.sx;							\
    setup.exp->EXPOP(u_opa, setup.sx.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_short_resize) {		\
    UTS::pce_type u = setup.sx;						\
    UTS::pce_type v = std::cos(setup.x);				\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP v;								\
    setup.exp->EXPOP(u_opa, v.getOrthogPolyApprox());			\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_short_const) {		\
    UTS::pce_type u = setup.sx;						\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP setup.a;							\
    setup.exp->EXPOP(u_opa, setup.a);					\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }									\
  TEUCHOS_UNIT_TEST( Stokhos_PCE, EXPOP##_short_const2) {		\
    UTS::pce_type u = setup.sx;						\
    UTS::opa_type u_opa = u.getOrthogPolyApprox();			\
    u OP setup.cx;							\
    setup.exp->EXPOP(u_opa, setup.cx.getOrthogPolyApprox());		\
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",	\
				   u_opa, "u_opa",			\
				   setup.rtol, setup.atol, out);	\
  }

  OPASSIGN_UNIT_TEST(+=, plusEqual)
  OPASSIGN_UNIT_TEST(-=, minusEqual)
  OPASSIGN_UNIT_TEST(*=, timesEqual)
  OPASSIGN_UNIT_TEST(/=, divideEqual)

  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::pce_type v = std::cos(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += v*w;
    setup.exp->times(t_opa, v.getOrthogPolyApprox(), w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_resize) {
    UTS::pce_type u = setup.cx;
    UTS::pce_type v = std::cos(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += v*w;
    setup.exp->times(t_opa, v.getOrthogPolyApprox(), w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_resize2) {
    UTS::pce_type u = setup.sx;
    UTS::pce_type v = std::cos(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += v*w;
    setup.exp->times(t_opa, v.getOrthogPolyApprox(), w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_const) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.a*w;
    setup.exp->times(t_opa, setup.a, w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_const2) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.cx*w;
    setup.exp->times(t_opa, setup.cx.getOrthogPolyApprox(), 
		     w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_short) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.sx*w;
    setup.exp->times(t_opa, setup.sx.getOrthogPolyApprox(), 
		     w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_short_const) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.sx*setup.a;
    setup.exp->times(t_opa, setup.sx.getOrthogPolyApprox(), 
		     setup.a);
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_short_const2) {
    UTS::pce_type u = std::sin(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.sx*setup.cx;
    setup.exp->times(t_opa, setup.sx.getOrthogPolyApprox(), 
		     setup.cx.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_resize_const) {
    UTS::pce_type u = setup.cx;
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.a*w;
    setup.exp->times(t_opa, setup.a, w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_resize_const2) {
    UTS::pce_type u = setup.cx;
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.cx*w;
    setup.exp->times(t_opa, setup.cx.getOrthogPolyApprox(), 
		     w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_short_resize_const) {
    UTS::pce_type u = setup.sx;
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.a*w;
    setup.exp->times(t_opa, setup.a, w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_PCE, saxpy_short_resize_const2) {
    UTS::pce_type u = setup.sx;
    UTS::pce_type w = std::exp(setup.x);
    UTS::opa_type u_opa = u.getOrthogPolyApprox();
    UTS::opa_type t_opa(setup.basis);
    u += setup.cx*w;
    setup.exp->times(t_opa, setup.cx.getOrthogPolyApprox(), 
		     w.getOrthogPolyApprox());
    setup.exp->plusEqual(u_opa, t_opa);
    success = Stokhos::comparePCEs(u.getOrthogPolyApprox(), "u",
				   u_opa, "u_opa",
				   setup.rtol, setup.atol, out);
  }
}
