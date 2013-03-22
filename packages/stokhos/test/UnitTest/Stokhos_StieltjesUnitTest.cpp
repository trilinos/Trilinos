// $Id$ 
// $Source$ 
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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"

// Quadrature functor to be passed into quadrature expansion for mapping
// from Stieltjes basis back to original PCE
struct stieltjes_pce_quad_func {
  stieltjes_pce_quad_func(const Stokhos::OrthogPolyApprox<int,double>& pce_,
			  const Stokhos::OrthogPolyBasis<int,double>& basis_) :
    pce(pce_), basis(basis_), vec(1) {}
  
  double operator() (const double& a) const {
    vec[0] = a;
    return pce.evaluate(vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable Teuchos::Array<double> vec;
};

// Class encapsulating setup of the Stieltjes basis for a given PCE
// u = Func(x), where Func is specified by a template-parameter
template <typename Func>
struct Stieltjes_PCE_Setup {
  typedef typename Func::OrdinalType OrdinalType;
  typedef typename Func::ValueType ValueType;
  ValueType rtol, atol;
  Func func;
  bool use_pce_quad_points;
  OrdinalType sz, st_sz;
  Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
  Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > exp;
  Teuchos::RCP<const Stokhos::StieltjesPCEBasis<OrdinalType,ValueType> > st_1d_basis;
  Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > st_basis;
  Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > st_quad;
  Stokhos::OrthogPolyApprox<OrdinalType,ValueType> u, v, u_st, v_st;
  
  Stieltjes_PCE_Setup(bool use_pce_quad_points_) : 
    func(), use_pce_quad_points(use_pce_quad_points_) 
  {
    rtol = 1e-8;
    atol = 1e-12;
    const OrdinalType d = 3;
    const OrdinalType p = 5;
    
    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
    for (OrdinalType i=0; i<d; i++)
      bases[i] = 
	Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
    basis =
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases));
    
    // Create approximation
    sz = basis->size();
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> x(basis);
    for (OrdinalType i=0; i<d; i++)
      x.term(i, 1) = 1.0;
    
    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis, 4*p));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    
    // Quadrature expansion
    exp = Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, quad));
    
    // Compute PCE via quadrature expansion
    u.reset(basis);
    v.reset(basis);
    func.eval(*exp, x, u);
    exp->times(v,u,u);
    
    // Compute Stieltjes basis
    st_1d_basis = 
      Teuchos::rcp(new Stokhos::StieltjesPCEBasis<OrdinalType,ValueType>(
		     p, Teuchos::rcp(&u,false), quad, use_pce_quad_points));
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > st_bases(1);
    st_bases[0] = st_1d_basis;
    st_basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(st_bases, 1e-15));
    st_sz = st_basis->size();
    u_st.reset(st_basis);
    v_st.reset(st_basis);
    u_st[0] = u.mean();
    u_st[1] = 1.0;
    
    // Tensor product quadrature
    st_quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(st_basis));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > st_Cijk =
      st_basis->computeTripleProductTensor();
    
    // Quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> st_exp(st_basis, 
								   st_Cijk,
								   st_quad);
    
    st_exp.times(v_st, u_st, u_st);
  }
  
};

//
// Stieltjes tests based on expansion of u = cos(x) where x is a U([-1,1])
// random variable
//
namespace StieltjesCosTest {
  
  template <typename Ordinal_Type, typename Value_Type>
  struct Stieltjes_Cos_Func {
    typedef Ordinal_Type OrdinalType;
    typedef Value_Type ValueType;
    static const bool is_even = true;
    void 
    eval(Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>& exp,
	 const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& x, 
	 Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& u) {
      exp.cos(u,x);
    }
  };
  Stieltjes_PCE_Setup< Stieltjes_Cos_Func<int,double> > setup(true);

  // Tests mapping from Stieltjes basis to original is correct
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, CosMap ) {
   Stokhos::OrthogPolyApprox<int,double> u2(setup.basis);
   setup.st_1d_basis->transformCoeffsFromStieltjes(setup.u_st.coeff(), 
						   u2.coeff());
   success = Stokhos::comparePCEs(setup.u, "u", u2, "u2", setup.rtol, setup.atol, out);
  }

  // Tests Stieltjes basis is orthogonal
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, CosOrthog ) {
    const Teuchos::Array<double>& norms = setup.st_1d_basis->norm_squared();
    const Teuchos::Array<double>& weights = setup.st_quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& values = 
      setup.st_quad->getBasisAtQuadPoints();
    Teuchos::SerialDenseMatrix<int,double> mat(setup.st_sz, 
					       setup.st_sz);
    for (int i=0; i<setup.st_sz; i++) {
      for (int j=0; j<setup.st_sz; j++) {
	for (unsigned int k=0; k<weights.size(); k++)
	  mat(i,j) += weights[k]*values[k][i]*values[k][j];
	mat(i,j) /= std::sqrt(norms[i]*norms[j]);
      }
      mat(i,i) -= 1.0;
    }
    success = mat.normInf() < setup.atol;
    if (!success) {
      out << "\n Error, mat.normInf() < atol = " << mat.normInf() 
	  << " < " << setup.atol << ": failed!\n";
      out << "mat = " << mat << std::endl;
    }
  }

  // Tests PCE computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, CosPCE ) {
    Stokhos::OrthogPolyApprox<int,double> v2(setup.basis);
    stieltjes_pce_quad_func quad_func(setup.v_st, *setup.st_basis);
    setup.exp->unary_op(quad_func, v2, setup.u);
    success = comparePCEs(setup.v, "v", v2, "v2", setup.rtol, setup.atol, out);
  }

  // Tests mean computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, CosMean ) {
    success = Teuchos::testRelErr("v.mean()", setup.v.mean(), 
				  "v_st.mean()", setup.v_st.mean(),
				  "rtol", setup.rtol, 
				  "rtol", setup.rtol, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));

  }

  // Tests mean standard deviation from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, CosStandardDeviation ) {
    success = Teuchos::testRelErr("v.standard_deviation()", 
				  setup.v.standard_deviation(), 
				  "v_st.standard_devaition()", 
				  setup.v_st.standard_deviation(),
				  "rtol", 1e-1, 
				  "rtol", 1e-1, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

}
  
//
// Stieltjes tests based on expansion of u = sin(x) where x is a U([-1,1])
// random variable
//
namespace StieltjesSinTest {

  template <typename Ordinal_Type, typename Value_Type>
  struct Stieltjes_Sin_Func {
    typedef Ordinal_Type OrdinalType;
    typedef Value_Type ValueType;
    static const bool is_even = false;
    void 
    eval(Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>& exp,
	 const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& x, 
	 Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& u) {
      exp.sin(u,x);
    }
  };
  Stieltjes_PCE_Setup< Stieltjes_Sin_Func<int,double> > setup(true);

  // Tests mapping from Stieltjes basis to original is correct
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, SinMap ) {
   Stokhos::OrthogPolyApprox<int,double> u2(setup.basis);
   setup.st_1d_basis->transformCoeffsFromStieltjes(setup.u_st.coeff(), 
						   u2.coeff());
   success = Stokhos::comparePCEs(setup.u, "u", u2, "u2", setup.rtol, setup.atol, out);
  }

  // Tests Stieltjes basis is orthogonal
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, SinOrthog ) {
    const Teuchos::Array<double>& norms = setup.st_1d_basis->norm_squared();
    const Teuchos::Array<double>& weights = setup.st_quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& values = 
      setup.st_quad->getBasisAtQuadPoints();
    Teuchos::SerialDenseMatrix<int,double> mat(setup.st_sz, 
					       setup.st_sz);
    for (int i=0; i<setup.st_sz; i++) {
      for (int j=0; j<setup.st_sz; j++) {
	for (unsigned int k=0; k<weights.size(); k++)
	  mat(i,j) += weights[k]*values[k][i]*values[k][j];
	mat(i,j) /= std::sqrt(norms[i]*norms[j]);
      }
      mat(i,i) -= 1.0;
    }
    success = mat.normInf() < setup.atol;
    if (!success) {
      out << "\n Error, mat.normInf() < atol = " << mat.normInf() 
	  << " < " << setup.atol << ": failed!\n";
      out << "mat = " << mat << std::endl;
    }
  }

  // Tests PCE computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, SinPCE ) {
    Stokhos::OrthogPolyApprox<int,double> v2(setup.basis);
    stieltjes_pce_quad_func quad_func(setup.v_st, *setup.st_basis);
    setup.exp->unary_op(quad_func, v2, setup.u);
    success = comparePCEs(setup.v, "v", v2, "v2", setup.rtol, setup.atol, out);
  }

  // Tests mean computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, SinMean ) {
    success = Teuchos::testRelErr("v.mean()", setup.v.mean(), 
				  "v_st.mean()", setup.v_st.mean(),
				  "rtol", setup.rtol, 
				  "rtol", setup.rtol, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

  // Tests mean standard deviation from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, SinStandardDeviation ) {
    success = Teuchos::testRelErr("v.standard_deviation()", 
				  setup.v.standard_deviation(), 
				  "v_st.standard_devaition()", 
				  setup.v_st.standard_deviation(),
				  "rtol", 1e-1, 
				  "rtol", 1e-1, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

}

//
// Stieltjes tests based on expansion of u = exp(x) where x is a U([-1,1])
// random variable.  For this test we don't use the PCE quad points and 
// instead use those generated for the Stieltjes basis
//
namespace StieltjesExpTest {

  template <typename Ordinal_Type, typename Value_Type>
  struct Stieltjes_Exp_Func {
    typedef Ordinal_Type OrdinalType;
    typedef Value_Type ValueType;
    static const bool is_even = false;
    void 
    eval(Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>& exp,
	 const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& x, 
	 Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& u) {
      exp.exp(u,x);
    }
  };
  Stieltjes_PCE_Setup< Stieltjes_Exp_Func<int,double> > setup(false);

  // Tests mapping from Stieltjes basis to original is correct
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, ExpMap ) {
   Stokhos::OrthogPolyApprox<int,double> u2(setup.basis);
   setup.st_1d_basis->transformCoeffsFromStieltjes(setup.u_st.coeff(), 
						   u2.coeff());
   success = Stokhos::comparePCEs(setup.u, "u", u2, "u2", setup.rtol, setup.atol, out);
  }

  // Tests Stieltjes basis is orthogonal
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, ExpOrthog ) {
    const Teuchos::Array<double>& norms = setup.st_1d_basis->norm_squared();
    const Teuchos::Array<double>& weights = setup.st_quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& values = 
      setup.st_quad->getBasisAtQuadPoints();
    Teuchos::SerialDenseMatrix<int,double> mat(setup.st_sz, 
					       setup.st_sz);
    for (int i=0; i<setup.st_sz; i++) {
      for (int j=0; j<setup.st_sz; j++) {
	for (unsigned int k=0; k<weights.size(); k++)
	  mat(i,j) += weights[k]*values[k][i]*values[k][j];
	mat(i,j) /= std::sqrt(norms[i]*norms[j]);
      }
      mat(i,i) -= 1.0;
    }
    success = mat.normInf() < setup.atol;
    if (!success) {
      out << "\n Error, mat.normInf() < atol = " << mat.normInf() 
	  << " < " << setup.atol << ": failed!\n";
      out << "mat = " << mat << std::endl;
    }
  }

  // Tests PCE computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, ExpPCE ) {
    Stokhos::OrthogPolyApprox<int,double> v2(setup.basis);
    stieltjes_pce_quad_func quad_func(setup.v_st, *setup.st_basis);
    setup.exp->unary_op(quad_func, v2, setup.u);
    success = comparePCEs(setup.v, "v", v2, "v2", setup.rtol, setup.atol, out);
  }

  // Tests mean computed from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, ExpMean ) {
    success = Teuchos::testRelErr("v.mean()", setup.v.mean(), 
				  "v_st.mean()", setup.v_st.mean(),
				  "rtol", setup.rtol, 
				  "rtol", setup.rtol, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

  // Tests mean standard deviation from Stieltjes basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_StieltjesPCEBasis, ExpStandardDeviation ) {
    success = Teuchos::testRelErr("v.standard_deviation()", 
				  setup.v.standard_deviation(), 
				  "v_st.standard_devaition()", 
				  setup.v_st.standard_deviation(),
				  "rtol", 1e-1, 
				  "rtol", 1e-1, 
                                  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
