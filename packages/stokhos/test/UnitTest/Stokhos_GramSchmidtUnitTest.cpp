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
#include "Stokhos_UnitTestHelpers.hpp"
#include "Stokhos_StieltjesPCEBasis.hpp"
#include "Stokhos_GramSchmidtBasis.hpp"
#include "Stokhos_UserDefinedQuadrature.hpp"

// Quadrature functor to be passed into quadrature expansion for mapping
// from Gram-Schmidt basis back to original PCE
struct gram_schmidt_pce_unary_quad_func {
  gram_schmidt_pce_unary_quad_func(
    const Stokhos::OrthogPolyApprox<int,double>& pce_,
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
struct gram_schmidt_pce_binary_quad_func {
  gram_schmidt_pce_binary_quad_func(
    const Stokhos::OrthogPolyApprox<int,double>& pce_,
    const Stokhos::OrthogPolyBasis<int,double>& basis_) :
    pce(pce_), basis(basis_), vec(2) {}
  
  double operator() (const double& a, const double& b) const {
    vec[0] = a;
    vec[1] = b;
    return pce.evaluate(vec);
  }
  const Stokhos::OrthogPolyApprox<int,double>& pce;
  const Stokhos::OrthogPolyBasis<int,double>& basis;
  mutable Teuchos::Array<double> vec;
};

// Class encapsulating setup of the Gram-Schmidt basis for a given PCE
template <typename OrdinalType, typename ValueType>
struct GramSchmidt_PCE_Setup {
  ValueType rtol, atol;
  OrdinalType sz, gs_sz;
  Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
  Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > exp;
  Teuchos::RCP<const  Stokhos::GramSchmidtBasis<OrdinalType,ValueType> > gs_basis;
  Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > gs_quad;
  Stokhos::OrthogPolyApprox<OrdinalType,ValueType> u, v, w, u_gs, v_gs, w_gs;
  
  GramSchmidt_PCE_Setup()
  {
    rtol = 1e-7;
    atol = 1e-12;
    const OrdinalType d = 3;
    const OrdinalType p = 5;
    
    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
    for (OrdinalType i=0; i<d; i++)
      bases[i] = 
	Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
    basis =
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases, 1e-15));
    
    // Create approximation
    sz = basis->size();
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> x(basis);
    for (OrdinalType i=0; i<d; i++)
      x.term(i, 1) = 1.0;
    
    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis));

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    
    // Quadrature expansion
    Teuchos::RCP<Teuchos::ParameterList> exp_params =
      Teuchos::rcp(new Teuchos::ParameterList);
    exp_params->set("Use Quadrature for Times", true);
    exp = Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, quad, exp_params));
    
    // Compute PCE via quadrature expansion
    u.reset(basis);
    v.reset(basis);
    w.reset(basis);
    exp->sin(u,x);
    exp->exp(v,x);
    exp->times(w,u,v);
    
    // Compute Stieltjes basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > st_bases(2);
    st_bases[0] = Teuchos::rcp(new Stokhos::StieltjesPCEBasis<OrdinalType,ValueType>(p, Teuchos::rcp(&u,false), quad, true));
    st_bases[1] = Teuchos::rcp(new Stokhos::StieltjesPCEBasis<OrdinalType,ValueType>(p, Teuchos::rcp(&v,false), quad, true));
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<OrdinalType,ValueType> > st_basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(st_bases, 1e-15));
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> u_st(st_basis), v_st(st_basis);
    u_st.term(0, 0) = u.mean();
    u_st.term(0, 1) = 1.0;
    v_st.term(0, 0) = v.mean();
    v_st.term(1, 1) = 1.0;

    // Use Gram-Schmidt to orthogonalize Stieltjes basis of u and v
    Teuchos::Array<ValueType> st_points_0;
    Teuchos::Array<ValueType> st_weights_0;
    Teuchos::Array< Teuchos::Array<ValueType> > st_values_0;
    st_bases[0]->getQuadPoints(p+1, st_points_0, st_weights_0, st_values_0);
    Teuchos::Array<ValueType> st_points_1;
    Teuchos::Array<ValueType> st_weights_1;
    Teuchos::Array< Teuchos::Array<ValueType> > st_values_1;
    st_bases[1]->getQuadPoints(p+1, st_points_1, st_weights_1, st_values_1);
    Teuchos::RCP< Teuchos::Array< Teuchos::Array<ValueType> > > st_points = 
      Teuchos::rcp(new Teuchos::Array< Teuchos::Array<ValueType> >(st_points_0.size()));
    for (int i=0; i<st_points_0.size(); i++) {
      (*st_points)[i].resize(2);
      (*st_points)[i][0] = st_points_0[i];
      (*st_points)[i][1] = st_points_1[i];
    }
    Teuchos::RCP< Teuchos::Array<ValueType> > st_weights = 
      Teuchos::rcp(new Teuchos::Array<ValueType>);
    *st_weights = st_weights_0;
    
    gs_basis = Teuchos::rcp(new Stokhos::GramSchmidtBasis<OrdinalType,ValueType>(st_basis, *st_points, *st_weights, 1e-15));
    gs_sz = gs_basis->size();

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk =
      gs_basis->computeTripleProductTensor();

    // Create quadrature for Gram-Schmidt basis using quad points and 
    // and weights from original basis mapped to Stieljtes basis
    Teuchos::RCP< const Teuchos::Array< Teuchos::Array<ValueType> > > points =
      st_points;
    Teuchos::RCP< const Teuchos::Array<ValueType> > weights = st_weights;
    gs_quad = Teuchos::rcp(new Stokhos::UserDefinedQuadrature<OrdinalType,ValueType>(gs_basis, points, weights));
    
    // Gram-Schmidt quadrature expansion
    Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> gs_exp(gs_basis, 
								   gs_Cijk,
								   gs_quad,
								   exp_params);
    
    u_gs.reset(gs_basis);
    v_gs.reset(gs_basis);
    w_gs.reset(gs_basis);
    
    // Map expansion in Stieltjes basis to Gram-Schmidt basis
    gs_basis->transformCoeffs(u_st.coeff(), u_gs.coeff());
    gs_basis->transformCoeffs(v_st.coeff(), v_gs.coeff());
    
    // Compute w_gs = u_gs*v_gs in Gram-Schmidt basis
    gs_exp.times(w_gs, u_gs, v_gs);
  }
  
};

namespace GramSchmidtTest {
  GramSchmidt_PCE_Setup<int,double> setup;

  // Tests mapping from Gram-Schmidt basis to original is correct
  TEUCHOS_UNIT_TEST( Stokhos_GramSchmidtBasis, Map ) {
   Stokhos::OrthogPolyApprox<int,double> u2(setup.basis);
   gram_schmidt_pce_binary_quad_func quad_func(setup.u_gs, *setup.gs_basis);
   setup.exp->binary_op(quad_func, u2, setup.u, setup.v);
   success = comparePCEs(setup.u, "u", u2, "u2", setup.rtol, setup.atol, out);
  }

  // Tests Gram-Schmidt basis is orthogonal
  TEUCHOS_UNIT_TEST( Stokhos_GramSchmidtBasis, Orthog ) {
    const Teuchos::Array<double>& norms = setup.gs_basis->norm_squared();
    const Teuchos::Array<double>& weights = setup.gs_quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& values = 
      setup.gs_quad->getBasisAtQuadPoints();
    Teuchos::SerialDenseMatrix<int,double> mat(setup.gs_sz, setup.gs_sz);
    for (int i=0; i<setup.gs_sz; i++) {
      for (int j=0; j<setup.gs_sz; j++) {
	for (int k=0; k<weights.size(); k++)
	  mat(i,j) += weights[k]*values[k][i]*values[k][j];
	mat(i,j) /= std::sqrt(norms[i]*norms[j]);
      }
      mat(i,i) -= 1.0;
    }
    double tol = 1e-4;
    success = mat.normInf() < tol;
    if (!success) {
      out << "\n Error, mat.normInf() < tol = " << mat.normInf() 
	  << " < " << tol << ": failed!\n";
      out << "mat = " << printMat(mat) << std::endl;
    }
  }

  // Tests PCE computed from Gram-Schmidt basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_GramSchmidtBasis, PCE ) {
    Stokhos::OrthogPolyApprox<int,double> w2(setup.basis);
    gram_schmidt_pce_binary_quad_func quad_func(setup.w_gs, *setup.gs_basis);
    setup.exp->binary_op(quad_func, w2, setup.u, setup.v);
    success = comparePCEs(setup.w, "w", w2, "w2", setup.rtol, setup.atol, out);
  }

  // Tests mean computed from Gram-Schmidt basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_GramSchmidtBasis, Mean ) {
    success = Teuchos::testRelErr("w.mean()", setup.w.mean(), 
				  "w_gs.mean()", setup.w_gs.mean(),
				  "rtol", setup.rtol, 
				  "rtol", setup.rtol, 
				  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

  // Tests mean standard deviation from Gram-Schmidt basis is same as original
  TEUCHOS_UNIT_TEST( Stokhos_GramSchmidtBasis, StandardDeviation ) {
    success = Teuchos::testRelErr("w.standard_deviation()", 
				  setup.w.standard_deviation(), 
				  "w_gs.standard_devaition()", 
				  setup.w_gs.standard_deviation(),
				  "rtol", 1e-3, 
				  "rtol", 1e-3, 
				  Teuchos::Ptr<std::ostream>(out.getOStream().get()));
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
