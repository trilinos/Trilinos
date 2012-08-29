// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"
#include "Stokhos_Sacado.hpp"
#include "Stokhos_ReducedBasisFactory.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

// storage options
enum Basis_Reduction_Method { LANCZOS, MONOMIAL_GS, LANCZOS_GS };
static const int num_basis_reduction_method = 3;
static const Basis_Reduction_Method basis_reduction_method_values[] = { 
  LANCZOS, MONOMIAL_GS, LANCZOS_GS };
static const char *basis_reduction_method_names[] = { 
  "Lanczos", "Monomial-GS", "Lanczos-GS" };

typedef Stokhos::LegendreBasis<int,double> basis_type;
typedef Sacado::ETPCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_type;

template <class ScalarType>
inline ScalarType f(const Teuchos::Array<ScalarType>& x,
		    double a, double b) {
  ScalarType y = a;
  int n = x.size();
  for (int i=0; i<n; i++)
    y += x[i] / (i+1.0);
  y = b + 1.0/y;
  return y;
}

template <class ScalarType>
inline ScalarType g(const Teuchos::Array<ScalarType>& x,
		    const ScalarType& y) {
  ScalarType z = y;
  //ScalarType z = 0.0;
  int n = x.size();
  for (int i=0; i<n; i++)
    z += x[i];
  z = std::exp(z);
  //z = y*z;
  return z;
}

double coeff_error(const pce_type& z, const pce_type& z2) {
  double err_z = 0.0;
  for (int i=0; i<z.size(); i++) {
    double ew = std::abs(z.coeff(i)-z2.coeff(i));
    if (ew > err_z) err_z = ew;
  }
  return err_z;
}

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a Gram-Schmidt-based dimension reduction example.\n");
 
    int d = 4;
    CLP.setOption("d", &d, "Stochastic dimension");

    int d2 = 1;
    CLP.setOption("d2", &d2, "Intermediate stochastic dimension");

    int p = 5;
    CLP.setOption("p", &p, "Polynomial order");

    int p2 = 5;
    CLP.setOption("p2", &p2, "Intermediate polynomial order");

    double pole = 10.0;
    CLP.setOption("pole", &pole, "Pole location");

    double shift = 0.0;
    CLP.setOption("shift", &shift, "Shift location");

    double rank_threshold = 1.0e-10;
    CLP.setOption("rank_threshold", &rank_threshold, "Rank threshold");

    double reduction_tolerance = 1.0e-12;
    CLP.setOption("reduction_tolerance", &reduction_tolerance, "Quadrature reduction tolerance");

    bool verbose = false;
    CLP.setOption("verbose", "quiet", &verbose, "Verbose output");

    Basis_Reduction_Method basis_reduction_method = MONOMIAL_GS;
    CLP.setOption("basis_reduction_method", &basis_reduction_method, 
		  num_basis_reduction_method, basis_reduction_method_values, 
		  basis_reduction_method_names, "Basis reduction method");

    bool project = true;
    CLP.setOption("project", "no-project", &project, "Use Projected Lanczos Method");

    bool use_stieltjes = false;
    CLP.setOption("stieltjes", "no-stieltjes", &use_stieltjes, "Use Old Stieltjes Method");

    CLP.parse( argc, argv );

    std::cout << "Summary of command line options:" << std::endl
	      << "\tbasis_reduction_method = " 
	      << basis_reduction_method_names[basis_reduction_method] 
	      << std::endl
	      << "\tproject                = " << project << std::endl
	      << "\tstieljtes              = " << use_stieltjes << std::endl
	      << "\tp                      = " << p << std::endl
	      << "\tp2                     = " << p2 << std::endl
	      << "\td                      = " << d << std::endl
	      << "\td2                     = " << d2 << std::endl
	      << "\tpole                   = " << pole << std::endl
	      << "\tshift                  = " << shift << std::endl
	      << "\trank_threshold         = " << rank_threshold << std::endl
	      << "\treduction_tolerance    = " << reduction_tolerance 
	      << std::endl
	      << "\tverbose                = " << verbose << std::endl
	      << std::endl << std::endl;

    // Create product basis
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] = Teuchos::rcp(new basis_type(p, true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    std::cout << "original basis size = " << basis->size() << std::endl;

    // Tensor product quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis,
    // 								 p+1));

    std::cout << "original quadrature size = " << quad->size() << std::endl;

    // Triple product tensor
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor(basis->size());
    
    // Quadrature expansion
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > quad_exp = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, 
								    Cijk, 
								    quad));

    // Create approximation
    Teuchos::Array<pce_type> x(d);
    for (int i=0; i<d; i++) {
      x[i].copyForWrite();
      x[i].reset(quad_exp);
      x[i].term(i,1) = 1.0;
    }
    Teuchos::Array<pce_type> x2(d2);
    for (int i=0; i<d2; i++) {
      x2[i] = x[i];
    }
    
    // Compute PCE via quadrature expansion
    pce_type y = f(x, pole, shift);
    pce_type z = g(x2, y);
    
    // Create new basis from (x2,y)
    Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > pces(d2+1);
    for (int i=0; i<d2; i++)
      pces[i] = x2[i].getOrthogPolyApprox();
    pces[d2] = y.getOrthogPolyApprox();
    Teuchos::ParameterList params;
    if (basis_reduction_method == MONOMIAL_GS)
      params.set("Reduced Basis Method", "Monomial Proj Gram-Schmidt");
    else if (basis_reduction_method == LANCZOS)
      params.set("Reduced Basis Method", "Product Lanczos");
    else if (basis_reduction_method == LANCZOS_GS)
      params.set("Reduced Basis Method", "Product Lanczos Gram-Schmidt");
    params.set("Verbose", verbose);
    params.set("Project", project);
    //params.set("Normalize", false);
    params.set("Use Old Stieltjes Method", use_stieltjes);
    params.set("Basis Reduction Method", "Column-pivoted QR");
    //params.set("Basis Reduction Method", "SVD");
    params.set("Orthogonalization Method", "Householder");
    //params.set("Orthogonalization Method", "Classical Gram-Schmidt");
    //params.set("Orthogonalization Method", "Modified Gram-Schmidt");
    params.set("Rank Threshold", rank_threshold);
    Teuchos::ParameterList& red_quad_params = 
      params.sublist("Reduced Quadrature");
    //red_quad_params.set("Reduced Quadrature Method", "Column-Pivoted QR");
    //red_quad_params.set("Reduced Quadrature Method", "Pseudo-Inverse");
    //red_quad_params.set("Reduced Quadrature Method", "None");
    red_quad_params.set("Reduced Quadrature Method", "L1 Minimization");
    red_quad_params.set("LP Solver", "qpOASES");
    red_quad_params.set("Reduction Tolerance", reduction_tolerance);
    red_quad_params.set("Verbose", verbose);
    Stokhos::ReducedBasisFactory<int,double> factory(params);
    Teuchos::RCP< Stokhos::ReducedPCEBasis<int,double> > gs_basis = 
      factory.createReducedBasis(p2, pces, quad, Cijk);
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad =
      gs_basis->getReducedQuadrature();

    std::cout << "reduced basis size = " << gs_basis->size() << std::endl;
    std::cout << "reduced quadrature size = " << gs_quad->size() << std::endl;
    
    // Triple product tensor & expansion
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk;
    Teuchos::RCP< Teuchos::ParameterList > gs_exp_params = 
      Teuchos::rcp(new Teuchos::ParameterList);
    if (basis_reduction_method == LANCZOS)
      gs_Cijk = gs_basis->computeTripleProductTensor(gs_basis->size());
    else {
      gs_Cijk = Teuchos::null;
      gs_exp_params->set("Use Quadrature for Times", true);
    }
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > gs_quad_exp = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		     gs_basis, gs_Cijk, gs_quad, gs_exp_params));

    // Create new expansions
    Teuchos::Array<pce_type> x2_gs(d2);
    for (int i=0; i<d2; i++) {
      x2_gs[i].copyForWrite();
      x2_gs[i].reset(gs_quad_exp);
      gs_basis->transformFromOriginalBasis(x2[i].coeff(), x2_gs[i].coeff());
      pce_type xx(quad_exp);
      gs_basis->transformToOriginalBasis(x2_gs[i].coeff(), xx.coeff());
      double x_err = coeff_error(x2[i], xx);
      std::cout << "x2[" << i << "] coeff error = " << x_err << std::endl;
    }
    pce_type y_gs(gs_quad_exp); 
    gs_basis->transformFromOriginalBasis(y.coeff(), y_gs.coeff());
    pce_type yy(quad_exp);
    gs_basis->transformToOriginalBasis(y_gs.coeff(), yy.coeff());
    double y_err = coeff_error(y, yy);
    std::cout << "y coeff error = " << y_err << std::endl;
    
    // Compute z_gs = g(x2_gs, y_gs) in Gram-Schmidt basis
    pce_type z_gs = g(x2_gs, y_gs);
    
    // Project z_gs back to original basis
    pce_type z2(quad_exp);
    gs_basis->transformToOriginalBasis(z_gs.coeff(), z2.coeff());

    if (verbose) {
      std::cout << "z = " << std::endl << z;
      std::cout << "z2 = " << std::endl << z2;
      std::cout << "z_gs = " << std::endl << z_gs;
    }
    double err_z = coeff_error(z, z2);
    
    std::cout.precision(12);
    std::cout.setf(std::ios::scientific);
    std::cout << "z.mean()       = " << z.mean() << std::endl
	      << "z2.mean()      = " << z2.mean() << std::endl
	      << "mean error     = " 
	      << std::abs(z.mean()-z2.mean())/std::abs(z.mean()) << std::endl
	      << "z.std_dev()    = " << z.standard_deviation() << std::endl
	      << "z2.std_dev()   = " << z2.standard_deviation() << std::endl
	      << "std_dev error  = " 
	      << std::abs(z.standard_deviation()-z2.standard_deviation())/std::abs(z.standard_deviation())
	      << std::endl
	      << "z coeff error  = " << err_z << std::endl;
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
