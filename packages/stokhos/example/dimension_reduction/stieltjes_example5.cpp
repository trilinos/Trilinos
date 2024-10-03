// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <iomanip>

#include "Stokhos.hpp"
#include "Stokhos_Sacado.hpp"

#include "Stokhos_ReducedBasisFactory.hpp"
#include "Stokhos_SDMUtils.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

template <typename ordinal_type, typename scalar_type>
void
f_func(const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& A,
       const Teuchos::SerialDenseVector<ordinal_type,scalar_type>& x,
       double shift,
       Teuchos::SerialDenseVector<ordinal_type,scalar_type>& f)
{
  ordinal_type n = A.numCols();
  Teuchos::SerialDenseVector<ordinal_type,scalar_type> y(n), s(n);
  s.putScalar(shift);
  y.assign(x);
  y -= s;
  f.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, A, y, 0.0);
}

template <typename ordinal_type, typename scalar_type>
scalar_type
g_func(const Teuchos::SerialDenseVector<ordinal_type,scalar_type>& f)
{
  return std::exp(-f.dot(f));
  //return -f.dot(f);
}

typedef Stokhos::LegendreBasis<int,double> basis_type;
typedef Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_type;
typedef Teuchos::SerialDenseMatrix<int,pce_type> SDM;
typedef Teuchos::SerialDenseVector<int,pce_type> SDV;
typedef Teuchos::SerialDenseMatrix<int,double> SDM0;
typedef Teuchos::SerialDenseVector<int,double> SDV0;

// measure transformation approaches
enum MT_METHOD { MT_STIELTJES, MT_LANCZOS, MT_GRAM_SCHMIDT };
const int num_mt_method = 3;
const MT_METHOD mt_method_values[] = { 
  MT_STIELTJES, MT_LANCZOS, MT_GRAM_SCHMIDT };
const char *mt_method_names[] = { "Stieltjes", "Lanczos", "Gram-Schmidt" };

double rel_err(double a, double b) {
  return std::abs(a-b)/std::abs(b);
}

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a Stieltjes-based dimension reduction example.\n");

    int n = 4;
    CLP.setOption("n", &n, "Number of random variables");

    int m = 2;
    CLP.setOption("m", &m, "Number of intermediate variables");
 
    int pmin = 1;
    CLP.setOption("pmin", &pmin, "Starting expansion order");

    int pmax = 7;
    CLP.setOption("pmax", &pmax, "Final expansion order");

    MT_METHOD mt_method = MT_LANCZOS;
    CLP.setOption("mt_method", &mt_method, 
		  num_mt_method, mt_method_values, mt_method_names, 
		  "Measure transformation method");

    bool normalize = false;
    CLP.setOption("normalize", "unnormalize", &normalize, 
		  "Normalize PC basis");

    bool project_integrals = false;
    CLP.setOption("project", "no_project", &project_integrals, 
		  "Project integrals");

#ifdef HAVE_STOKHOS_DAKOTA
    bool sparse_grid = true;
    CLP.setOption("sparse", "tensor", &sparse_grid, 
		  "Use Sparse or Tensor grid");
#else
    bool sparse_grid = false;
#endif

    double rank_threshold = 1.0e-12;
    CLP.setOption("rank_threshold", &rank_threshold, "Rank threshold");

    double reduction_tolerance = 1.0e-12;
    CLP.setOption("reduction_tolerance", &reduction_tolerance, "Quadrature reduction tolerance");

    bool verbose = false;
    CLP.setOption("verbose", "quient", &verbose, 
		  "Verbose output");
    
    CLP.parse( argc, argv );

    std::cout << "Summary of command line options:" << std::endl
	      << "\tm                   = " << m << std::endl
	      << "\tn                   = " << n << std::endl
	      << "\tmt_method           = " << mt_method_names[mt_method] << std::endl
	      << "\tnormalize           = " << normalize << std::endl
	      << "\tproject             = " << project_integrals << std::endl
	      << "\tsparse              = " << sparse_grid << std::endl
	      << "\trank_threshold      = " << rank_threshold << std::endl
	      << "\treduction_tolerance = " << reduction_tolerance << std::endl;

    int np = pmax-pmin+1;
    
    Teuchos::Array<double> mean(np), mean_st(np), std_dev(np), std_dev_st(np);
    Teuchos::Array<double> pt(np), pt_st(np);

    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(n); 

    Teuchos::Array<double> eval_pt(n, 0.56789);
    double pt_true = 0.0;

    SDM0 B0(n,n), Q0, R0;
    B0.random();
    Teuchos::Array<double> w(n, 1.0);
    Stokhos::QR_MGS(n, B0, w, Q0, R0);
    SDM0 A0(Teuchos::View, B0, m, n);
    
    // Loop over orders
    int idx = 0;
    for (int p=pmin; p<=pmax; p++) {

      std::cout << "p = " << p;
      
      // Create product basis
      for (int i=0; i<n; i++)
	bases[i] = Teuchos::rcp(new basis_type(p, normalize));
      Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
      std::cout << ", basis sz = " << basis->size();
      
      // Quadrature
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
#ifdef HAVE_STOKHOS_DAKOTA
      if (sparse_grid)
      	quad = 
	  Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis, p));
#endif
      if (!sparse_grid)
	quad = 
	  Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
      std::cout << ", quad sz = " << quad->size();

      // Triple product tensor
      Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
	basis->computeTripleProductTensor();
      
      // Quadrature expansion
      Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > quad_exp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, Cijk, quad));

      // Create approximation
      SDV x(n);
      for (int i=0; i<n; i++) {
	x[i].copyForWrite();
	x[i].reset(quad_exp);
	if (normalize)
	  x[i].term(i, 1) = 1.0/std::sqrt(3);
	else
	  x[i].term(i, 1) = 1.0;
      }

      // Create matrix
      SDM A(m,n);
      for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	  A(i,j) = A0(i,j);

      // Compute f = A*(x-1)
      SDV f(m);
      f_func(A, x, 1.0, f);

      //std::cout << "f = " << f << std::endl;

      // Compute g = exp(-f^T*f)
      pce_type g = g_func(f);

      // compute true point
      SDV0 x0(n), f0(m);
      for (int i=0; i<n; i++)
	x0[i] = x[i].evaluate(eval_pt);
      f_func(A0, x0, 1.0, f0);
      pt_true = g_func(f0);
	
      // Compute reduced basis
      Teuchos::ParameterList params;
      params.set("Verbose", verbose);
      if (mt_method == MT_GRAM_SCHMIDT)
	params.set("Reduced Basis Method", "Monomial Proj Gram-Schmidt");
      else if (mt_method == MT_LANCZOS)
	params.set("Reduced Basis Method", "Product Lanczos");
      else if (mt_method == MT_STIELTJES) {
	params.set("Reduced Basis Method", "Product Lanczos");
	params.set("Use Old Stieltjes Method", true);
      }
      params.set("Project", project_integrals);
      params.set("Normalize", normalize);
      params.set("Rank Threshold", rank_threshold);
      Teuchos::ParameterList& red_quad_params = 
	params.sublist("Reduced Quadrature");
      red_quad_params.set("Reduction Tolerance", reduction_tolerance);
      red_quad_params.set("Verbose", verbose);
      Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > pces(m);
      for (int i=0; i<m; i++)
	pces[i] = f[i].getOrthogPolyApprox();
      Stokhos::ReducedBasisFactory<int,double> factory(params);
      Teuchos::RCP< Stokhos::ReducedPCEBasis<int,double> > gs_basis = 
	factory.createReducedBasis(p, pces, quad, Cijk);
      Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad =
	gs_basis->getReducedQuadrature();
      Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk =
	Teuchos::null;
      Teuchos::RCP< Teuchos::ParameterList > gs_exp_params = 
	Teuchos::rcp(new Teuchos::ParameterList);
      gs_exp_params->set("Use Quadrature for Times", true);
      Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > st_quad_exp =
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		       gs_basis, 
		       gs_Cijk,
		       gs_quad,
		       gs_exp_params));
      std::cout << ", red. basis sz = " <<gs_basis->size();
      std::cout << ", red. quad sz = " << gs_quad->size();

      SDV f_st(m);
      for (int i=0; i<m; i++) {
	f_st(i).copyForWrite();
	f_st(i).reset(st_quad_exp);
	gs_basis->transformFromOriginalBasis(f(i).coeff(), f_st(i).coeff());
      }

      // Compute g_st = exp(-f_st^T*f_st)
      pce_type g_st = g_func(f_st);
      
      // Project g_st back to original basis
      pce_type g2(quad_exp);
      gs_basis->transformToOriginalBasis(g_st.coeff(), g2.coeff());

      // std::cout.precision(12);
      // std::cout << g;
      // std::cout << g2;
      // std::cout << g_st;
      mean[idx] = g.mean();
      mean_st[idx] = g2.mean();
      std_dev[idx] = g.standard_deviation();
      std_dev_st[idx] = g2.standard_deviation();
      pt[idx] = g.evaluate(eval_pt);
      pt_st[idx] = g2.evaluate(eval_pt);
      idx++;

      std::cout << std::endl;
    }

    idx = 0;
    int wi=10;
    std::cout << "Statistical error:" << std::endl;
    std::cout << "p  " 
	      << std::setw(wi) << "mean" << "  " 
	      << std::setw(wi) << "mean_st" << "  "
	      << std::setw(wi) << "std_dev" << "  "
	      << std::setw(wi) << "std_dev_st" << "  "
	      << std::setw(wi) << "point" << "  "
	      << std::setw(wi) << "point_st" << std::endl;
    for (int p=pmin; p<pmax; p++) {
      std::cout.precision(3);
      std::cout.setf(std::ios::scientific);
      std::cout << p << "  " 
		<< std::setw(wi) << rel_err(mean[idx], mean[np-1]) << "  "
		<< std::setw(wi) << rel_err(mean_st[idx], mean[np-1]) << "  "
		<< std::setw(wi) << rel_err(std_dev[idx], std_dev[np-1]) << "  "
		<< std::setw(wi) << rel_err(std_dev_st[idx], std_dev[np-1]) 
		<< "  "
		<< std::setw(wi) << rel_err(pt[idx], pt_true) << "  "
		<< std::setw(wi) << rel_err(pt_st[idx], pt_true) 
		<< std::endl;
      idx++;
    }
      
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
