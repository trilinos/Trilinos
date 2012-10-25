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
#include "Teuchos_TabularOutputter.hpp"

// quadrature methods
enum Quadrature_Method { TENSOR, SPARSE };
static const int num_quadrature_method = 2;
static const Quadrature_Method quadrature_method_values[] = { 
  TENSOR, SPARSE };
static const char *quadrature_method_names[] = { 
  "Tensor", "Sparse" };

// reduced basis methods
enum Reduced_Basis_Method { LANCZOS, MONOMIAL_PROJ_GS, MONOMIAL_PROJ_GS2, MONOMIAL_GS, LANCZOS_GS };
static const int num_reduced_basis_method = 5;
static const Reduced_Basis_Method reduced_basis_method_values[] = { 
  LANCZOS, MONOMIAL_PROJ_GS, MONOMIAL_PROJ_GS2, MONOMIAL_GS, LANCZOS_GS };
static const char *reduced_basis_method_names[] = { 
  "Lanczos", "Monomial-Proj-GS", "Monomial-Proj-GS2", "Monomial-GS", "Lanczos-GS" };

// orthogonalization methods
enum Orthogonalization_Method { HOUSEHOLDER, HOUSEHOLDER_NP, CGS, MGS, MGSRO, MGSNP, MGSNPRO, SVD };
static const int num_orthogonalization_method = 8;
static const Orthogonalization_Method orthogonalization_method_values[] = { 
  HOUSEHOLDER, HOUSEHOLDER_NP, CGS, MGS, MGSRO, MGSNP, MGSNPRO, SVD };
static const char *orthogonalization_method_names[] = { 
  "Householder", "Householder without Pivoting", "Classical Gram-Schmidt", "Modified Gram-Schmidt", "Modified Gram-Schmidt with Reorthogonalization", "Modified Gram-Schmidt without Pivoting", "Modified Gram-Schmidt without Pivoting with Reorthogonalization", "SVD" };

// quadrature reduction methods
enum Quadrature_Reduction_Method { NONE, QSQUARED, QSQUARED2, Q2 };
static const int num_quad_reduction_method = 4;
static const Quadrature_Reduction_Method quad_reduction_method_values[] = { 
  NONE, QSQUARED, QSQUARED2, Q2 };
static const char *quad_reduction_method_names[] = { 
  "None", "Q Squared", "Q Squared2", "Q2" };

// solver methods
enum Solver_Method { TRSM, GLPK, CLP, CLP_IP, QPOASES, BASIS_PURSUIT, ORTHOGONAL_MATCHING_PURSUIT };
static const int num_solver_method = 7;
static const Solver_Method solver_method_values[] = { 
  TRSM, GLPK, CLP, CLP_IP, QPOASES, BASIS_PURSUIT, ORTHOGONAL_MATCHING_PURSUIT };
static const char *solver_method_names[] = { 
  "TRSM", "GLPK", "Clp", "Clp-IP", "qpOASES", "Basis Pursuit", "Orthogonal Matching Pursuit" };

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
  int n = std::min(z.size(), z2.size());
  for (int i=0; i<n; i++) {
    double ew = std::abs(z.coeff(i)-z2.coeff(i));
    if (ew > err_z) err_z = ew;
  }
  return err_z;
}

double disc_orthog_error(const Stokhos::OrthogPolyBasis<int,double>& basis,
			 const Stokhos::Quadrature<int,double>& quad) {
  const Teuchos::Array<double>& weights = quad.getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& vals = 
    quad.getBasisAtQuadPoints();
  int nqp = quad.size();
  int npc = basis.size();
  double max_err = 0.0;

  // Loop over all basis function combinations
  for (int i=0; i<npc; ++i) {
    for (int j=0; j<npc; ++j) {

      // Compute inner product of ith and jth basis function
      double err = 0.0;
      for (int k=0; k<nqp; ++k)
	err += weights[k]*vals[k][i]*vals[k][j];

      // Subtract off what it should be
      if (i == j)
	err -= basis.norm_squared(i);

      // Accumulate max error
      if (std::abs(err) > max_err)
	max_err = std::abs(err);
    }
  }
  return max_err;
}

struct MyOptions {
  int d;
  int d2;
  int p_begin;
  int p_end;
  int p_truth;
  double pole;
  double shift;
  double rank_threshold;
  double rank_threshold2;
  double reduction_tolerance;
  bool verbose;
  Quadrature_Method quad_method;
  int level;
  Reduced_Basis_Method reduced_basis_method;
  Orthogonalization_Method orthogonalization_method;
  Quadrature_Reduction_Method quad_reduction_method;
  Solver_Method solver_method;
  Orthogonalization_Method quad_orthogonalization_method;
  bool eliminate_dependent_rows;
  bool use_Q;
  bool restrict_r;
  double objective_value;
  bool project;
  bool use_stieltjes;
};

class MyResults {
public:
  int basis_size;
  int quad_size;
  int reduced_basis_size;
  int reduced_quad_size;
  pce_type z;
  pce_type z_red;
  double mean_error;
  double std_dev_error;
  double coeff_error;
  double reduced_mean_error;
  double reduced_std_dev_error;
  double reduced_coeff_error;
  double disc_orthog_error;
};

void compute_pces(bool compute_z_red, int p, const MyOptions& options, 
		  MyResults& results)
{
  // Create product basis
  Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(options.d);
  for (int i=0; i<options.d; i++)
    bases[i] = Teuchos::rcp(new basis_type(p, true));
  Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

  results.basis_size = basis->size();
  
  // Quadrature
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad;
  if (options.quad_method == TENSOR)
    quad = 
      Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
  else if (options.quad_method == SPARSE) {
#ifdef HAVE_STOKHOS_DAKOTA
    if (options.level == -1)
      quad = 
	Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis));
    else
      quad = 
	Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(
		       basis, p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Sparse grid quadrature only supported when compiled with Dakota!");
#endif
    }

  results.quad_size = quad->size();

  // Triple product tensor
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
    basis->computeTripleProductTensor(basis->size());
    
  // Quadrature expansion
  Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > quad_exp = 
    Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(basis, 
								  Cijk, 
								  quad));

  // Create approximation
  Teuchos::Array<double> point(options.d, 1.0);
  Teuchos::Array<double> basis_vals(basis->size());
  basis->evaluateBases(point, basis_vals);
  Teuchos::Array<pce_type> x(options.d);
  for (int i=0; i<options.d; i++) {
    x[i].copyForWrite();
    x[i].reset(quad_exp);
    x[i].term(i,1) = 1.0 / basis_vals[i+1];
  }
  Teuchos::Array<pce_type> x2(options.d2);
  for (int i=0; i<options.d2; i++) {
    x2[i] = x[i];
  }
  
  // Compute PCE via quadrature expansion
  pce_type y = f(x, options.pole, options.shift);
  results.z = g(x2, y);

  if (!compute_z_red)
    return;
    
  // Create new basis from (x2,y)
  Teuchos::Array< Stokhos::OrthogPolyApprox<int,double> > pces(options.d2+1);
  for (int i=0; i<options.d2; i++)
    pces[i] = x2[i].getOrthogPolyApprox();
  pces[options.d2] = y.getOrthogPolyApprox();
  Teuchos::ParameterList params;
  if (options.reduced_basis_method == MONOMIAL_PROJ_GS)
    params.set("Reduced Basis Method", "Monomial Proj Gram-Schmidt");
  else if (options.reduced_basis_method == MONOMIAL_PROJ_GS2)
    params.set("Reduced Basis Method", "Monomial Proj Gram-Schmidt2");
  else if (options.reduced_basis_method == MONOMIAL_GS)
    params.set("Reduced Basis Method", "Monomial Gram-Schmidt");
  else if (options.reduced_basis_method == LANCZOS)
    params.set("Reduced Basis Method", "Product Lanczos");
  else if (options.reduced_basis_method == LANCZOS_GS)
    params.set("Reduced Basis Method", "Product Lanczos Gram-Schmidt");
  params.set("Verbose", options.verbose);
  params.set("Project", options.project);
  //params.set("Normalize", false);
  params.set("Use Old Stieltjes Method", options.use_stieltjes);
  params.set("Orthogonalization Method", 
	       orthogonalization_method_names[options.orthogonalization_method]);
  params.set("Rank Threshold", options.rank_threshold);
  Teuchos::ParameterList& red_quad_params = 
    params.sublist("Reduced Quadrature");
  red_quad_params.set(
    "Reduced Quadrature Method", 
    quad_reduction_method_names[options.quad_reduction_method]);
  red_quad_params.set(
    "Solver Method", solver_method_names[options.solver_method]);
  red_quad_params.set(
    "Eliminate Dependent Rows", options.eliminate_dependent_rows);
  red_quad_params.set("Write MPS File", false);
  red_quad_params.set("Reduction Tolerance", options.reduction_tolerance);
  red_quad_params.set("Verbose", options.verbose);
  red_quad_params.set("Objective Value", options.objective_value);
  red_quad_params.set("Q2 Rank Threshold", options.rank_threshold2);
  red_quad_params.set(
    "Orthogonalization Method", 
    orthogonalization_method_names[options.quad_orthogonalization_method]);
  red_quad_params.set("Use Q in LP", options.use_Q);
  red_quad_params.set("Restrict Rank", options.restrict_r);
  red_quad_params.set("Order Restriction", 2*p);
  Stokhos::ReducedBasisFactory<int,double> factory(params);
  Teuchos::RCP< Stokhos::ReducedPCEBasis<int,double> > gs_basis = 
    factory.createReducedBasis(p, pces, quad, Cijk);
  Teuchos::RCP<const Stokhos::Quadrature<int,double> > gs_quad =
    gs_basis->getReducedQuadrature();

  results.reduced_basis_size = gs_basis->size();
  results.reduced_quad_size = gs_quad->size();
    
  // Triple product tensor & expansion
  Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > gs_Cijk;
  Teuchos::RCP< Teuchos::ParameterList > gs_exp_params = 
    Teuchos::rcp(new Teuchos::ParameterList);
  if (options.reduced_basis_method == LANCZOS)
    gs_Cijk = gs_basis->computeTripleProductTensor(gs_basis->size());
  else {
    gs_Cijk = Teuchos::null;
    gs_exp_params->set("Use Quadrature for Times", true);
  }
  Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<int,double> > gs_quad_exp = 
    Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		   gs_basis, gs_Cijk, gs_quad, gs_exp_params));

  // Create new expansions
  Teuchos::Array<pce_type> x2_gs(options.d2);
  for (int i=0; i<options.d2; i++) {
    x2_gs[i].copyForWrite();
    x2_gs[i].reset(gs_quad_exp);
    gs_basis->transformFromOriginalBasis(x2[i].coeff(), x2_gs[i].coeff());
  }
  pce_type y_gs(gs_quad_exp); 
  gs_basis->transformFromOriginalBasis(y.coeff(), y_gs.coeff());
  
  // Compute z_gs = g(x2_gs, y_gs) in Gram-Schmidt basis
  pce_type z_gs = g(x2_gs, y_gs);
    
  // Project z_gs back to original basis
  pce_type z2(quad_exp);
  gs_basis->transformToOriginalBasis(z_gs.coeff(), z2.coeff());
  results.z_red = z2;

  // Compute discrete orthogonality error
  results.disc_orthog_error = disc_orthog_error(*gs_basis, *gs_quad);
}

int main(int argc, char **argv)
{
  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a Gram-Schmidt-based dimension reduction example.\n");
    MyOptions options;

    options.d = 4;
    CLP.setOption("d", &options.d, "Stochastic dimension");

    options.d2 = 1;
    CLP.setOption("d2", &options.d2, "Intermediate stochastic dimension");

    options.p_begin = 1;
    CLP.setOption("p_begin", &options.p_begin, "Starting polynomial order");

    options.p_end = 5;
    CLP.setOption("p_end", &options.p_end, "Ending polynomial order");

    options.p_truth = 7;
    CLP.setOption("p_truth", &options.p_truth, "Truth polynomial order");

    options.pole = 10.0;
    CLP.setOption("pole", &options.pole, "Pole location");

    options.shift = 0.0;
    CLP.setOption("shift", &options.shift, "Shift location");

    options.rank_threshold = 1.0e-120;
    CLP.setOption("rank_threshold", &options.rank_threshold, "Rank threshold");

    options.rank_threshold2 = 1.0e-120;
    CLP.setOption("rank_threshold2", &options.rank_threshold2, "Rank threshold for Q2");

    options.reduction_tolerance = 1.0e-12;
    CLP.setOption("reduction_tolerance", &options.reduction_tolerance, "Quadrature reduction tolerance");

    options.verbose = false;
    CLP.setOption("verbose", "quiet", &options.verbose, "Verbose output");

    options.quad_method = TENSOR;
    CLP.setOption("quadrature_method", &options.quad_method, 
		  num_quadrature_method, quadrature_method_values, 
		  quadrature_method_names, "Quadrature method");

    options.level = -1;
    CLP.setOption("level", &options.level, 
		  "Sparse grid level (set to -1 to use default)");

    options.reduced_basis_method = MONOMIAL_GS;
    CLP.setOption("reduced_basis_method", &options.reduced_basis_method, 
		  num_reduced_basis_method, reduced_basis_method_values, 
		  reduced_basis_method_names, "Reduced basis method");

    options.orthogonalization_method = MGSRO;
    CLP.setOption("orthogonalization_method", 
		  &options.orthogonalization_method, 
		  num_orthogonalization_method, 
		  orthogonalization_method_values, 
		  orthogonalization_method_names, 
		  "Orthogonalization method");

    options.quad_reduction_method = QSQUARED;
    CLP.setOption("reduced_quadrature_method", &options.quad_reduction_method, 
		  num_quad_reduction_method, quad_reduction_method_values, 
		  quad_reduction_method_names, "Reduced quadrature method");

    options.solver_method = TRSM;
    CLP.setOption("solver_method", &options.solver_method, num_solver_method, 
		  solver_method_values,  solver_method_names, 
		  "Reduced quadrature solver method");

    options.quad_orthogonalization_method = MGSRO;
    CLP.setOption("quad_orthogonalization_method", 
		  &options.quad_orthogonalization_method, 
		  num_orthogonalization_method, 
		  orthogonalization_method_values, 
		  orthogonalization_method_names, 
		  "Quadrature Orthogonalization method");

    options.eliminate_dependent_rows = true;
    CLP.setOption("cpqr", "no-cpqr", &options.eliminate_dependent_rows, 
		  "Eliminate dependent rows in quadrature constraint matrix");

    options.use_Q = true;
    CLP.setOption("use-Q", "no-use-Q", &options.use_Q, "Use Q in LP");

    options.restrict_r = false;
    CLP.setOption("restrict-rank", "no-restrict-rank", &options.restrict_r, 
		  "Restrict rank in LP");

    options.objective_value = 0.0;
    CLP.setOption("objective_value", &options.objective_value, 
		  "Value for LP objective function");

    options.project = true;
    CLP.setOption("project", "no-project", &options.project, 
		  "Use Projected Lanczos Method");

    options.use_stieltjes = false;
    CLP.setOption("stieltjes", "no-stieltjes", &options.use_stieltjes, 
		  "Use Old Stieltjes Method");

    CLP.parse( argc, argv );

    std::cout << "Summary of command line options:" << std::endl
	      << "\tquadrature_method           = " 
	      << quadrature_method_names[options.quad_method] 
	      << std::endl
	      << "\tlevel                       = " << options.level 
	      << std::endl
	      << "\treduced_basis_method        = " 
	      << reduced_basis_method_names[options.reduced_basis_method] 
	      << std::endl
	      << "\torthogonalization_method    = " 
	      << orthogonalization_method_names[options.orthogonalization_method] 
	      << std::endl
	      << "\tquadrature_reduction_method = " 
	      << quad_reduction_method_names[options.quad_reduction_method] 
	      << std::endl
	      << "\tsolver_method               = " 
	      << solver_method_names[options.solver_method] << std::endl
	      << "\tquad_orthogonalization_method = " 
	      << orthogonalization_method_names[options.quad_orthogonalization_method] 
	      << std::endl
	      << "\tcpqr                        = " << options.eliminate_dependent_rows 
	      << std::endl
	      << "\tuse-Q                       = " << options.use_Q << std::endl
	      << "\trestrict-rank               = " << options.restrict_r << std::endl
	      << "\tobjective_value             = " << options.objective_value << std::endl
	      << "\tproject                     = " << options.project << std::endl
	      << "\tstieljtes                   = " << options.use_stieltjes << std::endl
	      << "\tp_begin                     = " << options.p_begin << std::endl
	      << "\tp_end                       = " << options.p_end << std::endl
	      << "\tp_truth                     = " << options.p_truth << std::endl
	      << "\td                           = " << options.d << std::endl
	      << "\td2                          = " << options.d2 << std::endl
	      << "\tpole                        = " << options.pole << std::endl
	      << "\tshift                       = " << options.shift << std::endl
	      << "\trank_threshold              = " << options.rank_threshold << std::endl
	      << "\trank_threshold2             = " << options.rank_threshold2 << std::endl
	      << "\treduction_tolerance         = " << options.reduction_tolerance 
	      << std::endl
	      << "\tverbose                     = " << options.verbose << std::endl
	      << std::endl << std::endl;

    std::stringstream ss;
    typedef Teuchos::TabularOutputter TO;
    TO out(ss);
    out.setFieldTypePrecision(TO::DOUBLE, 5);
    out.pushFieldSpec("\"Order\"", TO::INT);
    out.pushFieldSpec("\"Basis Size\"", TO::INT);
    out.pushFieldSpec("\"Quad. Size\"", TO::INT);
    out.pushFieldSpec("\"Red. Basis Size\"", TO::INT);
    out.pushFieldSpec("\"Red. Quad. Size\"", TO::INT);
    out.pushFieldSpec("\"Coeff. Error\"", TO::DOUBLE);
    out.pushFieldSpec("\"Red. Coeff. Error\"", TO::DOUBLE);
    out.pushFieldSpec("\"Disc. Orthog. Error\"", TO::DOUBLE);
    out.outputHeader();

    MyResults results_truth;
    compute_pces(false, options.p_truth, options, results_truth);

    int n = options.p_end - options.p_begin + 1;
    Teuchos::Array<MyResults> results(n);
    for (int i=0; i<n; ++i) {
      int p = options.p_begin + i;
      compute_pces(true, p, options, results[i]);
      results[i].mean_error = 
	std::abs(results_truth.z.mean()-results[i].z.mean()) / 
	std::abs(results_truth.z.mean());
      results[i].std_dev_error = 
	std::abs(results_truth.z.standard_deviation()-results[i].z.standard_deviation()) / std::abs(results_truth.z.standard_deviation());
      results[i].coeff_error = coeff_error(results_truth.z, 
					   results[i].z);

      results[i].reduced_mean_error = 
	std::abs(results_truth.z.mean()-results[i].z_red.mean()) / 
	std::abs(results_truth.z.mean());
      results[i].reduced_std_dev_error = 
	std::abs(results_truth.z.standard_deviation()-results[i].z_red.standard_deviation()) / std::abs(results_truth.z.standard_deviation());
      results[i].reduced_coeff_error = coeff_error(results_truth.z, 
						   results[i].z_red);

      out.outputField(p);
      out.outputField(results[i].basis_size);
      out.outputField(results[i].quad_size);
      out.outputField(results[i].reduced_basis_size);
      out.outputField(results[i].reduced_quad_size);
      out.outputField(results[i].coeff_error);
      out.outputField(results[i].reduced_coeff_error);
      out.outputField(results[i].disc_orthog_error);
      out.nextRow();
    }
    std::cout << std::endl << ss.str() << std::endl;
    
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
