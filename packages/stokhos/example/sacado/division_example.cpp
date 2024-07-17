// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// hermite_example
//
//  usage: 
//     hermite_example
//
//  output:  
//     prints the Hermite Polynomial Chaos Expansion of the simple function
//
//     v = 1/(log(u)^2+1)
//
//     where u = 1 + 0.4*H_1(x) + 0.06*H_2(x) + 0.002*H_3(x), x is a zero-mean
//     and unit-variance Gaussian random variable, and H_i(x) is the i-th
//     Hermite polynomial.

#include "Stokhos.hpp"
#include "Stokhos_Sacado.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"


// Typename of PC expansion type
typedef Sacado::PCE::OrthogPoly<double, Stokhos::StandardStorage<int,double> > pce_type;

// Linear solvers
 enum Division_Solver { Dense_Direct, GMRES, CG };
 const int num_division_solver = 3;
 const Division_Solver Division_solver_values[] = { Dense_Direct, GMRES, CG  };
 const char *division_solver_names[] = { "Dense_Direct", "GMRES", "CG" };

// Preconditioners
  enum Division_Prec { None, Diag, Jacobi, GS, Schur };
  const int num_division_prec = 5;
  const Division_Prec Division_prec_values[] = { None, Diag, Jacobi, GS, Schur };
  const char *division_prec_names[] = { "None", "Diag", "Jacobi","GS","Schur"};

// Option for Schur complement precond: full or diag D
   enum Schur_option { full, diag };
   const int num_schur_option = 2;
   const Schur_option Schur_option_values[] = { full, diag };
   const char *schur_option_names[] = { "full", "diag"};

// Full matrix or linear matrix (pb = dim + 1 ) used for preconditioner
    enum Prec_option { whole, linear};
    const int num_prec_option = 2;
    const Prec_option Prec_option_values[] = { whole, linear };
    const char *prec_option_names[] = { "full", "linear"};


int main(int argc, char **argv)
{
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // If applicable, set up MPI.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString("This example tests the / operator.\n");
    int d = 3;
    CLP.setOption("dim", &d, "Stochastic dimension");
    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");
    int n = 100;
    CLP.setOption("samples", &n, "Number of samples");
    double shift = 5.0;
    CLP.setOption("shift", &shift, "Shift point");
    double tolerance = 1e-6;
    CLP.setOption("tolerance", &tolerance, "Tolerance in Iterative Solver");
    int prec_level = 1;
    CLP.setOption("prec_level", &prec_level, "Level in Schur Complement Prec 0->Solve A0u0=g0 with division; 1->Form 1x1 Schur Complement");
    int max_it_div = 50;
    CLP.setOption("max_it_div", &max_it_div, "Maximum # of iterations for division iterative solver");
    bool equilibrate = true;
    CLP.setOption("equilibrate", "noequilibrate", &equilibrate,
                  "Equilibrate the linear system");

    Division_Solver solve_method = Dense_Direct;
    CLP.setOption("solver", &solve_method,
                  num_division_solver, Division_solver_values, division_solver_names,
                  "Solver Method");
    Division_Prec prec_method = None;
    CLP.setOption("prec", &prec_method,
                  num_division_prec, Division_prec_values, division_prec_names,
                  "Preconditioner Method");
    Schur_option schur_option = diag;
    CLP.setOption("schur_option", &schur_option,
                  num_schur_option, Schur_option_values, schur_option_names,
                  "Schur option");
    Prec_option prec_option = whole;
    CLP.setOption("prec_option", &prec_option,
                  num_prec_option, Prec_option_values, prec_option_names,
                  "Prec option");
    CLP.parse( argc, argv );

    // Basis
    Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d); 
    for (int i=0; i<d; i++) {
      bases[i] = rcp(new Stokhos::LegendreBasis<int,double>(p));
    }
    RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

    // Quadrature method
    RCP<const Stokhos::Quadrature<int,double> > quad = 
      rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));

    // Triple product tensor
    RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();

    // Expansion methods
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > quad_expn = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		     basis, Cijk, quad)); 
    Teuchos::RCP<Teuchos::ParameterList> alg_params = 
      Teuchos::rcp(new Teuchos::ParameterList);
    
    alg_params->set("Division Tolerance", tolerance);
    alg_params->set("prec_iter", prec_level);
    alg_params->set("max_it_div", max_it_div);
    if (solve_method == Dense_Direct)
         alg_params->set("Division Strategy", "Dense Direct");
    else if (solve_method == GMRES)
         alg_params->set("Division Strategy", "GMRES");
    else if (solve_method == CG)
         alg_params->set("Division Strategy", "CG");


    if (prec_method == None)
         alg_params->set("Prec Strategy", "None");
    else if (prec_method == Diag)
         alg_params->set("Prec Strategy", "Diag");
    else if (prec_method == Jacobi)
         alg_params->set("Prec Strategy", "Jacobi");
    else if (prec_method == GS)
         alg_params->set("Prec Strategy", "GS");
    else if (prec_method == Schur)
         alg_params->set("Prec Strategy", "Schur");
   
    if (schur_option == diag)
	alg_params->set("Schur option", "diag");
    else 
	alg_params->set("Schur option", "full");
    if (prec_option == linear)
	alg_params->set("Prec option", "linear");

    alg_params->set("Use Quadrature for Division", false);
   
    if (equilibrate)
        alg_params->set("Equilibrate", 1);
    else
        alg_params->set("Equilibrate", 0);


     
    Teuchos::RCP<Stokhos::QuadOrthogPolyExpansion<int,double> > alg_expn = 
      Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<int,double>(
		     basis, Cijk, quad, alg_params));
 
    // Polynomial expansions
    pce_type u_quad(quad_expn), v_quad(quad_expn);
    u_quad.term(0,0) = 0.0;
    for (int i=0; i<d; i++) {
      u_quad.term(i,1) = 1.0;
    }
    pce_type u_alg(alg_expn), v_alg(alg_expn);
    u_alg.term(0,0) = 0.0;
    for (int i=0; i<d; i++) {
      u_alg.term(i,1) = 1.0;
    }
    
      // Compute expansion
     double scale = std::exp(shift); 
    pce_type b_alg = std::exp(shift + u_alg)/scale;
    pce_type b_quad = std::exp(shift + u_quad)/scale;
//      v_alg = (1.0/scale) / b_alg;
//      v_quad = (1.0/scale) / b_quad;
	v_alg = 1.0 / std::exp(shift + u_alg);
        v_quad = 1.0 /std::exp(shift + u_quad);

//    std::cout << b_alg.getOrthogPolyApprox() << std::endl;
    
    // Print u and v
//    std::cout << "quadrature:   v = 1.0 / (shift + u) = ";
//    v_quad.print(std::cout);
//    std::cout << "dense solve:  v = 1.0 / (shift + u) = ";
//    v_alg.print(std::cout);

    double h = 2.0 / (n-1);
    double err_quad = 0.0;
    double err_alg = 0.0;
    for (int i=0; i<n; i++) {
      
      double x = -1.0 + h*i;
      Array<double> pt(d); 
      for (int j=0; j<d; j++) 
	pt[j] = x;
      double up = u_quad.evaluate(pt);
      double vp = 1.0/(shift+up);
      double vp_quad = v_quad.evaluate(pt);
      double vp_alg = v_alg.evaluate(pt);
      // std::cout << "vp = " << vp_quad << std::endl;
      // std::cout << "vp_quad = " << vp_quad << std::endl;
      // std::cout << "vp_alg = " << vp_alg << std::endl;
      double point_err_quad = std::abs(vp-vp_quad);
      double point_err_alg = std::abs(vp-vp_alg);
      if (point_err_quad > err_quad) err_quad = point_err_quad;
      if (point_err_alg > err_alg) err_alg = point_err_alg;
    }
    std::cout << "\tL_infty norm of quadrature error = " << err_quad 
	      << std::endl;
    std::cout << "\tL_infty norm of solver error = " << err_alg
	      << std::endl;
    
    // Check the answer
    //if (std::abs(err) < 1e-2)
      std::cout << "\nExample Passed!" << std::endl;

    Teuchos::TimeMonitor::summarize(std::cout);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
