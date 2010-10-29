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

// ModelEvaluator implementing our problem
#include "twoD_diffusion_ME.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// NOX
#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_LinearSystem_Stratimikos.H"
#include "NOX_Epetra_LinearSystem_SGGS.hpp"
#include "NOX_Epetra_LinearSystem_SGJacobi.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Utilities
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"

// I/O utilities
#include "EpetraExt_VectorOut.h"

//The probability distribution of the random variables.
double uniform_weight(const double& x){
  return 1;
}

// Linear solvers
enum Krylov_Solver { AZTECOO, BELOS };
const int num_krylov_solver = 2;
const Krylov_Solver krylov_solver_values[] = { AZTECOO, BELOS };
const char *krylov_solver_names[] = { "AztecOO", "Belos" };

// SG solver approaches
enum SG_Solver { SG_KRYLOV, SG_GS, SG_JACOBI };
const int num_sg_solver = 3;
const SG_Solver sg_solver_values[] = { SG_KRYLOV, SG_GS, SG_JACOBI };
const char *sg_solver_names[] = { "Krylov", "Gauss-Seidel", "Jacobi" };

// Krylov methods
// enum Krylov_Method { GMRES, CG, FGMRES, RGMRES };
// const int num_krylov_method = 4;
// const Krylov_Method krylov_method_values[] = { GMRES, CG, FGMRES, RGMRES };
// const char *krylov_method_names[] = { "GMRES", "CG", "FGMRES", "RGMRES" };
enum Krylov_Method { GMRES, CG, FGMRES };
const int num_krylov_method = 3;
const Krylov_Method krylov_method_values[] = { GMRES, CG, FGMRES };
const char *krylov_method_names[] = { "GMRES", "CG", "FGMRES" };

// Krylov preconditioning approaches
enum SG_Prec { MEAN, GS, AGS, AJ, KP };
const int num_sg_prec = 5;
const SG_Prec sg_prec_values[] = { MEAN, GS, AGS, AJ, KP };
const char *sg_prec_names[] = { "Mean-Based", 
				"Gauss-Seidel", 
				"Approx-Gauss-Seidel", 
				"Approx-Jacobi", 
				"Kronecker-Product" };

// Random field types
enum SG_RF { UNIFORM, LOGNORMAL };
const int num_sg_rf = 2;
const SG_RF sg_rf_values[] = { UNIFORM, LOGNORMAL };
const char *sg_rf_names[] = { "Uniform", "Log-Normal" };

int main(int argc, char *argv[]) {   

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
  Teuchos::RCP<Epetra_Comm> Comm;
#ifdef HAVE_MPI
  Comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  Comm = Teuchos::rcp(new Epetra_SerialComm);
#endif

  int MyPID = Comm->MyPID();

  try {

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs a variety of stochastic Galerkin solvers.\n");

    int n = 32;
    CLP.setOption("num_mesh", &n, "Number of mesh points in each direction");

    bool symmetric = false;
    CLP.setOption("symmetric", "unsymmetric", &symmetric, 
		  "Symmetric discretization");

    SG_RF randField = UNIFORM;
    CLP.setOption("rand_field", &randField, 
		   num_sg_rf, sg_rf_values, sg_rf_names,
		  "Random field type");

    double mean = 0.2;
    CLP.setOption("mean", &mean, "Mean");

    double sigma = 0.1;
    CLP.setOption("std_dev", &sigma, "Standard deviation");

    double weightCut = 1.0;
    CLP.setOption("weight_cut", &weightCut, "Weight cut");

    int num_KL = 2;
    CLP.setOption("num_kl", &num_KL, "Number of KL terms");

    int p = 5;
    CLP.setOption("order", &p, "Polynomial order");

    bool normalize_basis = true;
    CLP.setOption("normalize", "unnormalize", &normalize_basis, 
		  "Normalize PC basis");
    
    SG_Solver solve_method = SG_KRYLOV;
    CLP.setOption("sg_solver", &solve_method, 
		  num_sg_solver, sg_solver_values, sg_solver_names, 
		  "SG solver method");

    Krylov_Method outer_krylov_method = GMRES;
    CLP.setOption("outer_krylov_method", &outer_krylov_method, 
		  num_krylov_method, krylov_method_values, krylov_method_names, 
		  "Outer Krylov method (for Krylov-based SG solver)");

    Krylov_Solver outer_krylov_solver = AZTECOO;
    CLP.setOption("outer_krylov_solver", &outer_krylov_solver, 
		  num_krylov_solver, krylov_solver_values, krylov_solver_names, 
		  "Outer linear solver");

    double outer_tol = 1e-12;
    CLP.setOption("outer_tol", &outer_tol, "Outer solver tolerance");

    int outer_its = 1000;
    CLP.setOption("outer_its", &outer_its, "Maximum outer iterations");

    Krylov_Method inner_krylov_method = GMRES;
    CLP.setOption("inner_krylov_method", &inner_krylov_method, 
		  num_krylov_method, krylov_method_values, krylov_method_names, 
		  "Inner Krylov method (for G-S, Jacobi, etc...)");

    Krylov_Solver inner_krylov_solver = AZTECOO;
    CLP.setOption("inner_krylov_solver", &inner_krylov_solver, 
		  num_krylov_solver, krylov_solver_values, krylov_solver_names, 
		  "Inner linear solver");

    double inner_tol = 3e-13;
    CLP.setOption("inner_tol", &inner_tol, "Inner solver tolerance");

    int inner_its = 1000;
    CLP.setOption("inner_its", &inner_its, "Maximum inner iterations");

    SG_Prec precMethod = AGS;
    CLP.setOption("sg_prec_method", &precMethod, 
		  num_sg_prec, sg_prec_values, sg_prec_names,
		  "Preconditioner method");

    double gs_prec_tol = 1e-1;
    CLP.setOption("gs_prec_tol", &gs_prec_tol, "Gauss-Seidel preconditioner tolerance");

    int gs_prec_its = 1;
    CLP.setOption("gs_prec_its", &gs_prec_its, "Maximum Gauss-Seidel preconditioner iterations");

    CLP.parse( argc, argv );

    if (MyPID == 0) {
      std::cout << "Summary of command line options:" << std::endl
		<< "\tnum_mesh            = " << n << std::endl
		<< "\tsymmetric           = " << symmetric << std::endl
		<< "\trand_field          = " << sg_rf_names[randField] 
		<< std::endl
		<< "\tmean                = " << mean << std::endl
		<< "\tstd_dev             = " << sigma << std::endl
		<< "\tweight_cut          = " << weightCut << std::endl
		<< "\tnum_kl              = " << num_KL << std::endl
		<< "\torder               = " << p << std::endl
		<< "\tnormalize_basis     = " << normalize_basis << std::endl
		<< "\tsg_solver           = " << sg_solver_names[solve_method] 
		<< std::endl
		<< "\touter_krylov_method = " 
		<< krylov_method_names[outer_krylov_method] << std::endl
		<< "\touter_krylov_solver = " 
		<< krylov_solver_names[outer_krylov_solver] << std::endl
		<< "\tinner_krylov_method = " 
		<< krylov_method_names[inner_krylov_method] << std::endl
		<< "\tinner_krylov_solver = " 
		<< krylov_solver_names[inner_krylov_solver] << std::endl
		<< "\tprec_method         = " << sg_prec_names[precMethod] 
		<< std::endl;
    }

    bool nonlinear_expansion;
    if (randField == UNIFORM)
      nonlinear_expansion = false;
    else if (randField == LOGNORMAL)
      nonlinear_expansion = true;
    bool scaleOP = true; 

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      if (randField == UNIFORM)
        bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p,normalize_basis));
      else if (randField == LOGNORMAL)      
        bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p,normalize_basis));

    //  bases[i] = Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("beta",p,&uniform_weight,-weightCut,weightCut,true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (nonlinear_expansion)
      Cijk = basis->computeTripleProductTensor(sz);
    else
      Cijk = basis->computeTripleProductTensor(num_KL+1);
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
									 Cijk));
    std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;
    
    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model =
      Teuchos::rcp(new twoD_diffusion_ME(Comm, n, num_KL, sigma, 
					 mean, basis, nonlinear_expansion,
					 symmetric));

    // Set up NOX parameters
    Teuchos::RCP<Teuchos::ParameterList> noxParams = 
      Teuchos::rcp(new Teuchos::ParameterList);

    // Set the nonlinear solver method
    noxParams->set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = noxParams->sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information", 
                    NOX::Utils::OuterIteration + 
                    NOX::Utils::OuterIterationStatusTest + 
                    NOX::Utils::InnerIteration +
		    //NOX::Utils::Parameters + 
		    NOX::Utils::Details + 
		    NOX::Utils::LinearSolverDetails +
                    NOX::Utils::Warning + 
                    NOX::Utils::Error);

    // Create printing utilities
    NOX::Utils utils(printParams);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = noxParams->sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = noxParams->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

    // Alternative linear solver list for Stratimikos
    Teuchos::ParameterList& stratLinSolParams = 
      newtonParams.sublist("Stratimikos Linear Solver");
    // Teuchos::ParameterList& noxStratParams = 
    //   stratLinSolParams.sublist("NOX Stratimikos Options");
    Teuchos::ParameterList& stratParams = 
      stratLinSolParams.sublist("Stratimikos");

    // Sublist for convergence tests
    Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
    statusParams.set("Test Type", "Combo");
    statusParams.set("Number of Tests", 2);
    statusParams.set("Combo Type", "OR");
    Teuchos::ParameterList& normF = statusParams.sublist("Test 0");
    normF.set("Test Type", "NormF");
    normF.set("Tolerance", outer_tol);
    normF.set("Scale Type", "Scaled");
    Teuchos::ParameterList& maxIters = statusParams.sublist("Test 1");
    maxIters.set("Test Type", "MaxIters");
    maxIters.set("Maximum Iterations", 1);

    // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> det_nox_interface = 
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));

     // Create NOX linear system object
    Teuchos::RCP<const Epetra_Vector> det_u = model->get_x_init();
    Teuchos::RCP<Epetra_Operator> det_A = model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> det_iReq = det_nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> det_iJac = det_nox_interface;
    Teuchos::ParameterList det_printParams;
    det_printParams.set("MyPID", MyPID); 
    det_printParams.set("Output Precision", 3);
    det_printParams.set("Output Processor", 0);
    det_printParams.set("Output Information", NOX::Utils::Error);
    Teuchos::ParameterList det_lsParams, det_ML;
    Teuchos::RCP<NOX::Epetra::LinearSystem> det_linsys;
    ML_Epetra::SetDefaults("SA", det_ML);
    det_ML.set("ML output", 0);
    det_ML.set("max levels",5);
    det_ML.set("increasing or decreasing","increasing");
    det_ML.set("aggregation: type", "Uncoupled");
    det_ML.set("smoother: type","ML symmetric Gauss-Seidel");
    det_ML.set("smoother: sweeps",2);
    det_ML.set("smoother: pre or post", "both");
    det_ML.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
    det_ML.set("coarse: type","Amesos-KLU");
#else
    det_ML.set("coarse: type","Jacobi");
#endif
    if (inner_krylov_solver == AZTECOO) {
      if (inner_krylov_method == GMRES)
	det_lsParams.set("Aztec Solver", "GMRES");  
      else if (inner_krylov_method == CG)
	det_lsParams.set("Aztec Solver", "CG");  
      det_lsParams.set("Max Iterations", inner_its);
      det_lsParams.set("Size of Krylov Subspace", 100);
      det_lsParams.set("Tolerance", inner_tol); 
      det_lsParams.set("Output Frequency", 0);
      det_lsParams.set("Preconditioner", "ML");    
      det_lsParams.set("Zero Initial Guess", true);    
      det_lsParams.sublist("ML") = det_ML;
      det_linsys =
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(det_printParams, 
							  det_lsParams,
							  det_iReq, 
							  det_iJac, 
							  det_A,
							  *det_u));
    }
    else if (inner_krylov_solver == BELOS) {
      Teuchos::ParameterList& stratParams = 
	det_lsParams.sublist("Stratimikos");
      stratParams.set("Linear Solver Type", "Belos");
      Teuchos::ParameterList& belosParams = 
	stratParams.sublist("Linear Solver Types").sublist("Belos");
      Teuchos::ParameterList* belosSolverParams;
      if (inner_krylov_method == GMRES || inner_krylov_method == FGMRES) {
	belosParams.set("Solver Type","Block GMRES");
	belosSolverParams = 
	  &(belosParams.sublist("Solver Types").sublist("Block GMRES"));
	if (outer_krylov_method == FGMRES)
	  belosSolverParams->set("Flexible Gmres", true);
      }
      else if (inner_krylov_method == CG) {
	belosParams.set("Solver Type","Block CG");
	belosSolverParams = 
	  &(belosParams.sublist("Solver Types").sublist("Block CG"));
      }
      // else if (inner_krylov_method == RGMRES) {
      // 	belosParams.set("Solver Type","GCRODR");
      // 	belosSolverParams = 
      // 	  &(belosParams.sublist("Solver Types").sublist("GCRODR"));
      // }
      belosSolverParams->set("Convergence Tolerance", inner_tol);
      belosSolverParams->set("Maximum Iterations", inner_its);
      belosSolverParams->set("Output Frequency",0);
      belosSolverParams->set("Output Style",1);
      belosSolverParams->set("Verbosity",0);
      Teuchos::ParameterList& verbParams = belosParams.sublist("VerboseObject");
      verbParams.set("Verbosity Level", "none");
      stratParams.set("Preconditioner Type", "ML");
      stratParams.sublist("Preconditioner Types").sublist("ML").sublist("ML Settings") = det_ML;
      det_linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		       det_printParams, det_lsParams, det_iJac, 
		       det_A, *det_u));
    }
   
    // Set up stochastic parameters
    Epetra_LocalMap p_sg_map(num_KL, 0, *Comm);
    Teuchos::Array<Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> >sg_p_init(1);
    sg_p_init[0]= Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, p_sg_map));
    for (int i=0; i<num_KL; i++) {
      sg_p_init[0]->term(i,0)[i] = 0.0;
      sg_p_init[0]->term(i,1)[i] = 1.0;
    }

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init = 
      Teuchos::rcp(new Stokhos::EpetraVectorOrthogPoly(basis, 
						       *(model->get_x_map())));
    sg_x_init->init(0.0);
    
    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& sgOpParams = 
      sgParams->sublist("SG Operator");
    Teuchos::ParameterList& sgPrecParams = 
      sgParams->sublist("SG Preconditioner");

    if (!nonlinear_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    sgOpParams.set("Operator Method", "Matrix Free");
    if (precMethod == MEAN)  {
      sgPrecParams.set("Preconditioner Method", "Mean-based");
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
	sgPrecParams.sublist("Mean Preconditioner Parameters");
      precParams = det_ML;
    }
    else if(precMethod == GS) {
      sgPrecParams.set("Preconditioner Method", "Gauss-Seidel");
      sgPrecParams.sublist("Deterministic Solver Parameters") = det_lsParams;
      sgPrecParams.set("Deterministic Solver", det_linsys);
      sgPrecParams.set("Max Iterations", gs_prec_its);
      sgPrecParams.set("Tolerance", gs_prec_tol);
    }
    else if (precMethod == AGS)  {
      sgPrecParams.set("Preconditioner Method", "Approximate Gauss-Seidel");
      if (outer_krylov_method == CG)
	sgPrecParams.set("Symmetric Gauss-Seidel", true);
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
	sgPrecParams.sublist("Mean Preconditioner Parameters");
      precParams = det_ML;
    }
    else if (precMethod == AJ)  {
      sgPrecParams.set("Preconditioner Method", "Approximate Jacobi");
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
        sgPrecParams.sublist("Mean Preconditioner Parameters");
      precParams = det_ML;
      Teuchos::ParameterList& jacobiOpParams =
	sgPrecParams.sublist("Jacobi SG Operator");
      jacobiOpParams.set("Only Use Linear Terms", true);
    }
    else if (precMethod == KP)  {
      sgPrecParams.set("Preconditioner Method", "Kronecker Product");
      sgPrecParams.set("Only Use Linear Terms", true);
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& meanPrecParams =
        sgPrecParams.sublist("Mean Preconditioner Parameters");
      meanPrecParams = det_ML;
      sgPrecParams.set("G Preconditioner Type", "Ifpack");
      Teuchos::ParameterList& GPrecParams =
        sgPrecParams.sublist("G Preconditioner Parameters");
      if (outer_krylov_method == GMRES || outer_krylov_method == FGMRES)
	GPrecParams.set("Ifpack Preconditioner", "ILUT");
      if (outer_krylov_method == CG)
	GPrecParams.set("Ifpack Preconditioner", "ICT");
      GPrecParams.set("Overlap", 1);
      GPrecParams.set("fact: drop tolerance", 1e-4);
      GPrecParams.set("fact: ilut level-of-fill", 1.0);
      GPrecParams.set("schwarz: combine mode", "Add");
    }
    else
      TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner method " << precMethod
			 << "." << std::endl);

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
                                                 expansion, Cijk, sgParams,
                                                 Comm, sg_x_init, sg_p_init,
						 scaleOP));

     // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface =
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));

    // Create NOX stochastic linear system object
    Teuchos::RCP<const Epetra_Vector> u = sg_model->get_x_init();
    Teuchos::RCP<const Epetra_Map> base_map = model->get_x_map();
    Teuchos::RCP<const Epetra_Map> sg_map = sg_model->get_x_map();
    Teuchos::RCP<Epetra_Operator> A = sg_model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;

    // Build linear solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
    if (solve_method==SG_KRYLOV) {
      Teuchos::RCP<Epetra_Operator> M = sg_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = 
	nox_interface;
      if (outer_krylov_solver == AZTECOO) {
	if (outer_krylov_method == GMRES)
	  lsParams.set("Aztec Solver", "GMRES");
	else if (outer_krylov_method == CG)
	  lsParams.set("Aztec Solver", "CG");
	lsParams.set("Max Iterations", outer_its);
	lsParams.set("Size of Krylov Subspace", 100);
	lsParams.set("Tolerance", outer_tol); 
	lsParams.set("Output Frequency", 1);
	lsParams.set("Preconditioner", "ML");
	lsParams.set("Zero Initial Guess", true);
	lsParams.set("Preconditioner", "User Defined");
	linsys = 
	  Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(
			 printParams, lsParams, iJac, A, iPrec, M, *u));
      }
      else if (outer_krylov_solver == BELOS){
	stratParams.set("Linear Solver Type", "Belos");
	Teuchos::ParameterList& belosParams = 
	  stratParams.sublist("Linear Solver Types").sublist("Belos");
	Teuchos::ParameterList* belosSolverParams;
	if (outer_krylov_method == GMRES || outer_krylov_method == FGMRES) {
	  belosParams.set("Solver Type","Block GMRES");
	  belosSolverParams = 
	    &(belosParams.sublist("Solver Types").sublist("Block GMRES"));
	  if (outer_krylov_method == FGMRES)
	    belosSolverParams->set("Flexible Gmres", true);
	}
	else if (outer_krylov_method == CG) {
	  belosParams.set("Solver Type","Block CG");
	  belosSolverParams = 
	    &(belosParams.sublist("Solver Types").sublist("Block CG"));
	}
	// else if (inner_krylov_method == RGMRES) {
	//   belosParams.set("Solver Type","GCRODR");
	//   belosSolverParams = 
	//     &(belosParams.sublist("Solver Types").sublist("GCRODR"));
	// }
	belosSolverParams->set("Convergence Tolerance", outer_tol);
	belosSolverParams->set("Maximum Iterations", outer_its);
	belosSolverParams->set("Output Frequency",1);
	belosSolverParams->set("Output Style",1);
	belosSolverParams->set("Verbosity",33);
	stratLinSolParams.set("Preconditioner", "User Defined");
	linsys = 
	  Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			 printParams, stratLinSolParams, iJac, A, iPrec, M, 
			 *u, true));
      }
    }
    else if (solve_method==SG_GS) {
      lsParams.sublist("Deterministic Solver Parameters") = det_lsParams;
      lsParams.set("Max Iterations", outer_its);
      lsParams.set("Tolerance", outer_tol);
      linsys =
	Teuchos::rcp(new NOX::Epetra::LinearSystemSGGS(
		       printParams, lsParams, det_linsys, Cijk,
		       iReq, iJac, A, base_map, sg_map));
    }
    else {
      lsParams.sublist("Deterministic Solver Parameters") = det_lsParams;
      lsParams.set("Max Iterations", outer_its);
      lsParams.set("Tolerance", outer_tol);
      Teuchos::ParameterList& jacobiOpParams =
	lsParams.sublist("Jacobi SG Operator");
      jacobiOpParams.set("Only Use Linear Terms", true);
      linsys =
	Teuchos::rcp(new NOX::Epetra::LinearSystemSGJacobi(
		       printParams, lsParams, det_linsys, Cijk,
		       iReq, iJac, A, base_map, sg_map));
    }

    // Build NOX group
    Teuchos::RCP<NOX::Epetra::Group> grp = 
      Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));
    
    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::Generic> statusTests =
      NOX::StatusTest::buildStatusTests(statusParams, utils);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(grp, statusTests, noxParams);

    // Solve the system
    NOX::StatusTest::StatusType status;
    {
      TEUCHOS_FUNC_TIME_MONITOR("Total Solve Time");
      status = solver->solve();
    }

    // Get final solution
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // Save final solution to file
    EpetraExt::VectorToMatrixMarketFile("nox_solver_stochastic_solution.mm", 
					finalSolution);

    // Save mean and variance to file
    Stokhos::EpetraVectorOrthogPoly sg_x_poly(basis, View, 
					      *(model->get_x_map()), 
					      finalSolution);
    Epetra_Vector mean(*(model->get_x_map()));
    Epetra_Vector std_dev(*(model->get_x_map()));
    sg_x_poly.computeMean(mean);
    sg_x_poly.computeStandardDeviation(std_dev);
    EpetraExt::VectorToMatrixMarketFile("mean_gal.mm", mean);
    EpetraExt::VectorToMatrixMarketFile("std_dev_gal.mm", std_dev);
      
    // Evaluate SG responses at SG parameters
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_model->createOutArgs();
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(1))));
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(Teuchos::rcp(&finalSolution,false));
    sg_outArgs.set_g(1, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print mean and standard deviation of response
    Stokhos::EpetraVectorOrthogPoly sg_g_poly(basis, View, 
					      *(model->get_g_map(0)), *sg_g);
    Epetra_Vector g_mean(*(model->get_g_map(0)));
    Epetra_Vector g_std_dev(*(model->get_g_map(0)));
    sg_g_poly.computeMean(g_mean);
    sg_g_poly.computeStandardDeviation(g_std_dev);
    // std::cout << "\nResponse Expansion = " << std::endl;
    // std::cout.precision(12);
    // sg_g_poly.print(std::cout);
    std::cout << "\nResponse Mean =      " << std::endl << g_mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << g_std_dev << std::endl;

    if (status == NOX::StatusTest::Converged) 
      utils.out() << "Example Passed!" << std::endl;

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
