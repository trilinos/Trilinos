// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Teuchos_Assert.hpp"

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
enum Krylov_Method { GMRES, CG, FGMRES, RGMRES };
const int num_krylov_method = 4;
const Krylov_Method krylov_method_values[] = { GMRES, CG, FGMRES, RGMRES };
const char *krylov_method_names[] = { "GMRES", "CG", "FGMRES", "RGMRES" };

// Krylov operator approaches
enum SG_Op { MATRIX_FREE, KL_MATRIX_FREE, KL_REDUCED_MATRIX_FREE, FULLY_ASSEMBLED };
const int num_sg_op = 4;
const SG_Op sg_op_values[] = { MATRIX_FREE, KL_MATRIX_FREE, KL_REDUCED_MATRIX_FREE, FULLY_ASSEMBLED };
const char *sg_op_names[] = { "Matrix-Free",
			      "KL-Matrix-Free",
			      "KL-Reduced-Matrix-Free",
			      "Fully-Assembled" };

// Krylov preconditioning approaches
enum SG_Prec { MEAN, GS, AGS, AJ, ASC, KP, NONE };
const int num_sg_prec = 7;
const SG_Prec sg_prec_values[] = { MEAN, GS, AGS, AJ, ASC, KP, NONE };
const char *sg_prec_names[] = { "Mean-Based", 
				"Gauss-Seidel", 
				"Approx-Gauss-Seidel", 
				"Approx-Jacobi", 
				"Approx-Schur-Complement", 
				"Kronecker-Product",
				"None" };

// Random field types
enum SG_RF { UNIFORM, RYS, LOGNORMAL };
const int num_sg_rf = 3;
const SG_RF sg_rf_values[] = { UNIFORM, RYS, LOGNORMAL };
const char *sg_rf_names[] = { "Uniform", "Rys", "Log-Normal" };

int main(int argc, char *argv[]) {   

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
  Teuchos::RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
  globalComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  globalComm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  int MyPID = globalComm->MyPID();

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

    int num_spatial_procs = -1;
    CLP.setOption("num_spatial_procs", &num_spatial_procs, "Number of spatial processors (set -1 for all available procs)");

    bool rebalance_stochastic_graph = false;
    CLP.setOption("rebalance", "no-rebalance", &rebalance_stochastic_graph, 
		  "Rebalance parallel stochastic graph (requires Isorropia)");

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

    int p = 3;
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

    SG_Op opMethod = MATRIX_FREE;
    CLP.setOption("sg_operator_method", &opMethod, 
		  num_sg_op, sg_op_values, sg_op_names,
		  "Operator method");

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
		<< "\tnum_spatial_procs   = " << num_spatial_procs << std::endl
		<< "\trebalance           = " << rebalance_stochastic_graph
		<< std::endl
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
		<< "\touter_tol           = " << outer_tol << std::endl
		<< "\touter_its           = " << outer_its << std::endl
		<< "\tinner_krylov_method = " 
		<< krylov_method_names[inner_krylov_method] << std::endl
		<< "\tinner_krylov_solver = " 
		<< krylov_solver_names[inner_krylov_solver] << std::endl
		<< "\tinner_tol           = " << inner_tol << std::endl
		<< "\tinner_its           = " << inner_its << std::endl
		<< "\tsg_operator_method  = " << sg_op_names[opMethod] 
		<< std::endl
		<< "\tsg_prec_method      = " << sg_prec_names[precMethod] 
		<< std::endl
		<< "\tgs_prec_tol         = " << gs_prec_tol << std::endl
		<< "\tgs_prec_its         = " << gs_prec_its << std::endl;
    }

    bool nonlinear_expansion = false;
    if (randField == UNIFORM || randField == RYS)
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
      else if (randField == RYS)
        bases[i] = Teuchos::rcp(new Stokhos::RysBasis<int,double>(p,weightCut,normalize_basis));
      else if (randField == LOGNORMAL)      
        bases[i] = Teuchos::rcp(new Stokhos::HermiteBasis<int,double>(p,normalize_basis));

    //  bases[i] = Teuchos::rcp(new Stokhos::DiscretizedStieltjesBasis<int,double>("beta",p,&uniform_weight,-weightCut,weightCut,true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (nonlinear_expansion)
      Cijk = basis->computeTripleProductTensor();
    else
      Cijk = basis->computeLinearTripleProductTensor();
    Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > expansion = 
      Teuchos::rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis,
									 Cijk));
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    Teuchos::ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    parallelParams.set("Rebalance Stochastic Graph", 
		       rebalance_stochastic_graph);
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
      Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
					     parallelParams));
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = 
      sg_parallel_data->getMultiComm();
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      sg_parallel_data->getSpatialComm();
    
    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model =
      Teuchos::rcp(new twoD_diffusion_ME(app_comm, n, num_KL, sigma, 
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
    
    Teuchos::ParameterList det_lsParams;
    Teuchos::ParameterList& det_stratParams = 
      det_lsParams.sublist("Stratimikos");
    if (inner_krylov_solver == AZTECOO) {
      det_stratParams.set("Linear Solver Type", "AztecOO");
      Teuchos::ParameterList& aztecOOParams = 
	det_stratParams.sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve");
      Teuchos::ParameterList& aztecOOSettings =
	aztecOOParams.sublist("AztecOO Settings");
      if (inner_krylov_method == GMRES) {
	aztecOOSettings.set("Aztec Solver","GMRES");
      }
      else if (inner_krylov_method == CG) {
	aztecOOSettings.set("Aztec Solver","CG");
      }
      aztecOOSettings.set("Output Frequency", 0);
      aztecOOSettings.set("Size of Krylov Subspace", 100);
      aztecOOParams.set("Max Iterations", inner_its);
      aztecOOParams.set("Tolerance", inner_tol);
      Teuchos::ParameterList& verbParams = 
	det_stratParams.sublist("Linear Solver Types").sublist("AztecOO").sublist("VerboseObject");
      verbParams.set("Verbosity Level", "none");
    }
    else if (inner_krylov_solver == BELOS) {
      det_stratParams.set("Linear Solver Type", "Belos");
      Teuchos::ParameterList& belosParams = 
	det_stratParams.sublist("Linear Solver Types").sublist("Belos");
      Teuchos::ParameterList* belosSolverParams = NULL;
      if (inner_krylov_method == GMRES || inner_krylov_method == FGMRES) {
	belosParams.set("Solver Type","Block GMRES");
	belosSolverParams = 
	  &(belosParams.sublist("Solver Types").sublist("Block GMRES"));
	if (inner_krylov_method == FGMRES)
	  belosSolverParams->set("Flexible Gmres", true);
      }
      else if (inner_krylov_method == CG) {
	belosParams.set("Solver Type","Block CG");
	belosSolverParams = 
	  &(belosParams.sublist("Solver Types").sublist("Block CG"));
      }
      else if (inner_krylov_method == RGMRES) {
      	belosParams.set("Solver Type","GCRODR");
      	belosSolverParams = 
      	  &(belosParams.sublist("Solver Types").sublist("GCRODR"));
      }
      belosSolverParams->set("Convergence Tolerance", inner_tol);
      belosSolverParams->set("Maximum Iterations", inner_its);
      belosSolverParams->set("Output Frequency",0);
      belosSolverParams->set("Output Style",1);
      belosSolverParams->set("Verbosity",0);
      Teuchos::ParameterList& verbParams = belosParams.sublist("VerboseObject");
      verbParams.set("Verbosity Level", "none");
    }
    det_stratParams.set("Preconditioner Type", "ML");
    Teuchos::ParameterList& det_ML = 
      det_stratParams.sublist("Preconditioner Types").sublist("ML").sublist("ML Settings");
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
    Teuchos::RCP<NOX::Epetra::LinearSystem> det_linsys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		     det_printParams, det_lsParams, det_iJac, 
		     det_A, *det_u));
    
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
    if (opMethod == MATRIX_FREE)
      sgOpParams.set("Operator Method", "Matrix Free");
    else if (opMethod == KL_MATRIX_FREE)
      sgOpParams.set("Operator Method", "KL Matrix Free");
    else if (opMethod == KL_REDUCED_MATRIX_FREE) {
      sgOpParams.set("Operator Method", "KL Reduced Matrix Free");
      if (randField == UNIFORM || randField == RYS)
	sgOpParams.set("Number of KL Terms", num_KL);
      else
	sgOpParams.set("Number of KL Terms", basis->size());
      sgOpParams.set("KL Tolerance", outer_tol);
      sgOpParams.set("Sparse 3 Tensor Drop Tolerance", outer_tol);
      sgOpParams.set("Do Error Tests", true);
    }
    else if (opMethod == FULLY_ASSEMBLED)
      sgOpParams.set("Operator Method", "Fully Assembled");
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown operator method " << opMethod
			 << "." << std::endl);
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
    else if (precMethod == ASC)  {
      sgPrecParams.set("Preconditioner Method", "Approximate Schur Complement");
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams =
	sgPrecParams.sublist("Mean Preconditioner Parameters");
      precParams = det_ML;
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
    else if (precMethod == NONE)  {
      sgPrecParams.set("Preconditioner Method", "None");
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner method " << precMethod
			 << "." << std::endl);

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
                                                 expansion, sg_parallel_data, 
						 sgParams, scaleOP));
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_model->createOutArgs();

    // Set up stochastic parameters
    // The current implementation of the model doesn't actually use these 
    // values, but is hard-coded to certain uncertainty models
    Teuchos::Array<double> point(num_KL, 1.0);
    Teuchos::Array<double> basis_vals(sz);
    basis->evaluateBases(point, basis_vals);
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_p_init =
      sg_model->create_p_sg(0);
    for (int i=0; i<num_KL; i++) {
      sg_p_init->term(i,0)[i] = 0.0;
      sg_p_init->term(i,1)[i] = 1.0 / basis_vals[i+1];
    }
    sg_model->set_p_sg_init(0, *sg_p_init);

    // Setup stochastic initial guess
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_init = 
      sg_model->create_x_sg();
    sg_x_init->init(0.0);
    sg_model->set_x_sg_init(*sg_x_init);

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
      bool has_M = 
	sg_outArgs.supports(EpetraExt::ModelEvaluator::OUT_ARG_WPrec);
      Teuchos::RCP<Epetra_Operator> M;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec;
      if (has_M) {
	M = sg_model->create_WPrec()->PrecOp;
	iPrec = nox_interface;
      }
      stratParams.set("Preconditioner Type", "None");
      if (outer_krylov_solver == AZTECOO) {
	stratParams.set("Linear Solver Type", "AztecOO");
	Teuchos::ParameterList& aztecOOParams = 
	  stratParams.sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve");
	Teuchos::ParameterList& aztecOOSettings =
	  aztecOOParams.sublist("AztecOO Settings");
	if (outer_krylov_method == GMRES) {
	  aztecOOSettings.set("Aztec Solver","GMRES");
	}
	else if (outer_krylov_method == CG) {
	  aztecOOSettings.set("Aztec Solver","CG");
	}
	aztecOOSettings.set("Output Frequency", 1);
	aztecOOSettings.set("Size of Krylov Subspace", 100);
	aztecOOParams.set("Max Iterations", outer_its);
	aztecOOParams.set("Tolerance", outer_tol);
	stratLinSolParams.set("Preconditioner", "User Defined");
	if (has_M)
	  linsys = 
	    Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			   printParams, stratLinSolParams, iJac, A, iPrec, M, 
			   *u, true));
	else
	  linsys = 
	    Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			   printParams, stratLinSolParams, iJac, A, *u));
      }
      else if (outer_krylov_solver == BELOS){
	stratParams.set("Linear Solver Type", "Belos");
	Teuchos::ParameterList& belosParams = 
	  stratParams.sublist("Linear Solver Types").sublist("Belos");
	Teuchos::ParameterList* belosSolverParams = NULL;
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
	else if (inner_krylov_method == RGMRES) {
	  belosParams.set("Solver Type","GCRODR");
	  belosSolverParams = 
	    &(belosParams.sublist("Solver Types").sublist("GCRODR"));
	}
	belosSolverParams->set("Convergence Tolerance", outer_tol);
	belosSolverParams->set("Maximum Iterations", outer_its);
	belosSolverParams->set("Output Frequency",1);
	belosSolverParams->set("Output Style",1);
	belosSolverParams->set("Verbosity",33);
	stratLinSolParams.set("Preconditioner", "User Defined");
	if (has_M)
	  linsys = 
	    Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			   printParams, stratLinSolParams, iJac, A, iPrec, M, 
			   *u, true));
	else
	  linsys = 
	    Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
			   printParams, stratLinSolParams, iJac, A, *u));
	  
      }
    }
    else if (solve_method==SG_GS) {
      lsParams.sublist("Deterministic Solver Parameters") = det_lsParams;
      lsParams.set("Max Iterations", outer_its);
      lsParams.set("Tolerance", outer_tol);
      linsys =
	Teuchos::rcp(new NOX::Epetra::LinearSystemSGGS(
		       printParams, lsParams, det_linsys, iReq, iJac, 
		       basis, sg_parallel_data, A, base_map, sg_map));
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
		       printParams, lsParams, det_linsys, iReq, iJac, 
		       basis, sg_parallel_data, A, base_map, sg_map));
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
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_x_poly = 
      sg_model->create_x_sg(View, &finalSolution);
    Epetra_Vector mean(*(model->get_x_map()));
    Epetra_Vector std_dev(*(model->get_x_map()));
    sg_x_poly->computeMean(mean);
    sg_x_poly->computeStandardDeviation(std_dev);
    EpetraExt::VectorToMatrixMarketFile("mean_gal.mm", mean);
    EpetraExt::VectorToMatrixMarketFile("std_dev_gal.mm", std_dev);
      
    // Evaluate SG responses at SG parameters
    Teuchos::RCP<const Epetra_Vector> sg_p = sg_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> sg_g = 
      Teuchos::rcp(new Epetra_Vector(*(sg_model->get_g_map(0))));
    sg_inArgs.set_p(1, sg_p);
    sg_inArgs.set_x(Teuchos::rcp(&finalSolution,false));
    sg_outArgs.set_g(0, sg_g);
    sg_model->evalModel(sg_inArgs, sg_outArgs);

    // Print mean and standard deviation of response
    Teuchos::RCP<Stokhos::EpetraVectorOrthogPoly> sg_g_poly =
      sg_model->create_g_sg(0, View, sg_g.get());
    Epetra_Vector g_mean(*(model->get_g_map(0)));
    Epetra_Vector g_std_dev(*(model->get_g_map(0)));
    sg_g_poly->computeMean(g_mean);
    sg_g_poly->computeStandardDeviation(g_std_dev);
    std::cout.precision(16);
    // std::cout << "\nResponse Expansion = " << std::endl;
    // std::cout.precision(12);
    // sg_g_poly->print(std::cout);
    std::cout << std::endl;
    std::cout << "Response Mean =      " << std::endl << g_mean << std::endl;
    std::cout << "Response Std. Dev. = " << std::endl << g_std_dev << std::endl;

    if (status == NOX::StatusTest::Converged && MyPID == 0) 
      utils.out() << "Example Passed!" << std::endl;

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
