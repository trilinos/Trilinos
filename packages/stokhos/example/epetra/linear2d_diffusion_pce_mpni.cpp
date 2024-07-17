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
#include "BelosTypes.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "NOX_Epetra_LinearSystem_MPBD.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

// I/O utilities
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_BlockUtility.h"

int main(int argc, char *argv[]) {
  int n = 32;                        // spatial discretization (per dimension)
  int num_KL = 2;                    // number of KL terms
  int p = 3;                         // polynomial order
  double mu = 0.1;                   // mean of exponential random field
  double s = 0.2;                    // std. dev. of exponential r.f.
  bool nonlinear_expansion = false;  // nonlinear expansion of diffusion coeff
                                     // (e.g., log-normal)
  bool symmetric = true;            // use symmetric formulation
  bool use_solver = true;
  std::string solver_type = "RGMRES";

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;

  try {

    Epetra_Object::SetTracebackMode(1);

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create a communicator for Epetra objects
    Teuchos::RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = Teuchos::rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Create Stochastic Galerkin basis and quadrature
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p, true));
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = 
      Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases,
		     1e-12));
    // Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
    //   Teuchos::rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > quad = 
      Teuchos::rcp(new Stokhos::SparseGridQuadrature<int,double>(basis));
    const Teuchos::Array< Teuchos::Array<double> >& quad_points = 
      quad->getQuadPoints();
    const Teuchos::Array<double>& quad_weights = quad->getQuadWeights();
    const Teuchos::Array< Teuchos::Array<double> >& quad_values =
      quad->getBasisAtQuadPoints();
    int sz = basis->size();
    int num_mp = quad_weights.size();
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create multi-point parallel distribution
    int num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    Teuchos::RCP<const EpetraExt::MultiComm> multi_comm = 
      Stokhos::buildMultiComm(*globalComm, num_mp, num_spatial_procs);
    Teuchos::RCP<const Epetra_Comm> mp_comm = 
      Stokhos::getStochasticComm(multi_comm);
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      Stokhos::getSpatialComm(multi_comm);
    Teuchos::RCP<const Epetra_Map> mp_block_map = 
      Teuchos::rcp(new Epetra_Map(num_mp, 0, *mp_comm));
    int num_my_mp = mp_block_map->NumMyElements();
    //mp_block_map->Print(std::cout);

    Teuchos::RCP<Teuchos::ParameterList> detPrecParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    detPrecParams->set("Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      detPrecParams->sublist("Preconditioner Parameters");
    precParams.set("default values", "SA");
    precParams.set("ML output", 0);
    precParams.set("max levels",5);
    precParams.set("increasing or decreasing","increasing");
    precParams.set("aggregation: type", "Uncoupled");
    precParams.set("smoother: type","ML symmetric Gauss-Seidel");
    precParams.set("smoother: sweeps",2);
    precParams.set("smoother: pre or post", "both");
    precParams.set("coarse: max size", 200);
#ifdef HAVE_ML_AMESOS
    precParams.set("coarse: type","Amesos-KLU");
#else
    precParams.set("coarse: type","Jacobi");
#endif

    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model = 
      Teuchos::rcp(new twoD_diffusion_ME(app_comm, n, num_KL, mu, s, basis, 
					 nonlinear_expansion, symmetric,
					 detPrecParams));
    
    // Setup multi-point algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> mpParams = 
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& mpPrecParams = 
      mpParams->sublist("MP Preconditioner");
    mpPrecParams.set("Preconditioner Method", "Block Diagonal");
    mpPrecParams.set("MP Preconditioner Type", "ML");
    Teuchos::ParameterList& pointPrecParams = 
      mpPrecParams.sublist("MP Preconditioner Parameters");
    pointPrecParams = precParams;

    // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::MPModelEvaluator> mp_model =
      Teuchos::rcp(new Stokhos::MPModelEvaluator(model, multi_comm, 
						 mp_block_map, mpParams));

    // Set up multi-point parameters
    Teuchos::RCP<Stokhos::ProductEpetraVector> mp_p_init =
      mp_model->create_p_mp(0);
    int my_mp_begin = mp_block_map->MinMyGID();
    for (int j=0; j<num_my_mp; j++) {
      for (int i=0; i<num_KL; i++) {
	(*mp_p_init)[j][i] = quad_points[j+my_mp_begin][i];
      }
    }
    mp_model->set_p_mp_init(0, *mp_p_init);

    // Setup multi-point initial guess
    Teuchos::RCP<Stokhos::ProductEpetraVector> mp_x_init = 
      mp_model->create_x_mp();
    mp_x_init->init(0.0);
    mp_model->set_x_mp_init(*mp_x_init);

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

    // Alternative linear solver list for Stratimikos
    Teuchos::ParameterList stratLinSolParams;
    Teuchos::ParameterList& stratParams = 
      stratLinSolParams.sublist("Stratimikos");

    // Sublist for linear solver for the Newton method
    stratParams.set("Linear Solver Type", "Belos");
    Teuchos::ParameterList& belosParams = 
      stratParams.sublist("Linear Solver Types").sublist("Belos");
    Teuchos::ParameterList* belosSolverParams = NULL;
    if (solver_type == "GMRES") {
      belosParams.set("Solver Type","Block GMRES");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("Block GMRES"));
    }
    else if (solver_type == "CG") {
      belosParams.set("Solver Type","Block CG");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("Block CG"));
    }
    else if (solver_type == "RGMRES") {
      belosParams.set("Solver Type","GCRODR");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("GCRODR"));
      belosSolverParams->set("Num Recycled Blocks", 20);
    }
    else if (solver_type == "RCG") {
      belosParams.set("Solver Type","RCG");
      belosSolverParams = 
	&(belosParams.sublist("Solver Types").sublist("RCG"));
      Teuchos::RCP<const Teuchos::ParameterList> ortho_params = 
	  Teuchos::rcp(new Teuchos::ParameterList);
      belosSolverParams->set("Num Recycled Blocks", 10);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			 "Unknown solver type " << solver_type);
    belosSolverParams->set("Convergence Tolerance", 1e-12);
    belosSolverParams->set("Maximum Iterations", 1000);
    if (solver_type != "CG")
      belosSolverParams->set("Num Blocks", 100);
    belosSolverParams->set("Block Size", 1);
    belosSolverParams->set("Output Frequency",1);
    belosSolverParams->set("Output Style",1);
    //belosSolverParams->set("Verbosity",33);
    belosSolverParams->set("Verbosity", 
			   Belos::Errors + 
			   Belos::Warnings +
			   Belos::StatusTestDetails);
    stratLinSolParams.set("Preconditioner", "User Defined");
    Teuchos::ParameterList& verboseParams = 
      belosParams.sublist("VerboseObject");
    verboseParams.set("Verbosity Level", "medium");

    // Sublist for convergence tests
    Teuchos::ParameterList& statusParams = noxParams->sublist("Status Tests");
    statusParams.set("Test Type", "Combo");
    statusParams.set("Number of Tests", 2);
    statusParams.set("Combo Type", "OR");
    Teuchos::ParameterList& normF = statusParams.sublist("Test 0");
    normF.set("Test Type", "NormF");
    normF.set("Tolerance", 1e-10);
    normF.set("Scale Type", "Scaled");
    Teuchos::ParameterList& maxIters = statusParams.sublist("Test 1");
    maxIters.set("Test Type", "MaxIters");
    maxIters.set("Maximum Iterations", 1);

    // Create NOX interface
    Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface = 
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(mp_model));

    // Create NOX linear system object
    Teuchos::RCP<const Epetra_Vector> u = mp_model->get_x_init();
    
    Teuchos::RCP<NOX::Epetra::LinearSystem> linsys;
    Teuchos::RCP<Epetra_Operator> A = mp_model->create_W();
    Teuchos::RCP<Epetra_Operator> M = mp_model->create_WPrec()->PrecOp;
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = nox_interface;
    if (use_solver) {
      Teuchos::ParameterList& outerSolParams = 
	newtonParams.sublist("Linear Solver");
      outerSolParams.sublist("Deterministic Solver Parameters") = 
	stratLinSolParams;
      outerSolParams.set("Preconditioner Strategy", "Mean");
      Teuchos::RCP<Epetra_Operator> inner_A = model->create_W();
      Teuchos::RCP<Epetra_Operator> inner_M = model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::ModelEvaluatorInterface> inner_nox_interface = 
	Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(model));
      Teuchos::RCP<NOX::Epetra::Interface::Required> inner_iReq = 
	inner_nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> inner_iJac = 
	inner_nox_interface;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> inner_iPrec = 
	inner_nox_interface;
      Teuchos::RCP<const Epetra_Vector> inner_u = model->get_x_init();
      Teuchos::RCP<NOX::Epetra::LinearSystem> inner_linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(
		       printParams, 
		       stratLinSolParams,
		       inner_iJac, inner_A, inner_iPrec, inner_M,
		       *inner_u, true));
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemMPBD(printParams, 
						       outerSolParams,
						       inner_linsys,
						       iReq, iJac, A,
						       model->get_x_map()));
    }
    else {
      newtonParams.sublist("Stratimikos Linear Solver") = stratLinSolParams;
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams, 
							      stratLinSolParams,
							      iJac, A, iPrec, M,
							      *u, true));
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
    NOX::StatusTest::StatusType status = solver->solve();

    // Get final solution
    const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    Teuchos::RCP<const Epetra_Vector> mp_x = 
      Teuchos::rcp(&finalSolution, false);

    // Evaluate final responses
    Teuchos::RCP<const Epetra_Vector> mp_p = mp_model->get_p_init(1);
    Teuchos::RCP<Epetra_Vector> mp_g = 
      Teuchos::rcp(new Epetra_Vector(*(mp_model->get_g_map(0))));
    EpetraExt::ModelEvaluator::InArgs mp_inArgs = mp_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs mp_outArgs = mp_model->createOutArgs();
    mp_inArgs.set_x(mp_x);
    mp_inArgs.set_p(1, mp_p);
    mp_outArgs.set_g(0, mp_g);
    mp_model->evalModel(mp_inArgs, mp_outArgs);

    // Import x and g
    Teuchos::RCP<Epetra_LocalMap> mp_local_map =
      Teuchos::rcp(new Epetra_LocalMap(num_mp, 0, *globalComm));
    Teuchos::RCP<const Epetra_Map> mp_local_x_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *(model->get_x_map()), *mp_local_map, *globalComm));
    Teuchos::RCP<const Epetra_Map> mp_local_g_map = 
      Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
		     *(model->get_g_map(0)), *mp_local_map, *globalComm));
    Epetra_Import mp_x_importer(*mp_local_x_map, mp_x->Map());
    Epetra_Import mp_g_importer(*mp_local_g_map, mp_g->Map());
    Epetra_Vector mp_local_x(*mp_local_x_map);
    Epetra_Vector mp_local_g(*mp_local_g_map);
    mp_local_x.Import(*mp_x, mp_x_importer, Add);
    mp_local_g.Import(*mp_g, mp_g_importer, Add);

    // Compute PC expansions
    Stokhos::ProductEpetraVector mp_x_vec(
      mp_local_map,
      model->get_x_map(), 
      mp_local_x_map, 
      multi_comm, View, mp_local_x);
    Stokhos::ProductEpetraVector mp_g_vec(
      mp_local_map,
      model->get_g_map(0), 
      mp_local_g_map,
      multi_comm, View, mp_local_g);
    Teuchos::RCP<Epetra_LocalMap> sg_block_map = 
      Teuchos::rcp(new Epetra_LocalMap(sz, 0, *globalComm));
    Stokhos::EpetraVectorOrthogPoly sg_x(basis, sg_block_map, 
					 model->get_x_map(), multi_comm);
    Stokhos::EpetraVectorOrthogPoly sg_g(basis, sg_block_map, 
					 model->get_g_map(0), multi_comm);
    for (int i=0; i<sz; i++) {
      sg_x[i].PutScalar(0.0);
      sg_g[i].PutScalar(0.0);
      for (int j=0; j<num_mp; j++) {
	sg_x[i].Update(quad_weights[j]*quad_values[j][i], mp_x_vec[j], 1.0);
	sg_g[i].Update(quad_weights[j]*quad_values[j][i], mp_g_vec[j], 1.0);
      }
    }

    // Save SG solution to file
    EpetraExt::VectorToMatrixMarketFile("ni_stochastic_solution.mm", 
					*(sg_x.getBlockVector()));

    // Save mean and variance to file
    Epetra_Vector mean(*(model->get_x_map()));
    Epetra_Vector std_dev(*(model->get_x_map()));
    sg_x.computeMean(mean);
    sg_x.computeStandardDeviation(std_dev);
    EpetraExt::VectorToMatrixMarketFile("mean_gal.mm", mean);

    // Print mean and standard deviation of responses
    Epetra_Vector g_mean(*(model->get_g_map(0)));
    Epetra_Vector g_std_dev(*(model->get_g_map(0)));
    sg_g.computeMean(g_mean);
    sg_g.computeStandardDeviation(g_std_dev);
    // std::cout << "\nResponse Expansion = " << std::endl;
    // std::cout.precision(12);
    // sg_g.print(std::cout);
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
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
