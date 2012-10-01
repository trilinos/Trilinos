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

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Timing utilities
#include "Teuchos_TimeMonitor.hpp"

// I/O utilities
#include "EpetraExt_VectorOut.h"

int main(int argc, char *argv[]) {
  int n = 32;                        // spatial discretization (per dimension)
  int num_KL = 2;                    // number of KL terms
  int p = 3;                         // polynomial order
  double mu = 0.1;                   // mean of exponential random field
  double s = 0.2;                    // std. dev. of exponential r.f.
  bool nonlinear_expansion = false;  // nonlinear expansion of diffusion coeff
                                     // (e.g., log-normal)
  bool matrix_free = true;           // use matrix-free stochastic operator
  bool symmetric = false;            // use symmetric formulation

  double g_mean_exp = 0.172988;      // expected response mean
  double g_std_dev_exp = 0.0380007;  // expected response std. dev.
  double g_tol = 1e-6;               // tolerance on determining success

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  int MyPID;

  try {

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
    
    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
    for (int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(p, true));
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
    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    int num_spatial_procs = -1;
    if (argc > 1)
      num_spatial_procs = std::atoi(argv[1]);
    Teuchos::ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data =
      Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
					     parallelParams));
    Teuchos::RCP<const EpetraExt::MultiComm> sg_comm = 
      sg_parallel_data->getMultiComm();
    Teuchos::RCP<const Epetra_Comm> app_comm = 
      sg_parallel_data->getSpatialComm();

    // Create application
    Teuchos::RCP<twoD_diffusion_ME> model =
      Teuchos::rcp(new twoD_diffusion_ME(app_comm, n, num_KL, mu, s, basis, 
					 nonlinear_expansion, symmetric));
    
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
    if (matrix_free) {
      sgOpParams.set("Operator Method", "Matrix Free");
      sgPrecParams.set("Preconditioner Method", "Approximate Gauss-Seidel");
      sgPrecParams.set("Symmetric Gauss-Seidel", symmetric);
      sgPrecParams.set("Mean Preconditioner Type", "ML");
      Teuchos::ParameterList& precParams = 
      	sgPrecParams.sublist("Mean Preconditioner Parameters");
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
    }
    else {
      sgOpParams.set("Operator Method", "Fully Assembled");
      sgPrecParams.set("Preconditioner Method", "None");
    }

   // Create stochastic Galerkin model evaluator
    Teuchos::RCP<Stokhos::SGModelEvaluator> sg_model =
      Teuchos::rcp(new Stokhos::SGModelEvaluator(model, basis, Teuchos::null,
                                                 expansion, sg_parallel_data, 
						 sgParams));

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

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    if (symmetric)
      lsParams.set("Aztec Solver", "CG");
    else
      lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 1000);
    lsParams.set("Size of Krylov Subspace", 100);
    lsParams.set("Tolerance", 1e-12); 
    lsParams.set("Output Frequency", 1);
    if (matrix_free)
      lsParams.set("Preconditioner", "User Defined");
    else {
      lsParams.set("Preconditioner", "ML");
      Teuchos::ParameterList& precParams = 
	lsParams.sublist("ML");
      ML_Epetra::SetDefaults("DD", precParams);
      lsParams.set("Write Linear System", false);
    }

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
       Teuchos::rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));

    // Create NOX linear system object
    Teuchos::RCP<const Epetra_Vector> u = sg_model->get_x_init();
    Teuchos::RCP<Epetra_Operator> A = sg_model->create_W();
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
    if (matrix_free) {
      Teuchos::RCP<Epetra_Operator> M = sg_model->create_WPrec()->PrecOp;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = nox_interface;
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iJac, A, iPrec, M,
							  *u));
    }
    else {
      linsys = 
	Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							  iReq, iJac, A, 
							  *u));
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

    // Save final solution to file
    EpetraExt::VectorToMatrixMarketFile("nox_stochastic_solution.mm", 
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
    EpetraExt::ModelEvaluator::InArgs sg_inArgs = sg_model->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs sg_outArgs = 
      sg_model->createOutArgs();
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

    // Determine if example passed
    bool passed = false;
    if (status == NOX::StatusTest::Converged &&
	std::abs(g_mean[0]-g_mean_exp) < g_tol &&
	std::abs(g_std_dev[0]-g_std_dev_exp) < g_tol)
      passed = true;
    if (MyPID == 0) {
      if (passed)
	std::cout << "Example Passed!" << std::endl;
      else
	std::cout << "Example Failed!" << std::endl;
    }

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
