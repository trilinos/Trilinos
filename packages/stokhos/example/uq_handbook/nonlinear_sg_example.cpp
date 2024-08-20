// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

// NOX
#include "NOX.H"
#include "NOX_Epetra.H"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"

// Utilities
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

// Our model
#include "SimpleME.hpp"

// Function to create NOX's parameter list
Teuchos::RCP<Teuchos::ParameterList>
create_nox_parameter_list(int MyPID) {
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
                  //NOX::Utils::Details +
                  NOX::Utils::LinearSolverDetails +
                  NOX::Utils::Warning +
                  NOX::Utils::Error);

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
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 100);
  lsParams.set("Size of Krylov Subspace", 100);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 10);
  lsParams.set("Preconditioner", "Ifpack");

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
  maxIters.set("Maximum Iterations", 10);

  return noxParams;
}

Teuchos::RCP<NOX::Solver::Generic>
create_nox_solver(int MyPID,
                  const Teuchos::RCP<EpetraExt::ModelEvaluator>& sg_model) {
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set up NOX parameters (implemented above)
  RCP<ParameterList> noxParams = create_nox_parameter_list(MyPID);

  // Create printing utilities
  ParameterList& printParams = noxParams->sublist("Printing");
  NOX::Utils utils(printParams);

  // Create NOX interface
  RCP<NOX::Epetra::ModelEvaluatorInterface> nox_interface =
    rcp(new NOX::Epetra::ModelEvaluatorInterface(sg_model));

  // Create NOX linear system object
  RCP<const Epetra_Vector> u = sg_model->get_x_init();
  RCP<Epetra_Operator> A = sg_model->create_W();
  RCP<NOX::Epetra::Interface::Required> iReq = nox_interface;
  RCP<NOX::Epetra::Interface::Jacobian> iJac = nox_interface;
  RCP<NOX::Epetra::LinearSystemAztecOO> linsys;
  RCP<Epetra_Operator> M = sg_model->create_WPrec()->PrecOp;
  RCP<NOX::Epetra::Interface::Preconditioner> iPrec = nox_interface;
  ParameterList& dirParams = noxParams->sublist("Direction");
  ParameterList& newtonParams = dirParams.sublist("Newton");
  ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Preconditioner", "User Defined");
  linsys = rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                                                    iJac, A, iPrec, M, *u));

  // Build NOX group
  RCP<NOX::Epetra::Group> grp =
    rcp(new NOX::Epetra::Group(printParams, iReq, *u, linsys));

  // Create the Solver convergence test
  ParameterList& statusParams = noxParams->sublist("Status Tests");
  RCP<NOX::StatusTest::Generic> statusTests =
    NOX::StatusTest::buildStatusTests(statusParams, utils);

  // Create the solver
  RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grp, statusTests, noxParams);

  return solver;
}

const Epetra_Vector& get_final_solution(const NOX::Solver::Generic& solver) {
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution =
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
  return finalSolution;
}

int main(int argc, char *argv[]) {
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::ParameterList;
  using Stokhos::OneDOrthogPolyBasis;
  using Stokhos::LegendreBasis;
  using Stokhos::CompletePolynomialBasis;
  using Stokhos::OrthogPolyExpansion;
  using Stokhos::AlgebraicOrthogPolyExpansion;
  using Stokhos::ParallelData;
  using Stokhos::SGModelEvaluator;
  using Stokhos::EpetraVectorOrthogPoly;

  // Initialize MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int MyPID;
  bool success = true;
  try {

    // Create a communicator for Epetra objects
    RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Create Stochastic Galerkin basis and expansion
    const int p = 5;
    Array< RCP<const OneDOrthogPolyBasis<int,double> > > bases(1);
    bases[0] = rcp(new LegendreBasis<int,double>(p));
    RCP<const CompletePolynomialBasis<int,double> > basis =
      rcp(new CompletePolynomialBasis<int,double>(bases));
    RCP<Stokhos::Sparse3Tensor<int,double> > Cijk =
      basis->computeTripleProductTensor();
    RCP<OrthogPolyExpansion<int,double> > expansion =
      rcp(new AlgebraicOrthogPolyExpansion<int,double>(basis, Cijk));

    // Create stochastic parallel distribution
    int num_spatial_procs = -1;
    ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    RCP<ParallelData> sg_parallel_data =
      rcp(new ParallelData(basis, Cijk, globalComm, parallelParams));
    RCP<const Epetra_Comm> app_comm = sg_parallel_data->getSpatialComm();

    // Create application model evaluator
    RCP<EpetraExt::ModelEvaluator> model = rcp(new SimpleME(app_comm));

    // Setup stochastic Galerkin algorithmic parameters
    RCP<ParameterList> sgParams = rcp(new ParameterList);
    ParameterList& sgPrecParams = sgParams->sublist("SG Preconditioner");
    sgPrecParams.set("Preconditioner Method", "Mean-based");
    //sgPrecParams.set("Preconditioner Method", "Approximate Gauss-Seidel");
    //sgPrecParams.set("Preconditioner Method", "Approximate Schur Complement");
    sgPrecParams.set("Mean Preconditioner Type", "Ifpack");

    // Create stochastic Galerkin model evaluator
    RCP<SGModelEvaluator> sg_model =
      rcp(new SGModelEvaluator(model, basis, Teuchos::null, expansion,
                               sg_parallel_data, sgParams));

    // Stochastic Galerkin initial guess
    // Set the mean to the deterministic initial guess, higher-order terms
    // to zero
    RCP<EpetraVectorOrthogPoly> x_init_sg = sg_model->create_x_sg();
    x_init_sg->init(0.0);
    (*x_init_sg)[0] = *(model->get_x_init());
    sg_model->set_x_sg_init(*x_init_sg);

    // Stochastic Galerkin parameters
    // Linear expansion with the mean given by the deterministic initial
    // parameter values, linear terms equal to 1, and higher order terms
    // equal to zero.
    RCP<EpetraVectorOrthogPoly> p_init_sg = sg_model->create_p_sg(0);
    p_init_sg->init(0.0);
    (*p_init_sg)[0] = *(model->get_p_init(0));
    for (int i=0; i<model->get_p_map(0)->NumMyElements(); i++)
      (*p_init_sg)[i+1][i] = 1.0;
    sg_model->set_p_sg_init(0, *p_init_sg);
    std::cout << "Stochatic Galerkin parameter expansion = " << std::endl
              << *p_init_sg << std::endl;

    // Build nonlinear solver (implemented above)
    RCP<NOX::Solver::Generic> solver = create_nox_solver(MyPID, sg_model);

    // Solve the system
    NOX::StatusTest::StatusType status = solver->solve();

    // Get final solution
    const Epetra_Vector& finalSolution = get_final_solution(*solver);

    // Convert block Epetra_Vector to orthogonal polynomial of Epetra_Vector's
    RCP<Stokhos::EpetraVectorOrthogPoly> x_sg =
      sg_model->create_x_sg(View, &finalSolution);

    if (MyPID == 0)
      std::cout << "Final Solution = " << std::endl;
    std::cout << *x_sg << std::endl;

    if (status != NOX::StatusTest::Converged)
      success = false;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if (success && MyPID == 0)
    std::cout << "Example Passed!" << std::endl;

  if (!success)
    return 1;
  return 0;
}
