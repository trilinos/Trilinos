// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Finite Element PelletTransport Test Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 *
 * This problem is examined with a variety of time integration schemes in:
 * "Studies on the Convergence of Various Time-Integration Schemes for the
 * Radiation-Diffusion Problem," Curtis C. Ober & John N. Shadid, in prep.
 *
 * In this example, only a 1st-order fully implicit (backward Euler)
 * time integration scheme is considered currently.
 *
 * Values for time step size and finite spatial extent are specified in
 * the constructor initialization list in PelletTransport.C using
 * variables dt  and xmin,xmax, respectively.
 * The number of time steps to be taken is specified by variable
 * maxTimeSteps below.
 */

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"
#include "../test/utils/NOX_Epetra_DebugTools.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// For parsing command line
#include "Teuchos_CommandLineProcessor.hpp"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring
#include <vector>
#ifdef HAVE_NOX_EPETRAEXT       // Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#endif

// User's application specific files
#include "MaterialProps.H" // Interface file to NOX
#include "Problem_Interface.H" // Interface file to NOX
#include "PelletTransport.H"

using namespace std;

int main(int argc, char *argv[])
{
  std::cout << log(10.0) << std::endl;

  // Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::CommandLineProcessor clp( false );

  // Default run-time options that can be changed from the command line
  bool   verbose                = true    ;
  bool   restart                = false   ;
  int    NumGlobalElementsUO2   = 2001    ;
  int    NumGlobalElementsHe    = 2001    ;
  int    NumGlobalElementsClad  = 2001    ;
  // Physical parameters
  double qdot                   = 2.0e8   ;
  double xB                     = 0.02    ;
  double xminUO2                = 0.0     ;
  double xmaxUO2                = 0.0043  ;
  double xminHe                 = 0.0043  ;
  double xmaxHe                 = 0.00433 ;
  double xminClad               = 0.00433 ;
  double xmaxClad               = 0.00483 ;


  clp.setOption( "verbose", "no-verbose", &verbose, "Verbosity on or off." );
  clp.setOption( "restart", "no-restart", &restart, "Restart from previous solution." );
  clp.setOption( "nUO2",  &NumGlobalElementsUO2, "Number of elements in UO2 domain" );
  clp.setOption( "nHe",   &NumGlobalElementsHe, "Number of elements in He domain" );
  clp.setOption( "nClad", &NumGlobalElementsClad, "Number of elements in Cladding domain" );
  clp.setOption( "qdot", &qdot, "Source term coefficient, Qdot" );
  clp.setOption( "xB", &xB, "Non-stoichiometry at outer surface of fuel rod, xB" );
  clp.setOption( "xminUO2", &xminUO2, "Left boundary location of UO2 domain (m)" );
  clp.setOption( "xmaxUO2", &xmaxUO2, "Right boundary location of UO2 domain (m)" );
  clp.setOption( "xminHe", &xminHe, "Left boundary location of He domain (m)" );
  clp.setOption( "xmaxHe", &xmaxHe, "Right boundary location of He domain (m)" );
  clp.setOption( "xminClad", &xminClad, "Left boundary location of Cladding domain (m)" );
  clp.setOption( "xmaxClad", &xmaxClad, "Right boundary location of Cladding domain (m)" );

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL )
    return parse_return;

  // Create and reset the Timer
  Epetra_Time myTimer(Comm);
  double startWallTime = myTimer.WallTime();

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // The number of unknowns must be at least equal to the
  // number of processors.
  //if (NumGlobalNodes < NumProc)
  //{
  //  std::cout << "numGlobalNodes = " << NumGlobalNodes << " cannot be < number of processors = " << NumProc << std::endl;
  //  std::cout << "Test failed!" << std::endl;
  //  exit(1);
  //}

  // Create the PelletTransport problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the
  // function (F) and Jacobian evaluation routines.
  Teuchos::RCP<PelletTransport> Problem =
    Teuchos::rcp(new PelletTransport(
                   NumGlobalElementsUO2  , xminUO2  , xmaxUO2  ,
                   NumGlobalElementsHe   , xminHe   , xmaxHe   ,
                   NumGlobalElementsClad , xminClad , xmaxClad ,
                   Comm, restart) );

  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = Problem->getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  // Plug in the needed material models
  MaterialPropBase * propEnity;
  propEnity = MaterialPropFactory::factory().add_model( new MaterialProp_He()   );
  propEnity = MaterialPropFactory::factory().add_model( new MaterialProp_Clad() );
  propEnity = MaterialPropFactory::factory().add_model( new MaterialProp_UO2()  );
  dynamic_cast<MaterialProp_UO2 *>(propEnity)->set_qdot( qdot );
  dynamic_cast<MaterialProp_UO2 *>(propEnity)->set_xB( xB );



  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  if (verbose) {
    printParams.set("MyPID", MyPID);
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    printParams.set("Output Information",
                 NOX::Utils::OuterIteration +
                 NOX::Utils::OuterIterationStatusTest +
                 NOX::Utils::InnerIteration +
                 NOX::Utils::LinearSolverDetails +
                 NOX::Utils::Parameters +
                 NOX::Utils::Details +
                 NOX::Utils::Warning);
  }
  else
    printParams.set("Output Information", NOX::Utils::Error);

  // Create a print handle object
  NOX::Utils utils(printParams);

  // Sublist for line search
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "Full Step");

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 800);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 100);
  lsParams.set("Preconditioner", "AztecOO");

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RCP<Problem_Interface> interface =
    Teuchos::rcp(new Problem_Interface(*Problem));

#ifndef HAVE_NOX_EPETRAEXT
  utils.out() << "Cannot use Coloring without package epetraext !!!!" << std::endl;
  utils.out() << "Test passed!" << std::endl;
  exit(0);
#else
  // Create a timer for performance
  Epetra_Time fillTime(Comm);

  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType =
    EpetraExt::CrsGraph_MapColoring::JONES_PLASSMAN;
  bool colorParallel = true;
  int reordering = 0;
  int verbosityLevel = 0;
  bool distance1 = false;
  EpetraExt::CrsGraph_MapColoring tmpMapColoring(
    algType, reordering, distance1, verbosityLevel);
  Teuchos::RCP<Epetra_MapColoring> colorMap =
    Teuchos::rcp(&tmpMapColoring(*(Problem->getGraph())));
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  Teuchos::RCP< std::vector<Epetra_IntVector> > columns =
    Teuchos::rcp(&colorMapIndex(*(Problem->getGraph())));

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoring> FDC =
    Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams,
                               iReq,
                               noxSoln,
                               Problem->getGraph(),
                               colorMap,
                               columns,
                               colorParallel,
                               distance1));

  if (verbose)
    printf("\n[%d]\tTime to color Jacobian --> %e sec.  for %d colors.\n\n",
       MyPID,fillTime.ElapsedTime(), colorMap->NumColors());

  FDC->setDifferenceMethod(NOX::Epetra::FiniteDifference::Centered);

// -------------- End of block needed to use coloring --------------- */


  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = FDC;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                              iReq, iJac, FDC,
                              noxSoln));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
                    iReq,
                    noxSoln,
                    linSys));
  NOX::Epetra::Group& grp = *(grpPtr.get());

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-5, NOX::StatusTest::NormF::Unscaled));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(25));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> finiteval =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue());
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(absresid);
  combo->addStatusTest(maxiters);
  combo->addStatusTest(finiteval);

  // Create the method
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

  // Initialize time integration parameters
  int maxTimeSteps = 1;
  int timeStep = 0;
  double time = 100.;
  double dt = Problem->getdt();

  // Print initial solution
  char file_name[25];
  FILE *ifp;
  Epetra_Vector& xMesh = *(Problem->getMesh());
  int NumMyNodes = xMesh.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d_%d",MyPID,timeStep);
  ifp = fopen(file_name, "w");
  for( int i = 0; i < NumMyNodes; ++i )
    fprintf(ifp, "%d  %E  %E  %E\n", xMesh.Map().MinMyGID()+i,
        xMesh[i], (*soln)[2*i], (*soln)[2*i+1]);
  fclose(ifp);

  NOX::StatusTest::StatusType status = NOX::StatusTest::Unconverged;

  // Time integration loop
  while(timeStep < maxTimeSteps) {

    timeStep++;
    time += dt;

    if (verbose)
      utils.out() << "Time Step: " << timeStep << ",\tTime: " << time << std::endl;

    status = NOX::StatusTest::Unconverged;
    status = solver->solve();

    if (verbose)
      if (status != NOX::StatusTest::Converged)
    if (MyPID==0)
      utils.out() << "Nonlinear solver failed to converge!" << std::endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup =
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

    // Get the Epetra_Vector with the last residual from the solver
    const Epetra_Vector& finalResidual =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

    // End Nonlinear Solver **************************************

    if( 1 )
    {
      // Print solution
      (void) sprintf(file_name, "T_output.%03d_%05d",MyPID,timeStep);
      ifp = fopen(file_name, "w");
      for( int i = 0; i < NumMyNodes; ++i )
        fprintf(ifp, "%d  %E  %E\n", soln->Map().MinMyGID()+i,
                                     xMesh[i], finalSolution[2*i]);

      fclose(ifp);

      (void) sprintf(file_name, "x_output.%03d_%05d",MyPID,timeStep);
      ifp = fopen(file_name, "w");
      for( int i = 0; i < NumGlobalElementsUO2 + 1; ++i )
        fprintf(ifp, "%d  %E  %E\n", soln->Map().MinMyGID()+i,
                                     xMesh[i], finalSolution[2*i+1]);
      fclose(ifp);

      // Print residual
      (void) sprintf(file_name, "residual.%03d_%05d",MyPID,timeStep);
      ifp = fopen(file_name, "w");
      for( int i = 0; i < NumMyNodes; ++i )
        fprintf(ifp, "%d  %E  %E  %E\n", soln->Map().MinMyGID()+i,
                                     xMesh[i], finalResidual[2*i],
                                     finalResidual[2*i+1]);
      fclose(ifp);

      NOX::Epetra::DebugTools::writeVector( "restartVec", finalSolution, NOX::Epetra::DebugTools::MATRIX_MARKET );
    }

    utils.out() << "Time Step: " << timeStep << ",\tTime: " << time << ", Tmax = " << finalSolution[0] << ", Xmax = " << finalSolution[1] << std::endl;

    Problem->reset(finalSolution);
    grp.setX(finalSolution);
    solver->reset(grp.getX(), combo);
    grp.computeF();

  } // end time step while loop

  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(grp.getX())).getEpetraVector();
  NOX::Epetra::DebugTools::writeVector( "transient_restartVec", finalSolution, NOX::Epetra::DebugTools::MATRIX_MARKET );

  // Output the parameter list
  if (verbose) {
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << std::endl << "Final Parameters" << std::endl
       << "****************" << std::endl;
      solver->getList().print(utils.out());
      utils.out() << std::endl;
    }
  }

  // Output timing info
  if (verbose)
  {
    if(MyPID==0)
      utils.out() << "\nTimings :\n\tWallTime --> " << myTimer.WallTime() - startWallTime << " sec."
       << "\n\tElapsedTime --> " << myTimer.ElapsedTime()
       << " sec." << std::endl << std::endl;

  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return 0;
#endif // HAVE_NOX_EPETRAEXT
}
