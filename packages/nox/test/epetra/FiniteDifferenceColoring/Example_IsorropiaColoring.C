// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Finite Element Test Problem - Requires EpetraExt
/* Solves the nonlinear equation:
 *
 * d2u
 * --- - k * u    = 0   on      x=[0.0 to 0.5)
 * dx2
 *
 * d2u
 * --- - k * u**2 = 0   on      x=[0.5 to 1.0]
 * dx2
 *
 * subject to @ x=0, u=1
 */

#include "NOX_Common.H" // needed to determine HAVE_NOX_EPETRAEXT
#include "NOX_Epetra_FiniteDifferenceColoringWithUpdate.H"

#include <fstream>

// Ensure we have the required EpetraExt library built

#ifndef HAVE_NOX_EPETRAEXT

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{
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

  // Summarize test results
  if( Comm.MyPID() == 0 )
    std::cout << "Test failed! This test requires the EpetraExt library." << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return EXIT_FAILURE;
}

#else  // The real test

#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

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
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_IntVector.h"
#include "AztecOO.h"
#include "Teuchos_StandardCatchMacros.hpp"

#include <vector>

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "FiniteElementProblem.H"

#include "Isorropia_EpetraColorer.hpp"

using namespace std;


/*************************************/
Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> coloring_using_fdc(Teuchos::ParameterList & printParams, Teuchos::ParameterList & lsParams,Teuchos::RCP<Problem_Interface> &interface, NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<Epetra_MapColoring>&colorMap);
/**************************************/
Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> coloring_using_fdcwu(Teuchos::ParameterList & printParams, Teuchos::ParameterList & lsParams,Teuchos::RCP<Problem_Interface> &interface, NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<Epetra_MapColoring>&colorMap,Teuchos::RCP<Epetra_MapColoring>&updateColorMap
    );
/**************************************/
int solve_system(char * fpref,Epetra_Comm & Comm,Teuchos::ParameterList & printParams, Teuchos::RCP<Teuchos::ParameterList>&nlParamsPtr,Teuchos::RCP<NOX::Epetra::Interface::Required> &iReq,NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys,bool verbose);
/**************************************/
Teuchos::RCP<Epetra_CrsGraph> ExtractSubGraph(Teuchos::RCP<Epetra_CrsGraph> graph, int Nlid,const int* lids);
/**************************************/

int main(int argc, char *argv[])
{
  int ierr = 0;

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

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Get the number of elements from the command line
  int NumGlobalElements = 0;
  if ((argc > 2) && (verbose))
    NumGlobalElements = atoi(argv[2]) + 1;
  else if ((argc > 1) && (!verbose))
    NumGlobalElements = atoi(argv[1]) + 1;
  else
    NumGlobalElements = 101;

  bool success = false;

  try {
    // The number of unknowns must be at least equal to the
    // number of processors.
    if (NumGlobalElements < NumProc) {
      std::cout << "Error: numGlobalBlocks = " << NumGlobalElements
        << " cannot be < number of processors = " << NumProc << std::endl;
      throw;
    }

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the
    // function (RHS) and Jacobian evaluation routines.
    FiniteElementProblem Problem(NumGlobalElements, Comm);

    // Get the vector from the Problem
    Teuchos::RCP<Epetra_Vector> soln = Problem.getSolution();
    NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

    // Initialize Solution
    soln->PutScalar(1.0);

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
    printParams.set("MyPID", MyPID);
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
    if (verbose)
      printParams.set("Output Information",
          NOX::Utils::OuterIteration +
          NOX::Utils::OuterIterationStatusTest +
          NOX::Utils::LinearSolverDetails +
          NOX::Utils::InnerIteration +
          NOX::Utils::Parameters +
          NOX::Utils::Details +
          NOX::Utils::Warning);
    else
      printParams.set("Output Information", NOX::Utils::Error +
          NOX::Utils::TestDetails);


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
    lsParams.set("Output Frequency", 1);
    lsParams.set("Preconditioner", "Ifpack");
    lsParams.set("Max Age Of Prec", 5);

    // Create the interface between the test problem and the nonlinear solver
    // This is created by the user using inheritance of the abstract base class:
    // NLS_PetraGroupInterface
    Teuchos::RCP<Problem_Interface> interface =
      Teuchos::rcp(new Problem_Interface(Problem));

    const Epetra_CrsGraph & G=*Problem.getGraph();

    // Create the Epetra_RowMatrix using Finite Difference with Coloring
    Teuchos::ParameterList isorParamList;
    // Teuchos::ParameterList& zoltanParamList = isorParamList.sublist("ZOLTAN");
    // zoltanParamList.set("DISTANCE","2");
    Isorropia::Epetra::Colorer isorColorer(
        (Teuchos::RCP<const Epetra_CrsGraph>) Problem.getGraph(), isorParamList, false);

    Teuchos::RCP<Epetra_MapColoring> colorMap =  isorColorer.generateColMapColoring();

    // Build the update coloring
    Teuchos::RCP<Epetra_CrsGraph> subgraph=ExtractSubGraph(Problem.getGraph(),Problem.NumLocalNonLinearUnknowns(),Problem.getNonLinearUnknowns());
    Isorropia::Epetra::Colorer updateIsorColorer((Teuchos::RCP<const Epetra_CrsGraph>)subgraph, isorParamList, false);

    //  Teuchos::RCP<Epetra_MapColoring> updateColorMap_limited=updateIsorColorer.generateColMapColoring();
    Teuchos::RCP<Epetra_MapColoring> updateColorMap=updateIsorColorer.generateColMapColoring();

    // Explictly recolor dummies to color zero
    int Nlid=Problem.NumLocalNonLinearUnknowns();
    const int *lids=Problem.getNonLinearUnknowns();
    Epetra_IntVector rIDX(G.RowMap());
    for(int i=0;i<Nlid;i++)
      rIDX[lids[i]]=1;
    Epetra_IntVector *cIDX;
    if(G.RowMap().SameAs(G.ColMap()))
      cIDX=&rIDX;
    else{
      cIDX=new Epetra_IntVector(G.ColMap(),true);
      cIDX->Import(rIDX,*G.Importer(),Insert);
    }

    for(int i=0;i<G.NumMyCols();i++)
      if((*cIDX)[i]==0) (*updateColorMap)[i]=0;

    if(!G.RowMap().SameAs(G.ColMap()))
      delete cIDX;

    int base_colors=colorMap->MaxNumColors(),update_colors=updateColorMap->MaxNumColors();
    if(!MyPID)
      std::cout<<"First time colors = "<<base_colors<<" Update colors = "<<update_colors-1<<endl;


    // Use this constructor to create the graph numerically as a means of timing
    // the old way of looping without colors :
    //  NOX::Epetra::FiniteDifferenceColoring A(interface, soln,
    //                                          *colorMap, *columns);
    // Or use this as the standard way of using finite differencing with coloring
    // where the application is responsible for creating the matrix graph
    // beforehand, ie as is done in Problem.
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;

    NOX::Epetra::Vector noxSoln2(noxSoln);

    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
      coloring_using_fdc(printParams,lsParams,interface,noxSoln,Problem,colorMap);
    char name_fdc[] = "fdc";
    ierr=solve_system(name_fdc,Comm,printParams,nlParamsPtr,iReq,noxSoln,Problem,linSys,verbose);

    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys2 =
      coloring_using_fdcwu(printParams,lsParams,interface,noxSoln2,Problem,colorMap,updateColorMap);
    char name_fdcwu[] = "fdcwu"; 
    ierr=solve_system(name_fdcwu,Comm,printParams,nlParamsPtr,iReq,noxSoln2,Problem,linSys2,verbose);

    success = ierr==0;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}



/*************************************/
Teuchos::RCP<Epetra_CrsGraph> ExtractSubGraph(Teuchos::RCP<Epetra_CrsGraph> graph, int Nlid,const int* lids){
  Epetra_CrsGraph &G=*graph;
  int N=G.NumMyRows();
  int *indices,NumIndices;

  // Build index vectors for rows
  Epetra_IntVector rIDX(G.RowMap());
  for(int i=0;i<Nlid;i++)
    rIDX[lids[i]]=1;

  // Import/export to get column vector if needed
  Epetra_IntVector *cIDX;
  if(G.RowMap().SameAs(G.ColMap()))
    cIDX=&rIDX;
  else{
    cIDX=new Epetra_IntVector(G.ColMap(),true);
    cIDX->Import(rIDX,*G.Importer(),Insert);
  }

  // Allocate new matrix
  Epetra_CrsGraph *Sub=new Epetra_CrsGraph(Copy,G.RowMap(),G.ColMap(),0,false);

  // Do extraction
  for(int i=0;i<N;i++){
    int rgid=G.GRID(i);
    if(rIDX[i]!=1) continue;
    G.ExtractMyRowView(i,NumIndices,indices);
    for(int k=0;k<NumIndices;k++){
      if((*cIDX)[indices[k]]==1){
        int cgid=G.GCID(indices[k]);
        Sub->InsertGlobalIndices(rgid,1,&cgid);
      }
    }
  }
  Sub->FillComplete();
  Sub->OptimizeStorage();

  if(!G.RowMap().SameAs(G.ColMap()))
    delete cIDX;

  return Teuchos::RCP<Epetra_CrsGraph>(Sub);
}

/*************************************/
Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> coloring_using_fdc(Teuchos::ParameterList & printParams, Teuchos::ParameterList & lsParams,Teuchos::RCP<Problem_Interface> &interface, NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<Epetra_MapColoring>&colorMap){

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;

  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  Teuchos::RCP< std::vector<Epetra_IntVector> > columns =
    Teuchos::rcp(&colorMapIndex(*(Problem.getGraph())));

  Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoring> A =
    Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams,
          iReq, noxSoln,
          Problem.getGraph(),
          colorMap,
          columns,
          true));
  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
          iReq, iJac, A, noxSoln));
  return linSys;
}

/**************************************/
Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> coloring_using_fdcwu(Teuchos::ParameterList & printParams, Teuchos::ParameterList & lsParams,Teuchos::RCP<Problem_Interface> &interface, NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<Epetra_MapColoring>&colorMap,Teuchos::RCP<Epetra_MapColoring>&updateColorMap){

  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  // Use constructor w/ update
  Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoringWithUpdate> A =
    Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoringWithUpdate(printParams,
          iReq, noxSoln,
          Problem.getGraph(),
          colorMap,updateColorMap));


  // Create the linear system
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
          iReq, iJac, A, noxSoln));
  return linSys;
}



/*************************************/
int solve_system(char *fpref,Epetra_Comm & Comm,Teuchos::ParameterList & printParams, Teuchos::RCP<Teuchos::ParameterList>&nlParamsPtr,Teuchos::RCP<NOX::Epetra::Interface::Required> &iReq,NOX::Epetra::Vector &noxSoln,FiniteElementProblem& Problem,Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys,bool verbose){
  int MyPID=Comm.MyPID();
  int i;

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr =
    Teuchos::rcp(new NOX::Epetra::Group(printParams,
          iReq,
          noxSoln,
          linSys));
  NOX::Epetra::Group& grp = *(grpPtr.get());

  // ATOL vector if using NOX::StatusTest::WRMS
  Teuchos::RCP<Epetra_Vector> soln = Problem.getSolution();
  NOX::Epetra::Vector weights(soln);
  weights.scale(1.0e-8);

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormF> relresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
  Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create the method
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = solver->solve();

  if (verbose) {
    if (status != NOX::StatusTest::Converged)
      if (MyPID==0)
        cout << "Nonlinear solver failed to converge!" << std::endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (verbose) {
    NOX::Utils utils(printParams);
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      std::cout << std::endl << "Final Parameters" << std::endl
        << "****************" << std::endl;
      solver->getList().print(cout);
      std::cout << std::endl;
    }
  }

  // Print solution
  char file_name[80];
  FILE *ifp;
  int NumMyElements = soln->Map().NumMyElements();
  (void) sprintf(file_name, "%s.output.%d",fpref,MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %22.16e\n", soln->Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);

  // Summarize test results
  if (Comm.MyPID() == 0) {
    if (status == NOX::StatusTest::Converged)
      std::cout << "Test passed!" << std::endl;
    else
      std::cout << "Test failed!" << std::endl;
  }
  return 0;
}
#endif
