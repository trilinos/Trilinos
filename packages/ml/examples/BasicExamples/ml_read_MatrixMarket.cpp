
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Teuchos_CommandLineProcessor.hpp"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_MultiVectorIn.h"
#include "AztecOO.h"

// includes required by ML
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"
#include <fstream>

using namespace Teuchos;

void ML_Read_Matrix_Dimensions(const char *filename, int *numGlobalRows,Epetra_Comm &Comm);
// Small function to handle exiting gracefully.  All pids must call.

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  int mypid = Comm.MyPID();
  Teuchos::CommandLineProcessor clp(false);
  std::string matrixFile = ""; clp.setOption("matrix", &matrixFile, "Matrix market file containing matrix.  [REQUIRED]");
  std::string xmlFile = ""; clp.setOption("xml", &xmlFile, "XML file containing ML options. [OPTIONAL]");
  std::string nullspaceFile = ""; clp.setOption("nullspace", &nullspaceFile, "File containing nullspace modes. [OPTIONAL]");
  std::string coordFile = ""; clp.setOption("coord", &coordFile, "File containing coordinate vectors. [OPTIONAL]");
  std::string rhsFile = ""; clp.setOption("rhs", &rhsFile, "File containing right-hand side vector.  [OPTIONAL]");
  std::string mapFile = ""; clp.setOption("map", &mapFile, "File containing matrix rowmap.  [OPTIONAL]");
  int numPDEs = 1; clp.setOption("npdes", &numPDEs, "Number of PDEs. [Default=1]");
  std::string krylovSolver = "gmres"; clp.setOption("krylov", &krylovSolver, "outer Krylov solver.");
  int output=10; clp.setOption("output", &output, "how often to print residual history.");
  int maxits=100; clp.setOption("maxits", &maxits, "maximum number of Krylov iterations.  [OPTIONAL]");
  double tol=1e-8; clp.setOption("tol", &tol, "Krylov solver tolerance.  [OPTIONAL]");
  //bool   printTimings     = true;              clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  int indexBase = 1;
  Epetra_Map *RowMap=NULL;
  int numGlobalRows = -999;
  if (mapFile != "") {
    if (!mypid) std::cout << "reading rowmap from " << mapFile << std::endl;
    EpetraExt::MatrixMarketFileToMap(mapFile.c_str(),Comm,RowMap);
  } else {
    if (Comm.NumProc() > 1) {
      // If parallel and the rowmap hasn't been given explicitly,
      // get the matrix dimensions and create a row map that
      // will not break aggregation procedure.  (On a processor, the
      // number of local rows must be divisible by #dof per node.)
      // The main idea is that the dof's associated with a node should all
      // reside on the same processor.
      if (!mypid) std::cout << "creating compatible rowmap" << std::endl;
      ML_Read_Matrix_Dimensions(matrixFile.c_str(), &numGlobalRows, Comm);
      int numNodes = numGlobalRows / numPDEs;
      if ((numGlobalRows - numNodes * numPDEs) != 0 && !mypid)
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Number of matrix rows is not divisible by #dofs");
      int numMyNodes;
      int nproc = Comm.NumProc();
      if (Comm.MyPID() < nproc-1) numMyNodes = numNodes / nproc;
      else numMyNodes = numNodes - (numNodes/nproc) * (nproc-1);
      RowMap = new Epetra_Map(numGlobalRows,numMyNodes*numPDEs,indexBase,Comm);
    }
  }


  int errcode=0;
  Epetra_CrsMatrix *A;
  Epetra_MultiVector *nullspaceVector=0;
  Epetra_MultiVector *pRHS=0;

  if (!mypid) std::cout << "reading matrix from " << matrixFile << std::endl;
  if (RowMap) errcode =EpetraExt::MatrixMarketFileToCrsMatrix(matrixFile.c_str(), *RowMap, A);
  else        errcode =EpetraExt::MatrixMarketFileToCrsMatrix(matrixFile.c_str(), Comm, A);
  TEUCHOS_TEST_FOR_EXCEPTION(errcode, std::runtime_error, "error reading file " + matrixFile);

  errcode=0;
  if (rhsFile != "") {
    if (!mypid) std::cout << "reading rhs vector from " << rhsFile << std::endl;
    errcode =EpetraExt::MatrixMarketFileToMultiVector(rhsFile.c_str(), A->DomainMap(),pRHS);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(errcode, std::runtime_error, "error reading file " + rhsFile);

  errcode=0;
  if (nullspaceFile != "") {
    if (!mypid) std::cout << "reading rigid body modes from " << nullspaceFile << std::endl;
    errcode =EpetraExt::MatrixMarketFileToMultiVector(nullspaceFile.c_str(), A->DomainMap(),nullspaceVector);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(errcode, std::runtime_error, "error reading file " + nullspaceFile);

  Epetra_MultiVector *coordVector=0;
  Epetra_Map coordMap(A->NumGlobalRows()/numPDEs,indexBase,Comm);
  errcode=0;
  if (coordFile != "") {
    if (!mypid) std::cout << "reading coordinates from " << coordFile;
    errcode =EpetraExt::MatrixMarketFileToMultiVector(coordFile.c_str(), coordMap, coordVector);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(errcode, std::runtime_error, "error reading file " + coordFile);

  // ML expects the nullspace modes in a single double vector.

  double *nullspacePointer=0;
  if (nullspaceVector) {
    int MyLength = nullspaceVector->MyLength();
    nullspacePointer = new double[nullspaceVector->NumVectors() * MyLength];
    for (int k=0; k < nullspaceVector->NumVectors(); k++)
      for (int j=0; j < MyLength; j++)
         nullspacePointer[k*MyLength + j] = (*nullspaceVector)[k][j];
  }

  // ML expects coordinates in separate double vectors.

  double* mv=0;
  int stride;
  errcode=0;
  if (coordVector) {
    errcode = coordVector->ExtractView(&mv,&stride);
  }
  TEUCHOS_TEST_FOR_EXCEPTION(errcode,std::runtime_error,"error extracting pointers to coordinates");

  double **xyz=0;
  if (coordVector) {
    xyz = new double*[coordVector->NumVectors()];
    for (int i=0; i<coordVector->NumVectors(); i++) xyz[i] = mv+i*stride;
  }

  // =========================== begin of ML part ===========================

  // create a parameter list for ML options
  ParameterList MLList;

  ML_Epetra::SetDefaults("SA", MLList);
  //MLList.set("ML output",10);
  MLList.set("PDE equations",numPDEs);
  //MLList.set("smoother: type","symmetric Gauss-Seidel");
  if (coordVector) {
    if (xyz[0]) MLList.set("x-coordinates",xyz[0]);
    if (xyz[1]) MLList.set("y-coordinates",xyz[1]);
    if (coordVector->NumVectors() == 3)
      MLList.set("z-coordinates",xyz[2]);
    MLList.set("null space: type","from coordinates");
  }
  if (nullspacePointer) {
    MLList.set("null space: type","pre-computed");
    MLList.set("null space: dimension",nullspaceVector->NumVectors());
    MLList.set("null space: vectors",nullspacePointer);
  }

  // Read in XML options
  if (xmlFile != "")
    ML_Epetra::ReadXML(xmlFile,MLList,Comm);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================

  const Epetra_Map Map = A->DomainMap();
  bool populateRhs=false;
  if (pRHS==0) {
    pRHS = new Epetra_Vector(Map);
    populateRhs=true;
  }
  Epetra_MultiVector RHS = *pRHS;

  //Epetra_Vector LHS(Copy,RHS);
  Epetra_Vector LHS(Map);
  LHS.PutScalar(0.0);
  double mynorm;
  Epetra_Vector trueX(Map);
  if (populateRhs) {
    trueX.SetSeed(846930886);
    trueX.Random();
    trueX.Norm2(&mynorm);
    //trueX.Scale(1.0/mynorm);
    A->Multiply(false,trueX,RHS);
  }

  Epetra_LinearProblem Problem(A,&LHS,&RHS);
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_scaling, AZ_none);
  if (krylovSolver=="gmres")
    solver.SetAztecOption(AZ_solver, AZ_gmres);
  else if (krylovSolver=="cg")
    solver.SetAztecOption(AZ_solver, AZ_cg);
  else if (krylovSolver=="fixedpt")
    solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Bad Krylov solver option");
  //solver.SetAztecOption(AZ_kspace, 100);
  solver.SetAztecOption(AZ_output, output);
  //solver.SetAztecOption(AZ_conv, AZ_noscaled);
  solver.SetPrecOperator(MLPrec);

  solver.Iterate(maxits, tol);

  //Calculate a final residual
  Epetra_Vector workvec(Map);
  A->Multiply(false,LHS,workvec);
  workvec.Update(1.0,RHS,-1.0);
  RHS.Norm2(&mynorm);
  workvec.Scale(1./mynorm);
  workvec.Norm2(&mynorm);
  if (Comm.MyPID() == 0) std::cout << "||r||_2 = " << mynorm << std::endl;
  //Calculate a relative error

  // delete the preconditioner. Do it BEFORE calling MPI_Finalize
  delete MLPrec;
  delete [] nullspacePointer;
  delete nullspaceVector;
  delete pRHS;
  if (coordVector) delete coordVector;
  if (xyz) delete [] xyz;

  delete A;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);

} //main

//** ************************************************************************ **

void ML_Read_Matrix_Dimensions(const char *filename, int *numGlobalRows, Epetra_Comm &Comm)
{
    const int lineLength = 1025;
    char line[lineLength+1], token1[50], token2[50], token3[50], token4[50], token5[50];
    FILE *fid = fopen(filename,"r");
    int N, NZ;
    TEUCHOS_TEST_FOR_EXCEPTION(fgets(line, lineLength, fid)==0,std::runtime_error,"error opening matrix file");
    TEUCHOS_TEST_FOR_EXCEPTION(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0,
                                      std::runtime_error,"error opening matrix file");
    TEUCHOS_TEST_FOR_EXCEPTION((strcmp(token1, "%%MatrixMarket") || strcmp(token2, "matrix") ||
                               strcmp(token3, "coordinate") || strcmp(token4, "real") ||
                               strcmp(token5, "general")),
                               std::runtime_error,"error opening matrix file");
    // Next, strip off header lines (which start with "%")
    do {
      TEUCHOS_TEST_FOR_EXCEPTION(fgets(line, lineLength, fid)==0,
                                 std::runtime_error,"error opening matrix file");
    } while (line[0] == '%');

    // Next get problem dimensions: M, N, NZ
      TEUCHOS_TEST_FOR_EXCEPTION(sscanf(line, "%d %d %d", numGlobalRows, &N, &NZ)==0,
                                 std::runtime_error,"error opening matrix file");
} //ML_Read_Matrix_Dimensions()
