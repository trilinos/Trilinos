
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

//usage: ./driver.exe -f A.dat -rbm rbm.dat -i options.xml

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
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

#include "TumiStuff.h"

using namespace Teuchos;

void ML_Read_Matrix_Dimensions(char *filename, int *numGlobalRows,Epetra_Comm &Comm);
// Small function to handle exiting gracefully.  All pids must call.
void ML_Exit(int mypid, int code, char *fmt, ...);

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
  int i=1;
  int mypid = Comm.MyPID();
  char matrixFile[80] = "\0";
  char rbmFile[80] = "\0";
  char coordFile[80] = "\0";
  char xmlFile[80] = "\0";
  char rhsFile[80] = "\0";
  int numPDEs = 1;

  while (i<argc) {
    if      (strncmp(argv[i],"-f",2) == 0) {strncpy(matrixFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-rbm",4) == 0) {strncpy(rbmFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-c",2) == 0) {strncpy(coordFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-i",2) == 0) {strncpy(xmlFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-b",2) == 0) {strncpy(rhsFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-npdes",6) == 0) {numPDEs = (int) strtol(argv[i+1],NULL,10); i+=2;}
    else
      ML_Exit(mypid,EXIT_FAILURE,"Unrecognized option %s\nUsage:  ./ml_read_MatrixMarket.exe -f matrix_file [-rbm rigid_body_mode_file] [-i xml_input_file]\n",argv[i]);
  }

  int indexBase = 1;
  Epetra_Map *RowMap=NULL;
  int numGlobalRows = -999;
  if (Comm.NumProc() > 1) {
    // In parallel, get matrix dimension and create row map that
    // will not break aggregation procedure.  (On a processor, the
    // number of local rows must be divisible by #dof per node.)
    // The main idea is that the dof's associated with a node should all
    // reside on the same processor.
    ML_Read_Matrix_Dimensions(matrixFile, &numGlobalRows, Comm);
    int numNodes = numGlobalRows / numPDEs;
    if ((numGlobalRows - numNodes * numPDEs) != 0)
      ML_Exit(mypid,EXIT_FAILURE,"Number of matrix rows is not divisible by #dofs");
    int numMyNodes;
    int nproc = Comm.NumProc();
    if (Comm.MyPID() < nproc-1) numMyNodes = numNodes / nproc;
    else numMyNodes = numNodes - (numNodes/nproc) * (nproc-1);
    RowMap = new Epetra_Map(numGlobalRows,numMyNodes*numPDEs,indexBase,Comm);
  }


  int errcode=0;
  Epetra_CrsMatrix *A;
  Epetra_MultiVector *rbmVector=0;
  Epetra_MultiVector *pRHS=0;

  if (!mypid) printf("reading matrix from %s\n",matrixFile);
  if (RowMap) errcode =EpetraExt::MatrixMarketFileToCrsMatrix(matrixFile, *RowMap, A);
  else        errcode =EpetraExt::MatrixMarketFileToCrsMatrix(matrixFile, Comm, A);
  if (errcode) ML_Exit(mypid,EXIT_FAILURE,"error reading file %s", matrixFile);

  errcode=0;
  if (strncmp(rhsFile,"\0",2) != 0) {
    if (!mypid) printf("reading rhs vector from %s\n",rhsFile);
    errcode =EpetraExt::MatrixMarketFileToMultiVector(rhsFile, A->DomainMap(),pRHS);
  }
  if (errcode) ML_Exit(mypid,EXIT_FAILURE,"error reading file %s", rhsFile);

  errcode=0;
  if (strncmp(rbmFile,"\0",2) != 0) {
    if (!mypid) printf("reading rigid body modes from %s\n",rbmFile);
    errcode =EpetraExt::MatrixMarketFileToMultiVector(rbmFile, A->DomainMap(),rbmVector);
  }
  if (errcode) ML_Exit(mypid,EXIT_FAILURE,"error reading file %s", rbmFile);

  Epetra_MultiVector *coordVector=0;
  Epetra_Map coordMap(A->NumGlobalRows()/numPDEs,indexBase,Comm);
  errcode=0;
  if (strncmp(coordFile,"\0",2) != 0) {
    if (!mypid) printf("reading coordinates from %s\n",coordFile);
    errcode =EpetraExt::MatrixMarketFileToMultiVector(coordFile, coordMap, coordVector);
  }
  if (errcode) ML_Exit(mypid,EXIT_FAILURE,"error reading file %s", coordFile);

  // ML expects the rigid body modes in a single double vector.
  
  double *rbmPointer=0;
  if (rbmVector) {
    int MyLength = rbmVector->MyLength();
    rbmPointer = new double[rbmVector->NumVectors() * MyLength];
    for (int k=0; k < rbmVector->NumVectors(); k++)
      for (int j=0; j < MyLength; j++)
         rbmPointer[k*MyLength + j] = (*rbmVector)[k][j];
  }

  // ML expects coordinates in separate double vectors.
  
  double* mv=0;
  int stride;
  errcode=0;
  if (coordVector) {
    errcode = coordVector->ExtractView(&mv,&stride);
  }
  if (errcode) ML_Exit(mypid,EXIT_FAILURE,"error extracting pointers to coordinates\n");

  double **xyz=0;
  if (coordVector) {
    xyz = new double*[coordVector->NumVectors()];
    for (int i=0; i<coordVector->NumVectors(); i++) xyz[i] = mv+i*stride;
  }

  // =========================== begin of ML part ===========================

  MatrixVecChecker(*A);
  MPI_Finalize();
  return(EXIT_SUCCESS);

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
  if (rbmPointer) {
    MLList.set("null space: type","pre-computed");
    MLList.set("null space: dimension",rbmVector->NumVectors());
    MLList.set("null space: vectors",rbmPointer);
  }

  // Read in XML options
  if (strncmp(xmlFile,"\0",2) != 0)
    ML_Epetra::ReadXML(xmlFile,MLList,Comm);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================

  const Epetra_Map Map = A->DomainMap();
  if (pRHS==0) pRHS = new Epetra_Vector(Map);
  Epetra_MultiVector RHS = *pRHS;

  Epetra_MultiVector LHS(RHS);
  //Epetra_Vector LHS(Map);
  LHS.PutScalar(0.0);
  double mynorm;
  RHS.Norm2(&mynorm);
  cout << "rhs norm = " << mynorm << std::endl;
  /*
  Epetra_Vector trueX(Map);
  trueX.SetSeed(90201);
  trueX.Random();
  trueX.Norm2(&mynorm);
  trueX.Scale(1.0/mynorm);
  A->Multiply(false,trueX,RHS);
  */

  Epetra_LinearProblem Problem(A,&LHS,&RHS);
  AztecOO solver(Problem);

  solver.SetAztecOption(AZ_scaling, AZ_none);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  //solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_kspace, 50);
  solver.SetAztecOption(AZ_output, 10);
  solver.SetAztecOption(AZ_conv, AZ_noscaled);
  solver.SetPrecOperator(MLPrec);

  solver.Iterate(100, 1e-8);

  //string finalSolFile = MLList.get("solution file","finalSol.mm");
  //EpetraExt::MultiVectorToMatrixMarketFile(finalSolFile.c_str(), LHS, "Final solution");

  //Calculate a final residual
  Epetra_Vector workvec(Map);
  A->Multiply(false,LHS,workvec);
  workvec.Update(1.0,RHS,-1.0);
  RHS.Norm2(&mynorm);
  workvec.Scale(1./mynorm);
  workvec.Norm2(&mynorm);
  if (Comm.MyPID() == 0) cout << "||r||_2 = " << mynorm << endl;
  //Calculate a relative error
  /*
  workvec.Update(1.0,trueX,-1.0,LHS,0.0);
  trueX.Norm2(&mynorm);
  workvec.Scale(1./mynorm);
  workvec.Norm2(&mynorm);
  if (Comm.MyPID() == 0) cout << "||x_true - x_calc||_2 / ||x_true||  = " << mynorm << endl;
  */

  // delete the preconditioner. Do it BEFORE calling MPI_Finalize
  delete MLPrec;
  delete [] rbmPointer;
  delete rbmVector;
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

void ML_Read_Matrix_Dimensions(char *filename, int *numGlobalRows, Epetra_Comm &Comm)
{
    char line[35], token1[35], token2[35], token3[35], token4[35], token5[35];
    int lineLength = 1025;
    FILE *fid = fopen(filename,"r");
    int N, NZ;
    if(fgets(line, lineLength, fid)==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),EXIT_FAILURE,"error opening matrix file");
    }
    if(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),EXIT_FAILURE,"error reading matrix file header" );
    }
    if (strcmp(token1, "%%MatrixMarket") || strcmp(token2, "matrix") ||
        strcmp(token3, "coordinate") || strcmp(token4, "real") ||
        strcmp(token5, "general"))
    {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),EXIT_FAILURE,"error reading matrix file header");
    }
    // Next, strip off header lines (which start with "%")
    do {
      if(fgets(line, lineLength, fid)==0) {
        if (fid!=0) fclose(fid);
        ML_Exit(Comm.MyPID(),EXIT_FAILURE,"error reading matrix file comments" );
      }
    } while (line[0] == '%');

    // Next get problem dimensions: M, N, NZ
    if(sscanf(line, "%d %d %d", numGlobalRows, &N, &NZ)==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),EXIT_FAILURE,"error reading matrix file dimensions" );
    }
} //ML_Read_Matrix_Dimensions()

//** ************************************************************************ **

void ML_Exit(int mypid, int code, char *fmt,...)
{
  char msg[800];
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg,fmt, ap);
  va_end(ap);
  if (!mypid && msg != 0)
    printf("%s\n",msg);
#ifdef ML_MPI
  MPI_Finalize();
#endif
  exit(code);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) */
