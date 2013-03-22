// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

// usage: ./driver.exe -f A.dat -rbm rbm.dat -i options.xml

#include <MueLu_ConfigDefs.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_MultiVectorIn.h>

#include "MatrixVectorChecker.hpp"

using namespace Teuchos;

void ML_Read_Matrix_Dimensions(char *filename, int *numGlobalRows,Epetra_Comm &Comm);
// Small function to handle exiting gracefully.  All pids must call.
void ML_Exit(int mypid, int code, const char *fmt, ...);

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
  char matrixFile[80] = "mat_example.mm\0";
  char rbmFile[80] = "\0";
  char coordFile[80] = "\0";
  // char xmlFile[80] = "\0";
  char rhsFile[80] = "\0";
  int numPDEs = 1;

  while (i<argc) {
    if      (strncmp(argv[i],"-f",2) == 0) {strncpy(matrixFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-rbm",4) == 0) {strncpy(rbmFile,argv[i+1],80); i+=2;}
    else if (strncmp(argv[i],"-c",2) == 0) {strncpy(coordFile,argv[i+1],80); i+=2;}
    // else if (strncmp(argv[i],"-i",2) == 0) {strncpy(xmlFile,argv[i+1],80); i+=2;}
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
    for (int ii=0; ii<coordVector->NumVectors(); ii++) xyz[ii] = mv+ii*stride;
  }


  //TODO  MatrixVectorChecker(*A);

  delete [] rbmPointer;
  delete rbmVector;
  delete pRHS;
  if (coordVector) delete coordVector;
  if (xyz) delete [] xyz;

  delete A;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  std::cout << "TEST PASSED" << std::endl;
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

void ML_Exit(int mypid, int code, const char *fmt,...)
{
  char msg[800];
  va_list ap;
  va_start(ap, fmt);
  vsprintf(msg,fmt, ap);
  va_end(ap);
  if (!mypid && msg != 0)
    printf("%s\n",msg);
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  std::cout << "TEST FAILED" << std::endl;
  exit(code);
}
