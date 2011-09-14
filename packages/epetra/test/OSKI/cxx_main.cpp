//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


int main(int argc, char **argv){
  return 0;
}

#ifdef THIS_SHOULD_NOT_BE_DEFINED

#include "Epetra_Time.h"
#include "Epetra_OskiMatrix.h"
#include "Epetra_OskiVector.h"
#include "Epetra_OskiUtils.h"
#include "Epetra_OskiPermutation.h"

#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"

using namespace Teuchos;

// Turn on timing
//#define ML_SCALING

// prototypes defined after main
void ML_Exit(int mypid, const char *msg, int code);
void ML_Print_Help();
void ML_Read_Matrix_Dimensions(const char *filename, int *numGlobalRows, Epetra_Comm &Comm);

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  int mypid = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
  int mypid = 0;
#endif

  // Read XML input deck
  ParameterList masterList;
  if (argc > 1) {
    if (strncmp("-h",argv[1],2) == 0) {
      cout << "help" << endl;
      ML_Print_Help();
      ML_Exit(mypid,0,EXIT_SUCCESS);
    }
    else {
      int i=0,j;
      FILE* fid = fopen(argv[1],"r");
      if (fid) {
        i++;
        fclose(fid);
      }
      Comm.SumAll(&i, &j, 1);
      if (j!=Comm.NumProc()) {
        cout << "Could not open input file." << endl;
        ML_Print_Help();
        ML_Exit(mypid,0,EXIT_FAILURE);
      }
      FileInputSource fileSrc(argv[1]);
      XMLObject fileXML = fileSrc.getObject();
      XMLParameterListReader ListReader;
      masterList = ListReader.toParameterList(fileXML);
    }
  } else {
    cout << "No input file specified." << endl;
    ML_Print_Help();
    ML_Exit(mypid,0,EXIT_SUCCESS);
}

  ParameterList *fileList, *AztecOOList;
  try {fileList = &(masterList.sublist("data files",true));}
  catch(...) {ML_Exit(mypid,"Missing \"data files\" sublist.",EXIT_FAILURE);}
  try {AztecOOList = &(masterList.sublist("AztecOO"));}
  catch(...) {ML_Exit(mypid,"Missing \"AztecOO\" sublist.",EXIT_FAILURE);}

#ifdef ML_SCALING
   const int ntimers=4;
   enum {total, probBuild, precBuild, solve};
   ml_DblLoc timeVec[ntimers], maxTime[ntimers], minTime[ntimers];

  for (int i=0; i<ntimers; i++) timeVec[i].rank = Comm.MyPID();
  timeVec[total].value = MPI_Wtime();
#endif
  string matrixfile = fileList->get("matrix input file","A.dat");
  const char *datafile = matrixfile.c_str();
  int numGlobalRows;
  ML_Read_Matrix_Dimensions(datafile, &numGlobalRows, Comm);

#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime();
#endif

  // ===================================================== //
  // READ IN MATRICES FROM FILE                            //
  // ===================================================== //

  if (!mypid) printf("reading %s\n",datafile); fflush(stdout);
  Epetra_CrsMatrix *Amat=NULL;
  //Epetra_Map *RowMap=NULL;
  int errCode=0;
  //if (RowMap) errCode=EpetraExt::MatrixMarketFileToCrsMatrix(datafile, *RowMap, Amat);
  //else        errCode=EpetraExt::MatrixMarketFileToCrsMatrix(datafile, Comm, Amat);
  errCode=EpetraExt::MatrixMarketFileToCrsMatrix(datafile, Comm, Amat);
  if (errCode) ML_Exit(mypid,"error reading matrix", EXIT_FAILURE);
  Amat->OptimizeStorage();
    
  Epetra_Vector LHS(Amat->RowMap()); LHS.Random();
  Epetra_Vector RHS(Amat->RowMap()); RHS.PutScalar(0.0);
  Epetra_LinearProblem Problem(Amat, &LHS, &RHS);

#ifdef ML_SCALING
  timeVec[probBuild].value = MPI_Wtime() - timeVec[probBuild].value;
#endif

  // =========================== build preconditioner ===========================
  
#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime();
#endif

  // no preconditioner right now

#ifdef ML_SCALING
  timeVec[precBuild].value = MPI_Wtime() - timeVec[precBuild].value;
#endif

  // =========================== outer solver =============================

#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime();
#endif
  solver.SetParameters(*AztecOOList);
  int maxits = AztecOOList->get("Aztec iterations",250);
  double tol = AztecOOList->get("Aztec tolerance",1e-10);
#ifdef ML_SCALING
  timeVec[solve].value = MPI_Wtime() - timeVec[solve].value;
#endif

  // compute the real residual
  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  delete Amat;
  //delete RowMap;

#ifdef ML_SCALING
  timeVec[total].value = MPI_Wtime() - timeVec[total].value;

  //avg
  double dupTime[ntimers],avgTime[ntimers];
  for (int i=0; i<ntimers; i++) dupTime[i] = timeVec[i].value;
  MPI_Reduce(dupTime,avgTime,ntimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  for (int i=0; i<ntimers; i++) avgTime[i] = avgTime[i]/Comm.NumProc();
  //min
  MPI_Reduce(timeVec,minTime,ntimers,MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
  //max
  MPI_Reduce(timeVec,maxTime,ntimers,MPI_DOUBLE_INT,MPI_MAXLOC,0,MPI_COMM_WORLD);

  if (Comm.MyPID() == 0) {
    printf("timing :  max (pid)  min (pid)  avg\n");
    printf("Problem build         :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[probBuild].value,maxTime[probBuild].rank,
             minTime[probBuild].value,minTime[probBuild].rank,
             avgTime[probBuild]);
    printf("Preconditioner build  :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[precBuild].value,maxTime[precBuild].rank,
             minTime[precBuild].value,minTime[precBuild].rank,
             avgTime[precBuild]);
    printf("Solve                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[solve].value,maxTime[solve].rank,
             minTime[solve].value,minTime[solve].rank,
             avgTime[solve]);
    printf("Total                 :   %2.3e (%d)  %2.3e (%d)  %2.3e \n",
             maxTime[total].value,maxTime[total].rank,
             minTime[total].value,minTime[total].rank,
             avgTime[total]);
  }
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif


  return(EXIT_SUCCESS);
} //main

void ML_Exit(int mypid, const char *msg, int code)
{
  if (!mypid && msg != 0)
    printf("%s\n",msg);
#ifdef ML_MPI
  MPI_Finalize();
#endif
  exit(code);
}

void ML_Print_Help()
{
  printf("Usage: ml_scaling.exe [XML input file]\n");
  printf("  The XML input file must have three sublists:\n");
  printf("    1) \"data files\"               -- input file names\n");
  printf("    2) \"AztecOO\"                  -- outer Krylov options\n");
}

void ML_Read_Matrix_Dimensions(const char *filename, int *numGlobalRows, Epetra_Comm &Comm)
{
    char line[35], token1[35], token2[35], token3[35], token4[35], token5[35];
    int lineLength = 1025;
    FILE *fid = fopen(filename,"r");
    int N, NZ;
    if(fgets(line, lineLength, fid)==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),"error opening matrix file", EXIT_FAILURE);
    }
    if(sscanf(line, "%s %s %s %s %s", token1, token2, token3, token4, token5 )==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),"error reading matrix file header", EXIT_FAILURE);
    }
    if (strcmp(token1, "%%MatrixMarket") || strcmp(token2, "matrix") ||
        strcmp(token3, "coordinate") || strcmp(token4, "real") ||
        strcmp(token5, "general"))
    {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),"error reading matrix file header", EXIT_FAILURE);
    }
    // Next, strip off header lines (which start with "%")
    do {
      if(fgets(line, lineLength, fid)==0) {
        if (fid!=0) fclose(fid);
        ML_Exit(Comm.MyPID(),"error reading matrix file comments", EXIT_FAILURE);
      }
    } while (line[0] == '%');

    // Next get problem dimensions: M, N, NZ
    if(sscanf(line, "%d %d %d", numGlobalRows, &N, &NZ)==0) {
      if (fid!=0) fclose(fid);
      ML_Exit(Comm.MyPID(),"error reading matrix file dimensions", EXIT_FAILURE);
    }
} //ML_Read_Matrix_Dimensions()

#endif
