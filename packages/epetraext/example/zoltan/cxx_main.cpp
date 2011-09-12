//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

// EpetraExt::CrsGraph_Zoltan Example routine
#include <Epetra_ConfigDefs.h>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

#include "Trilinos_Util.h"

#include "EpetraExt_Zoltan_CrsGraph.h"
#include "EpetraExt_SymmRCM_CrsGraph.h"
//#include "EpetraExt_ZoltanOrder_CrsGraph.h"
#include "EpetraExt_LPTrans_From_GraphTrans.h"
#include "EpetraExt_Version.h"

#define perror(str) { fprintf(stderr,"%s\n",str);  exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);  exit(ierr); }

int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>2) if (argv[2][0]=='-' && argv[2][1]=='v') verbose = true;
  if (argc<2) perror("error: enter name of data file on cmd line");

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  if (!verbose) Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose) {
    cout << EpetraExt::EpetraExt_Version() << endl << endl;
    cout << Comm << endl << flush;
  }

  Comm.Barrier();
  bool verbose_all = verbose;

  if (verbose) verbose = (MyPID==0);

  //Read in Matrix File and distribute
  int NumGlobalEqs;
  int NumLocalEqs;
  int NumNZs;
  int *NumRowNZs;
  double *Values;
  double *Xprime;
  double *B;
  double *X;
  int *Bindx;
  int NUpdate;
  int *Update;

  Trilinos_Util_read_hb(argv[1], Comm.MyPID(), &NumGlobalEqs, &NumNZs, &Values, &Bindx, &Xprime, &B, &X );
  Trilinos_Util_distrib_msr_matrix( Comm, &NumGlobalEqs, &NumNZs, &NUpdate, &Update, &Values, &Bindx, &Xprime, &B, &X );

  NumLocalEqs = NUpdate;
  NumRowNZs = new int[NumLocalEqs];
  for( int i = 0; i < NumLocalEqs; ++i ) NumRowNZs[i] = Bindx[i+1]-Bindx[i]+1;

  Epetra_Map Map(NumGlobalEqs,NumLocalEqs,Update,0,Comm);

  if( verbose ) cout << "Building Epetra_CrsMatrix" << endl;

  Epetra_CrsMatrix A(Copy, Map, NumRowNZs );

  //Add individual rows
  double *RowVals;
  int *ColInds;
  int NumEntries;
  for( int i = 0; i < NumLocalEqs; ++i )
  {
    RowVals = Values + Bindx[i];
    ColInds = Bindx + Bindx[i];
    NumEntries = Bindx[i+1] - Bindx[i];
    ierr = A.InsertGlobalValues( Update[i], NumEntries, RowVals, ColInds );
    if( ierr ) { printf("Row %d:", Update[i] ); perror1("Error Putting Row: ",ierr); }
    ierr = A.InsertGlobalValues( Update[i], 1, Values+i, Update+i);
    if( ierr ) { perror1("Error Putting Diag: ",ierr); }
  }

  ierr = A.FillComplete();
  if( ierr ) perror1("Error in FillComplete",ierr);

  Epetra_Vector XX(Copy,Map,X);
  Epetra_Vector BB(Copy,Map,B);

  Epetra_LinearProblem Prob(&A,&XX,&BB);

  //Generate Zoltan Load Balanced Version of Linear Problem
  if( verbose ) cout << "Creating Zoltan Partitioning Transform!\n";

  EpetraExt::Zoltan_CrsGraph * ZoltanTrans = new EpetraExt::Zoltan_CrsGraph();
  EpetraExt::LinearProblem_GraphTrans * ZoltanLPTrans =
    new EpetraExt::LinearProblem_GraphTrans( 
         *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(ZoltanTrans)) );

  if( verbose ) cout << "Creating Load Balanced Linear Problem\n";
  Epetra_LinearProblem &BalancedProb = (*ZoltanLPTrans)(Prob);

  EpetraExt::CrsGraph_SymmRCM * RCMTrans = new EpetraExt::CrsGraph_SymmRCM();
  EpetraExt::LinearProblem_GraphTrans * RCMLPTrans =
    new EpetraExt::LinearProblem_GraphTrans( 
         *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(RCMTrans)) );

  if( verbose ) cout << "Creating SymmRCMed Linear Problem\n";
  Epetra_LinearProblem &RCMProb = (*RCMLPTrans)(BalancedProb);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

