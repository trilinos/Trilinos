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

// AMD Test routine
#include <EpetraExt_ConfigDefs.h>
#include "EpetraExt_Version.h"

#ifdef HAVE_MPI
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
#include "EpetraExt_LPTrans_From_GraphTrans.h"
#include "../epetra_test_err.h"

#ifdef HAVE_EXPERIMENTAL
#include "EpetraExt_AMD_CrsGraph.h"
#endif

#include "Teuchos_RCP.hpp"

#define perror(str) { fprintf(stderr,"%s\n",str);  exit(-1); }
#define perror1(str,ierr) { fprintf(stderr,"%s %d\n",str,ierr);  exit(ierr); }

int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef HAVE_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = true;
  bool failure = false;

#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  bool verbose_all = verbose;

  if (verbose) verbose = (MyPID==0);

  if (verbose) cout << EpetraExt::EpetraExt_Version() << endl << endl;

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

  if( verbose ) cout << endl << "Building Epetra_CrsMatrix" << endl;

  Epetra_CrsMatrix A(Copy, Map, NumRowNZs );

  //Add individual rows
  double *RowVals;
  int *ColInds;
  for( int i = 0; i < NumLocalEqs; ++i )
  {
    RowVals = Values + Bindx[i];
    ColInds = Bindx + Bindx[i];
    int NumEntries = Bindx[i+1] - Bindx[i];
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

  // Generate Zoltan Load Balanced Version of Linear Problem
  if( verbose ) cout << "Creating Zoltan Partitioning Transform!\n";

  Teuchos::RCP<EpetraExt::Zoltan_CrsGraph> ZoltanTrans = Teuchos::rcp( new EpetraExt::Zoltan_CrsGraph() );
  if (ZoltanTrans==Teuchos::null) failure = true;

  Teuchos::RCP<EpetraExt::LinearProblem_GraphTrans> ZoltanLPTrans =
    Teuchos::rcp( new EpetraExt::LinearProblem_GraphTrans( 
         *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(&*ZoltanTrans)) ) );

  if( verbose ) cout << "Creating Load Balanced Linear Problem\n";
  Epetra_LinearProblem &BalancedProb = (*ZoltanLPTrans)(Prob);

  // Running this transform fwd() and rvs()
  failure = failure || !ZoltanLPTrans->fwd();
  failure = failure || !ZoltanLPTrans->rvs();

#ifdef HAVE_EXPERIMENTAL

  Teuchos::RCP<EpetraExt::CrsGraph_AMD> AMDTrans = Teuchos::rcp( new EpetraExt::CrsGraph_AMD() );
  if (AMDTrans==Teuchos::null) failure = true;

  Teuchos::RCP<EpetraExt::LinearProblem_GraphTrans> AMDLPTrans =
    Teuchos::rcp( new EpetraExt::LinearProblem_GraphTrans( 
         *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(&*AMDTrans)) ) );

  if( verbose ) cout << "Creating AMD Linear Problem\n";
  Epetra_LinearProblem &AMDProb = (*AMDLPTrans)(BalancedProb);

  // Running this transform fwd() and rvs()
  failure = failure || !AMDLPTrans->fwd();
  failure = failure || !AMDLPTrans->rvs();

#endif
 
  // Clean up
  if (NumRowNZs) delete [] NumRowNZs;
  if (Values) free(Values);
  if (Bindx) free(Bindx);
  if (Xprime) free(Xprime);
  if (B) free(B);
  if (X) free(X);
  if (Update) free(Update);

  Comm.Barrier();

  if ( ierr != 0 || failure ) {
    if (verbose)
      std::cout << "End Result: TEST FAILED" << std::endl;
    return -1;
  }
  //
  // Default return value
  //
  if (verbose)
    std::cout << "End Result: TEST PASSED" << std::endl;
  return 0;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}

