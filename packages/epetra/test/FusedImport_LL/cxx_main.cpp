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


#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_Flops.h"

#ifdef EPETRA_MPI

#include "Epetra_MpiComm.h"
#include "mpi.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

// prototypes

int check(Epetra_CrsMatrix& A, int NumMyRows1, int NumGlobalRows1, int NumMyNonzeros1,
	  int NumGlobalNonzeros1, int * MyGlobalElements, bool verbose);

int power_method(bool TransA, Epetra_CrsMatrix& A, 
		 Epetra_Vector& q,
		 Epetra_Vector& z, 
		 Epetra_Vector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose);

int check_graph_sharing(Epetra_Comm& Comm);


double test_with_matvec(const Epetra_CrsMatrix &A, const Epetra_CrsMatrix &B){
  const Epetra_Map & Xmap = A.DomainMap();
  const Epetra_Map & Ymap = A.RangeMap();

  Epetra_Vector X(Xmap), Y1(Ymap), Y2(Ymap);

  X.SetSeed(24601);
  X.Random();

  A.Apply(X,Y1);

  B.Apply(X,Y2);

  Y1.Update(-1.0,Y2,1.0);
  double norm;
  Y1.Norm2(&norm);

  return norm;
}

int build_matrix_unfused(const Epetra_CrsMatrix & SourceMatrix, Epetra_Import & RowImporter, Epetra_CrsMatrix *&A){
  int rv=0;
  rv=A->Import(SourceMatrix, RowImporter, Insert);
  if(rv) {cerr<<"build_matrix_unfused: Import failed"<<endl; return rv;}

  rv=A->FillComplete(SourceMatrix.DomainMap(), SourceMatrix.RangeMap()); 
  return rv;
}



int build_test_matrix(Epetra_MpiComm & Comm, int test_number, Epetra_CrsMatrix *&A){
  int NumProc = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if(test_number==1){
    // Case 1: Tridiagonal
    //    int NumMyEquations = 10000;
    int NumMyEquations = 100;

    long long NumGlobalEquations = (NumMyEquations * NumProc) + EPETRA_MIN(NumProc,3);
    if(MyPID < 3)  NumMyEquations++;
    
    // Construct a Map that puts approximately the same Number of equations on each processor
    Epetra_Map Map(NumGlobalEquations, NumMyEquations, (long long)0, Comm);
    
    // Get update list and number of local equations from newly created Map
    long long* MyGlobalElements = new long long[Map.NumMyElements()];
    Map.MyGlobalElements(MyGlobalElements);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int* NumNz = new int[NumMyEquations];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (int i = 0; i < NumMyEquations; i++)
      if((MyGlobalElements[i] == 0) || (MyGlobalElements[i] == NumGlobalEquations - 1))
	NumNz[i] = 1;
      else
	NumNz[i] = 2;
    
    // Create a Epetra_Matrix     
    A=new Epetra_CrsMatrix(Copy, Map, NumNz);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    double* Values = new double[2];
    Values[0] = -1.0; 
    Values[1] = -1.0;
    long long* Indices = new long long[2];
    double two = 2.0;
    int NumEntries;
    
    for (int i = 0; i < NumMyEquations; i++) {
      if(MyGlobalElements[i] == 0) {
	Indices[0] = 1;
	NumEntries = 1;
      }
      else if (MyGlobalElements[i] == NumGlobalEquations-1) {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
      else {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
      A->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
      A->InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
    }
    
    A->FillComplete();
    
    // Cleanup
    delete [] MyGlobalElements;
    delete [] NumNz;
    delete [] Values;
    delete [] Indices;
    
  }
  return 0; 
}




int main(int argc, char *argv[])
{
  int ierr = 0, total_err=0;

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose && MyPID==0)
    cout << Epetra_Version() << std::endl << std::endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
		    << " is alive."<<endl;

  // Redefine verbose to only print on PE 0
  if(verbose && rank!=0) verbose = false;

  // Matrix & Map pointers
  Epetra_CrsMatrix *A, *B, *C;
  Epetra_Map* Map1;
  Epetra_Import* Import1;
  double diff_tol=1e-12;

#define ENABLE_TEST_1
#define ENABLE_TEST_2
#define ENABLE_TEST_3
#define ENABLE_TEST_4

  /////////////////////////////////////////////////////////
  // Test #1: Tridiagonal Matrix; Migrate to Proc 0
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_1
  {
    ierr=build_test_matrix(Comm,1,A); 
    long long num_global = A->RowMap().NumGlobalElements64();
    
    // New map with all on Proc1 + importer
    if(MyPID==0) Map1=new Epetra_Map(num_global,num_global,(long long) 0,Comm);
    else         Map1=new Epetra_Map(num_global,0,(long long)0,Comm);
    Import1 = new Epetra_Import(*Map1,A->RowMap());

    // Execute fused constructor
    B=new Epetra_CrsMatrix(*A,*Import1,&A->RangeMap());
    
    double diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #1 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }
    
    delete A; delete B; delete Map1; delete Import1;
  }
#endif


  /////////////////////////////////////////////////////////
  // Test #2: Tridiagonal Matrix; Locally Reversed Map
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_2
  {
    ierr=build_test_matrix(Comm,1,A); 
    int num_local = A->RowMap().NumMyElements();
    
    std::vector<long long> MyGIDS(num_local);
    for(int i=0; i<num_local; i++)
      MyGIDS[i] = A->RowMap().GID64(num_local-i-1);

    // New map with all on Proc1 + importer
    Map1=new Epetra_Map((long long)-1,num_local,&MyGIDS[0],0,Comm);
    Import1 = new Epetra_Import(*Map1,A->RowMap());

    // Execute fused constructor
    B=new Epetra_CrsMatrix(*A,*Import1,&A->RangeMap());

    double diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #2 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1;
  }
#endif

  /////////////////////////////////////////////////////////
  // Test #3: Tridiagonal Matrix; Globally Reversed Map
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_3
  {
    ierr=build_test_matrix(Comm,1,A); 
    int num_local  = A->RowMap().NumMyElements();
    long long num_global = A->RowMap().NumGlobalElements64();
    int num_scansum = 0;

    Comm.ScanSum(&num_local,&num_scansum,1);

    std::vector<long long> MyGIDS(num_local);
    for(int i=0; i<num_local; i++)
      MyGIDS[i] = num_global - num_scansum + num_local - i - 1;

    Map1=new Epetra_Map((long long)-1,num_local,&MyGIDS[0],(long long)0,Comm);
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    
    // Execute fused constructor
    B=new Epetra_CrsMatrix(*A,*Import1,&A->RangeMap());
    
    double diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #3 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1;

  }
#endif


  /////////////////////////////////////////////////////////
  // Test #4: Tridiagonal Matrix; MMM style halo import
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_4
  {
    ierr=build_test_matrix(Comm,1,A); 

    // Assume we always own the diagonal
    int num_local = A->NumMyCols()-A->NumMyRows();
    std::vector<long long> MyGIDS(num_local);

    for(int i=0, idx=0; i<A->NumMyCols(); i++)
      if(A->LRID(A->GCID64(i)) == -1){
	MyGIDS[idx] = A->GCID(i);
	idx++;
      }

    Map1=new Epetra_Map((long long)-1,num_local,&MyGIDS[0],(long long)0,Comm);
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    
    // Execute fused constructor
    B=new Epetra_CrsMatrix(*A,*Import1,&A->RangeMap());

    // Build unfused matrix to compare
    C=new Epetra_CrsMatrix(Copy,*Map1,0);
    build_matrix_unfused(*A,*Import1,C);

    double diff=test_with_matvec(*B,*C);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #4 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete C; delete Map1; delete Import1;
  }
#endif


  // Final output for OK
  if(MyPID==0 && total_err==0)
    cout<<"FusedImport: All tests PASSED."<<endl;

    // Cleanup  
  MPI_Finalize();

  return total_err ;
}



#else
int main(){
  
  return 0;
}
#endif
