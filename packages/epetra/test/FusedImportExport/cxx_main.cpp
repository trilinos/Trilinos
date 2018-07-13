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
#include "Epetra_Export.h"
#include "Epetra_Util.h"
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
  const Epetra_Map & Xamap  = A.DomainMap();
  const Epetra_Map & Yamap  = A.RangeMap();
  const Epetra_Map & Xbmap  = B.DomainMap();
  const Epetra_Map & Ybmap  = B.RangeMap();

  Epetra_Vector Xa(Xamap), Xb(Xbmap), Ya(Yamap), Yb(Ybmap), Diff(Yamap);

  Xa.SetSeed(24601);
  Xa.Random();

  // Handle domain map change
  if(!Xamap.SameAs(Xbmap)) {
    Epetra_Import Ximport(Xbmap,Xamap);
    Xb.Import(Xa,Ximport,Insert);
  }
  else {
    Xb=Xa;
  }

  // Do the multiplies
  A.Apply(Xa,Ya);
  B.Apply(Xb,Yb);

  // Handle Rangemap change
  if(!Yamap.SameAs(Ybmap)) {
    Epetra_Import Yimport(Yamap,Ybmap);
    Diff.Import(Yb,Yimport,Insert);
  }
  else {
    Diff=Yb;
  }

  // Check solution
  Diff.Update(-1.0,Ya,1.0);
  double norm;
  Diff.Norm2(&norm);

  return norm;
}


// B here is the "reduced" matrix.  Square matrices w/ Row=Domain=Range only.
double test_with_matvec_reduced_maps(const Epetra_CrsMatrix &A, const Epetra_CrsMatrix &B, const Epetra_Map & Bfullmap){
  const Epetra_Map & Amap  = A.DomainMap();
  Epetra_Vector Xa(Amap), Ya(Amap), Diff(Amap);
  const Epetra_Map *Bmap  = Bfullmap.NumMyElements() > 0 ? &B.DomainMap() : 0;
  Epetra_Vector *Xb = Bmap ? new Epetra_Vector(*Bmap) : 0;
  Epetra_Vector *Yb = Bmap ? new Epetra_Vector(*Bmap) : 0;

  Epetra_Vector Xb_alias(View,Bfullmap, Bmap ? Xb->Values(): 0);
  Epetra_Vector Yb_alias(View,Bfullmap, Bmap ? Yb->Values(): 0);

  Epetra_Import Ximport(Bfullmap,Amap);

  // Set the input vector
  Xa.SetSeed(24601);
  Xa.Random();
  Xb_alias.Import(Xa,Ximport,Insert);

  // Do the multiplies
  A.Apply(Xa,Ya);
  if(Bmap) B.Apply(*Xb,*Yb);

  // Check solution
  Epetra_Import Yimport(Amap,Bfullmap);
  Diff.Import(Yb_alias,Yimport,Insert);


  Diff.Update(-1.0,Ya,1.0);
  double norm;
  Diff.Norm2(&norm);

  delete Xb; delete Yb;
  return norm;
}



int build_matrix_unfused(const Epetra_CrsMatrix & SourceMatrix, Epetra_Import & RowImporter, Epetra_CrsMatrix *&A){
  int rv=0;
  rv=A->Import(SourceMatrix, RowImporter, Insert);
  if(rv) {cerr<<"build_matrix_unfused: Import failed"<<endl; return rv;}

  rv=A->FillComplete(SourceMatrix.DomainMap(), SourceMatrix.RangeMap());
  return rv;
}

int build_matrix_unfused(const Epetra_CrsMatrix & SourceMatrix, Epetra_Export & RowExporter, Epetra_CrsMatrix *&A){
  int rv=0;
  rv=A->Export(SourceMatrix, RowExporter, Insert);
  if(rv) {cerr<<"build_matrix_unfused: Export failed"<<endl; return rv;}

  rv=A->FillComplete(SourceMatrix.DomainMap(), SourceMatrix.RangeMap());
  return rv;
}



void build_test_matrix(Epetra_MpiComm & Comm, int test_number, Epetra_CrsMatrix *&A){
  int NumProc = Comm.NumProc();
  int MyPID   = Comm.MyPID();

  if(test_number==1){
    // Case 1: Tridiagonal
    int NumMyEquations = 100;

    int NumGlobalEquations = (NumMyEquations * NumProc) + EPETRA_MIN(NumProc,3);
    if(MyPID < 3)  NumMyEquations++;

    // Construct a Map that puts approximately the same Number of equations on each processor
    Epetra_Map Map(NumGlobalEquations, NumMyEquations, 0, Comm);

    // Get update list and number of local equations from newly created Map
    int* MyGlobalElements = new int[Map.NumMyElements()];
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
    int* Indices = new int[2];
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
}


void build_test_map(const Epetra_Map & oldMap, Epetra_Map *& newMap) {
  int NumProc = oldMap.Comm().NumProc();
  int MyPID   = oldMap.Comm().MyPID();

  int num_global = oldMap.NumGlobalElements();
  if(NumProc<3) {
    // Dump everything onto -proc 0
    int num_local = MyPID==0 ? num_global : 0;
    newMap = new Epetra_Map(num_global,num_local,0,oldMap.Comm());
  }
  else {
    // Split everything between procs 0 and 2 (leave proc 1 empty)
    int num_local=0;
    if(MyPID==0) num_local = num_global/2;
    else if(MyPID==2) num_local =  num_global - ((int)num_global/2);
    newMap = new Epetra_Map(num_global,num_local,0,oldMap.Comm());
  }
}


int main(int argc, char *argv[])
{
  int total_err=0;

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
  Epetra_Export* Export1;
  double diff_tol=1e-12;

#define ENABLE_TEST_1
#define ENABLE_TEST_2
#define ENABLE_TEST_3
#define ENABLE_TEST_4
#define ENABLE_TEST_5
#define ENABLE_TEST_6

  /////////////////////////////////////////////////////////
  // Test #1: Tridiagonal Matrix; Migrate to Proc 0
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_1
  {
    double diff;
    build_test_matrix(Comm,1,A);
    int num_global = A->RowMap().NumGlobalElements();

    // New map with all on Proc1
    if(MyPID==0) Map1=new Epetra_Map(num_global,num_global,0,Comm);
    else         Map1=new Epetra_Map(num_global,0,0,Comm);

    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #1 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #1 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1; delete Export1;
  }
#endif


  /////////////////////////////////////////////////////////
  // Test #2: Tridiagonal Matrix; Locally Reversed Map
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_2
  {
    double diff;
    build_test_matrix(Comm,1,A);
    int num_local = A->RowMap().NumMyElements();

    std::vector<int> MyGIDS(num_local);
    for(int i=0; i<num_local; i++)
      MyGIDS[i] = A->RowMap().GID(num_local-i-1);

    // New map with all on Proc1
    Map1=new Epetra_Map(-1,num_local,&MyGIDS[0],0,Comm);

    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #2 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #2 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1; delete Export1;
  }
#endif

  /////////////////////////////////////////////////////////
  // Test #3: Tridiagonal Matrix; Globally Reversed Map
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_3
  {
    double diff;
    build_test_matrix(Comm,1,A);
    int num_local  = A->RowMap().NumMyElements();
    int num_global = A->RowMap().NumGlobalElements();
    int num_scansum = 0;

    Comm.ScanSum(&num_local,&num_scansum,1);

    // New Map
    std::vector<int> MyGIDS(num_local);
    for(int i=0; i<num_local; i++)
      MyGIDS[i] = num_global - num_scansum + num_local - i - 1;
    Map1=new Epetra_Map(-1,num_local,&MyGIDS[0],0,Comm);


    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #3 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,0,&A->RangeMap());

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #3 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1; delete Export1;
  }
#endif


  /////////////////////////////////////////////////////////
  // Test #4: Tridiagonal Matrix; MMM style halo import
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_4
  {
    double diff;
    build_test_matrix(Comm,1,A);

    // Assume we always own the diagonal
    int num_local = A->NumMyCols()-A->NumMyRows();
    std::vector<int> MyGIDS(num_local);

    for(int i=0, idx=0; i<A->NumMyCols(); i++)
      if(A->LRID(A->GCID(i)) == -1){
	MyGIDS[idx] = A->GCID(i);
	idx++;
      }

    // New map
    const int * MyGIDS_ptr = Epetra_Util_data_ptr(MyGIDS);
    Map1=new Epetra_Map(-1,num_local,MyGIDS_ptr,0,Comm);


    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,0,&A->RangeMap());

    // Build unfused matrix to compare
    C=new Epetra_CrsMatrix(Copy,*Map1,0);
    build_matrix_unfused(*A,*Import1,C);

    diff=test_with_matvec(*B,*C);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #4 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,0,&A->RangeMap());

    diff=test_with_matvec(*B,*C);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #4 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete C; delete Map1; delete Import1; delete Export1;
  }
#endif


  /////////////////////////////////////////////////////////
  // Test 5: Tridiagonal Matrix; Migrate to Proc 0, Replace Maps
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_5
  {
    double diff;
    build_test_matrix(Comm,1,A);

    // New map with all on Procs 0 and 2
    build_test_map(A->RowMap(),Map1);

    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,Map1,Map1);

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #5 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,Map1,Map1);

    diff=test_with_matvec(*A,*B);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #5 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1; delete Export1;
  }
#endif


  /////////////////////////////////////////////////////////
  // Test 6: Tridiagonal Matrix; Migrate to Proc 0, Replace Comm
  /////////////////////////////////////////////////////////
#ifdef ENABLE_TEST_6
  {
    double diff;
    build_test_matrix(Comm,1,A);

    // New map with all on Procs 0 and 2
    build_test_map(A->RowMap(),Map1);

    // Execute fused import constructor
    Import1 = new Epetra_Import(*Map1,A->RowMap());
    B=new Epetra_CrsMatrix(*A,*Import1,Map1,Map1,true);

    diff=test_with_matvec_reduced_maps(*A,*B,*Map1);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedImport: Test #6 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    // Execute fused export constructor
    delete B;
    Export1 = new Epetra_Export(A->RowMap(),*Map1);
    B=new Epetra_CrsMatrix(*A,*Export1,Map1,Map1,true);

    diff=test_with_matvec_reduced_maps(*A,*B,*Map1);
    if(diff > diff_tol){
      if(MyPID==0) cout<<"FusedExport: Test #6 FAILED with norm diff = "<<diff<<"."<<endl;
      total_err--;
    }

    delete A; delete B; delete Map1; delete Import1; delete Export1;
  }
#endif


  // Final output for OK
  if(MyPID==0 && total_err==0)
    cout<<"FusedImportExport: All tests PASSED."<<endl;

  // Cleanup
  MPI_Finalize();

  return total_err ;
}



#else
int main(){

  return 0;
}
#endif
