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


//  Epetra_LinearProblemRedistor Test routine

//
//  This file demonstrates a problem in some destructor somewhere - or 
//  perhaps the proble could be in this test code.  
//
//  In any case, the symptom is a Segmentation fualt 
// 


#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"

#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_LinearProblemRedistor.h"
#include "Trilinos_Util.h"
#ifdef HAVE_COMM_ASSERT_EQUAL
#include "Comm_assert_equal.h"
#endif

int checkResults( bool trans, 
		  Epetra_LinearProblemRedistor * redistor, 
		  Epetra_LinearProblem * OrigProb,
		  Epetra_LinearProblem * RedistProb,
		  bool verbose) ;
  
int main(int argc, char *argv[]) {

{
  int i;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // Uncomment to debug in parallel int tmp; if (comm.MyPID()==0) cin >> tmp; comm.Barrier();

  bool verbose = false;
  bool veryVerbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  if (!verbose) comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  if (verbose) cout << comm << endl << flush;

  bool verbose1 = verbose;
  if (verbose) verbose = (comm.MyPID()==0);

  int nx = 4;
  int ny = comm.NumProc()*nx; // Scale y grid with number of processors
  
  // Create funky stencil to make sure the matrix is non-symmetric (transpose non-trivial):

  // (i-1,j-1) (i-1,j  )
  // (i  ,j-1) (i  ,j  ) (i  ,j+1)
  // (i+1,j-1) (i+1,j  )
  
  int npoints = 2;

  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

  Epetra_Map * map;
  Epetra_CrsMatrix * A;
  Epetra_Vector * x, * b, * xexact;
	
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  if (verbose) 
    cout << "npoints = " << npoints << " nx = " << nx << " ny = " << ny  << endl ; 

  if (verbose && nx<6 ) {
    cout << *A << endl;
    cout << "B       = " << endl << *b << endl;
  }
  // Construct linear problem object
  
  Epetra_LinearProblem origProblem(A, x, b);
  Epetra_LinearProblem *redistProblem;
  
  Epetra_Time timer(comm);
  
  // Construct redistor object, use all processors and replicate full problem on each

  double start = timer.ElapsedTime();
  Epetra_LinearProblemRedistor redistor(&origProblem, comm.NumProc(), true);
  if (verbose) cout << "\nTime to construct redistor  = " 
		    << timer.ElapsedTime() - start << endl;

  bool ConstructTranspose = true; 
  bool MakeDataContiguous = true;
  
  //
  //  Commenting out the following define avoids the segmentation fault
  //
#define CREATEREDISTPROBLEM
#ifdef CREATEREDISTPROBLEM
  start = timer.ElapsedTime();
  redistor.CreateRedistProblem(ConstructTranspose, MakeDataContiguous, redistProblem);
  if (verbose) cout << "\nTime to create redistributed problem = " 
		    << timer.ElapsedTime() - start << endl;
  
#endif
#if 0
  // Now test output of redistor by performing matvecs

  int ierr = 0;
  ierr += checkResults( ConstructTranspose, &redistor, &origProblem, 
			redistProblem, verbose);
#endif
  cout << " Here we are just before the return()  "  << endl ; 
}
 return(0) ;
  
}

//
//  checkResults computes Ax = b1 and either Rx=b2 or R^T x=b2 (depending on trans) where
//  A = origProblem and R = redistProblem.  
//  checkResults returns 0 (OK) if norm(b1-b2)/norm(b1) < size(b1) * 1e-15; 
//
//  I guess that we need the redistor so that we can move x and b2 around.  
//
const double error_tolerance = 1e-15 ; 

int checkResults( bool trans, 
		  Epetra_LinearProblemRedistor * redistor, 
		  Epetra_LinearProblem * A,
		  Epetra_LinearProblem * R,
		  bool verbose) {
  
  int m = A->GetRHS()->MyLength();
  int n = A->GetLHS()->MyLength();
  assert( m == n ) ; 

  Epetra_MultiVector *x = A->GetLHS() ; 
  Epetra_MultiVector x1( *x ) ; 
  //  Epetra_MultiVector Difference( x1 ) ; 
  Epetra_MultiVector *b = A->GetRHS();
  Epetra_RowMatrix *matrixA = A->GetMatrix();
  assert( matrixA != 0 ) ; 
  int iam = matrixA->Comm().MyPID();
  
  //  Epetra_Time timer(A->Comm());
  //  double start = timer.ElapsedTime();

  matrixA->Multiply(trans, *b, x1) ;   // x = Ab 

  int M,N,nz;
  int *ptr, *ind;
  double *val, *rhs, *lhs;
  int Nrhs, ldrhs, ldlhs;

  redistor->ExtractHbData( M, N, nz, ptr, ind, 
		    val, Nrhs, rhs, ldrhs, 
		    lhs, ldlhs);

  assert( M == N ) ; 
  if ( verbose ) {
    cout << " iam = " << iam 
	 << " m = " << m  << " n = " << n  << " M = " << M << endl ;  
    
    cout << " iam = " << iam << " ptr = " << ptr[0] << " "   << ptr[1] << " "  << ptr[2] << " "  << ptr[3] << " "  << ptr[4] << " " << ptr[5] << endl ;
    
    cout << " iam = " << iam << " ind = " << ind[0] << " "   << ind[1] << " "  << ind[2] << " "  << ind[3] << " "  << ind[4] << " " << ind[5] << endl ;
    
    cout << " iam = " << iam << " val = " << val[0] << " "   << val[1] << " "  << val[2] << " "  << val[3] << " "  << val[4] << " " << val[5] << endl ;
  }
  //  Create a serial map in case we end up needing it 
  //  If it is created inside the else block below it would have to
  //  be with a call to new().
  int NumMyElements_ = 0 ;
  if (matrixA->Comm().MyPID()==0) NumMyElements_ = n;
  Epetra_Map SerialMap( n, NumMyElements_, 0, matrixA->Comm() );

  //  These are unnecessary and useless
  //  Epetra_Vector serial_A_rhs( SerialMap ) ; 
  //  Epetra_Vector serial_A_lhs( SerialMap ) ;

  //  Epetra_Export exporter( matrixA->BlockRowMap(), SerialMap) ;

  //
  //  In each process, we will compute Rb putting the answer into LHS
  //


  for ( int k = 0 ; k < Nrhs; k ++ ) { 
    for ( int i = 0 ; i < M ; i ++ ) { 
      lhs[ i + k * ldlhs ] = 0.0; 
    }
    for ( int i = 0 ; i < M ; i++ ) { 
      for ( int l = ptr[i]; l < ptr[i+1]; l++ ) {
	int j = ind[l] ; 
	if ( verbose && N < 40 ) {
	  cout << " i = " << i << " j = " << j ;
	  cout << " l = " << l << " val[l] = " << val[l] ;
	  cout << " rhs = " << rhs[ j + k * ldrhs ] << endl ;
	}
	lhs[ i + k * ldrhs ] += val[l] * rhs[ j + k * ldrhs ] ;
      }
    }

    if ( verbose && N < 40 ) { 
      cout << " lhs = " ; 
      for ( int j = 0 ; j < N ; j++ ) cout << " " << lhs[j] ; 
      cout << endl ; 
      cout << " rhs = " ; 
      for ( int j = 0 ; j < N ; j++ ) cout << " " << rhs[j] ; 
      cout << endl ; 
    }

    const Epetra_Comm &comm = matrixA->Comm() ; 
#ifdef HAVE_COMM_ASSERT_EQUAL
    //
    //  Here we double check to make sure that lhs and rhs are 
    //  replicated.  
    //
    for ( int j = 0 ; j < N ; j++ ) { 
      assert( Comm_assert_equal( &comm, lhs[ j + k * ldrhs ] ) ) ; 
      assert( Comm_assert_equal( &comm, rhs[ j + k * ldrhs ] ) ) ; 
    } 
#endif
  }

  //
  //  Now we have to redistribue them back 
  //
  redistor->UpdateOriginalLHS( A->GetLHS() ) ; 
  //
  //  Now we want to compare x and x1 which have been computed as follows:
  //  x = Rb  
  //  x1 = Ab 
  //

  double Norm_x1, Norm_diff ;
  EPETRA_CHK_ERR( x1.Norm2( &Norm_x1 ) ) ; 

  //  cout << " x1 = " << x1 << endl ; 
  //  cout << " *x = " << *x << endl ; 

  x1.Update( -1.0, *x, 1.0 ) ; 
  EPETRA_CHK_ERR( x1.Norm2( &Norm_diff ) ) ; 

  //  cout << " diff, i.e. updated x1 = " << x1 << endl ; 

  int ierr = 0;

  if ( verbose ) {
    cout << " Norm_diff = "  << Norm_diff << endl ; 
    cout << " Norm_x1 = "  << Norm_x1 << endl ; 
  }

  if ( Norm_diff / Norm_x1 > n * error_tolerance ) ierr++ ; 

  if (ierr!=0 && verbose) cerr << "Status: Test failed" << endl;
  else if (verbose) cerr << "Status: Test passed" << endl;

  return(ierr);
}
