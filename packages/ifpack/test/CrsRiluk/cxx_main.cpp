//@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "Ifpack_ConfigDefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Object.h"
#include "Ifpack_Version.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


// prototype
int power_method(Epetra_CrsMatrix& A, 
		Epetra_Vector& q,
		Epetra_Vector& z, 
		Epetra_Vector& resid, 
		double * lambda, int niters, double tolerance,
		bool verbose);


int main(int argc, char *argv[])
{
	int ierr = 0, i, j;

#ifdef EPETRA_MPI

	// Initialize MPI

	MPI_Init(&argc,&argv);
	int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if (verbose && MyPID==0)
    cout << Ifpack_Version() << endl << endl;

  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyEquations = 10000;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyEquations++;
  int IndexBase = 0;
  int ElementSize = 7;
  bool DistributedGlobal = (NumGlobalEquations>NumMyEquations);

  // Construct a Map that puts approximately the same Number of equations on each processor

  Epetra_Map Map(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int * MyGlobalElements = new int[Map.NumMyElements()];
  Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor

  int * NumNz = new int[NumMyEquations];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (i=0; i<NumMyEquations; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalEquations-1)
      NumNz[i] = 1;
    else
      NumNz[i] = 2;

  // Create a Epetra_Matrix

  Epetra_CrsMatrix A(Copy, Map, NumNz);
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
  int NumEntries;
  
  for (i=0; i<NumMyEquations; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (MyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
     assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i)>0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // Create vectors for Power method

  Epetra_Vector q(Map);
  Epetra_Vector z(Map);
  Epetra_Vector resid(Map);

  // variable needed for iteration
  double lambda = 0.0;
  int niters = 10000;
  double tolerance = 1.0e-3;

  // Iterate
  Epetra_Time timer(Comm);
  ierr += power_method(A, q, z, resid, &lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();
  double total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for first solve = " << MFLOPs << endl<< endl;

  // Increase diagonal dominance
  if (verbose) cout << "\n\nIncreasing the magnitude of first diagonal term and solving again\n\n"
		    << endl;

  
  if (A.MyGlobalRow(0)) {
    int numvals = A.NumGlobalEntries(0);
    double * Rowvals = new double [numvals];
    int    * Rowinds = new int    [numvals];
    A.ExtractGlobalRowCopy(0, numvals, numvals, Rowvals, Rowinds); // Get A[0,0]

    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;
    
    A.ReplaceGlobalValues(0, numvals, Rowvals, Rowinds);
  }
  // Iterate (again)
  lambda = 0.0;
  timer.ResetStartTime();
  A.ResetFlops(); q.ResetFlops(); z.ResetFlops(); resid.ResetFlops();
  ierr += power_method(A, q, z, resid, &lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = A.Flops() + q.Flops() + z.Flops() + resid.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) cout << "\n\nTotal MFLOPs for second solve = " << MFLOPs << endl<< endl;


  // Release all objects
  delete [] NumNz;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;

			

  if (verbose1) {
    // Test ostream << operator (if verbose1)
    // Construct a Map that puts 2 equations on each PE
    
    int NumMyElements1 = 2;
    int NumMyEquations1 = NumMyElements1;
    int NumGlobalEquations1 = NumMyEquations1*NumProc;

    Epetra_Map Map1(-1, NumMyElements1, 0, Comm);
    
    // Get update list and number of local equations from newly created Map
    int * MyGlobalElements1 = new int[Map1.NumMyElements()];
    Map1.MyGlobalElements(MyGlobalElements1);
    
    // Create an integer vector NumNz that is used to build the Petra Matrix.
    // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation on this processor
    
    int * NumNz1 = new int[NumMyEquations1];
    
    // We are building a tridiagonal matrix where each row has (-1 2 -1)
    // So we need 2 off-diagonal terms (except for the first and last equation)
    
    for (i=0; i<NumMyEquations1; i++)
      if (MyGlobalElements1[i]==0 || MyGlobalElements1[i] == NumGlobalEquations1-1)
	NumNz1[i] = 1;
      else
	NumNz1[i] = 2;
    
    // Create a Epetra_Matrix
    
    Epetra_CrsMatrix A1(Copy, Map1, NumNz1);
    
    // Add  rows one-at-a-time
    // Need some vectors to help
    // Off diagonal Values will always be -1
    
    
    double *Values1 = new double[2];
    Values1[0] = -1.0; Values1[1] = -1.0;
    int *Indices1 = new int[2];
    double two1 = 2.0;
    int NumEntries1;
    
    for (i=0; i<NumMyEquations1; i++)
      {
	if (MyGlobalElements1[i]==0)
	  {
	    Indices1[0] = 1;
	    NumEntries1 = 1;
	  }
	else if (MyGlobalElements1[i] == NumGlobalEquations1-1)
	  {
	    Indices1[0] = NumGlobalEquations1-2;
	    NumEntries1 = 1;
	  }
	else
	  {
	    Indices1[0] = MyGlobalElements1[i]-1;
	    Indices1[1] = MyGlobalElements1[i]+1;
	    NumEntries1 = 2;
	  }
	assert(A1.InsertGlobalValues(MyGlobalElements1[i], NumEntries1, Values1, Indices1)==0);
	assert(A1.InsertGlobalValues(MyGlobalElements1[i], 1, &two1, MyGlobalElements1+i)>0); // Put in the diagonal entry
      }
    
    // Finish up
    assert(A1.TransformToLocal()==0);
    
    if (verbose) cout << "\n\nPrint out tridiagonal matrix, each part on each processor.\n\n" << endl;
    cout << A1 << endl;
    
  // Release all objects
  delete [] NumNz1;
  delete [] Values1;
  delete [] Indices1;
  delete [] MyGlobalElements1;

  }
			
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(Epetra_CrsMatrix& A, 
		 Epetra_Vector& q,
		 Epetra_Vector& z, 
		 Epetra_Vector& resid, 
		 double * lambda, int niters, double tolerance,
		 bool verbose) {  

  // Fill z with random Numbers
  z.Random();

  // variable needed for iteration
  double normz, residual;

  int ierr = 1;

  for (int iter = 0; iter < niters; iter++)
    {
      z.Norm2(&normz); // Compute 2-norm of z
      q.Scale(1.0/normz, z);
      A.Multiply(false, q, z); // Compute z = A*q
      q.Dot(z, lambda); // Approximate maximum eigenvaluE
      if (iter%100==0 || iter+1==niters)
	{
	  resid.Update(1.0, z, -(*lambda), q, 0.0); // Compute A*q - lambda*q
	  resid.Norm2(&residual);
	  if (verbose) cout << "Iter = " << iter << "  Lambda = " << *lambda 
			     << "  Residual of A*q - lambda*q = " << residual << endl;
	} 
      if (residual < tolerance) {
	ierr = 0;
	break;
      }
    }
  return(ierr);
}
