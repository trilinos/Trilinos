//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Version.h"

// prototype
int power_method(Epetra_CrsMatrix& A, double & lambda, int niters, double tolerance,
		 bool verbose);

 
int main(int argc, char *argv[])
{
  int ierr = 0, i;

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

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  bool verbose = (MyPID==0);

  if (verbose)
    cout << Epetra_Version() << endl << endl;

  cout << Comm << endl;

  // Get the number of local equations from the command line
  if (argc!=2)
   {
     if (verbose) 
       cout << "Usage: " << argv[0] << " number_of_equations" << endl;
    exit(1);
   }
  int NumGlobalElements = atoi(argv[1]);

  if (NumGlobalElements < NumProc)
      {
     if (verbose)
       cout << "numGlobalBlocks = " << NumGlobalElements 
	    << " cannot be < number of processors = " << NumProc << endl;
     exit(1);
      }

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.

  Epetra_Map Map(NumGlobalElements, 0, Comm);
  
  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  int * MyGlobalElements = new int[NumMyElements];
    Map.MyGlobalElements(MyGlobalElements);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int * NumNz = new int[NumMyElements];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (i=0; i<NumMyElements; i++)
    if (MyGlobalElements[i]==0 || MyGlobalElements[i] == NumGlobalElements-1)
      NumNz[i] = 2;
    else
      NumNz[i] = 3;

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
  
  for (i=0; i<NumMyElements; i++)
    {
    if (MyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
      {
	Indices[0] = NumGlobalElements-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = MyGlobalElements[i]-1;
	Indices[1] = MyGlobalElements[i]+1;
	NumEntries = 2;
      }
     ierr = A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
     assert(ierr==0);
     // Put in the diagonal entry
     ierr = A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i);
     assert(ierr==0);
    }
   
  // Finish up
  ierr = A.TransformToLocal();
  assert(ierr==0);

  // Create vectors for Power method


  // variable needed for iteration
  double lambda = 0.0;
  int niters = NumGlobalElements*10;
  double tolerance = 1.0e-3;

  // Iterate
  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  Epetra_Time timer(Comm);
  ierr += power_method(A, lambda, niters, tolerance, verbose);
  double elapsed_time = timer.ElapsedTime();
  double total_flops =counter.Flops();
  double MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
    cout << "\n\nTotal MFLOPs for first solve = " << MFLOPs << endl<< endl;

  // Increase diagonal dominance
  if (verbose) 
    cout << "\nIncreasing magnitude of first diagonal term, solving again\n\n"
		    << endl;

  if (A.MyGlobalRow(0)) {
    int numvals = A.NumGlobalEntries(0);
    double * Rowvals = new double [numvals];
    int    * Rowinds = new int    [numvals];
    A.ExtractGlobalRowCopy(0, numvals, numvals, Rowvals, Rowinds); // Get A[0,0]
    for (i=0; i<numvals; i++) if (Rowinds[i] == 0) Rowvals[i] *= 10.0;

    A.ReplaceGlobalValues(0, numvals, Rowvals, Rowinds);
    delete [] Rowvals;
    delete [] Rowinds;
  }
 
  // Iterate (again)
  lambda = 0.0;
  timer.ResetStartTime();
  counter.ResetFlops();
  ierr += power_method(A, lambda, niters, tolerance, verbose);
  elapsed_time = timer.ElapsedTime();
  total_flops = counter.Flops();
  MFLOPs = total_flops/elapsed_time/1000000.0;

  if (verbose) 
    cout << "\n\nTotal MFLOPs for second solve = " << MFLOPs << endl<< endl;


  // Release all objects
  delete [] NumNz;
  delete [] Values;
  delete [] Indices;
  delete [] MyGlobalElements;
			
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}

int power_method(Epetra_CrsMatrix& A, double &lambda, int niters, 
		 double tolerance, bool verbose) {  

  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());
  Epetra_Vector resid(A.RowMap());

  Epetra_Flops * counter = A.GetFlopCounter();
  if (counter!=0) {
    q.SetFlopCounter(A);
    z.SetFlopCounter(A);
    resid.SetFlopCounter(A);
  }

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
      q.Dot(z, &lambda); // Approximate maximum eigenvalue
      if (iter%100==0 || iter+1==niters)
	{
	  resid.Update(1.0, z, -lambda, q, 0.0); // Compute A*q - lambda*q
	  resid.Norm2(&residual);
	  if (verbose) cout << "Iter = " << iter << "  Lambda = " << lambda 
			    << "  Residual of A*q - lambda*q = " 
			    << residual << endl;
	} 
      if (residual < tolerance) {
	ierr = 0;
	break;
      }
    }
  return(ierr);
}
