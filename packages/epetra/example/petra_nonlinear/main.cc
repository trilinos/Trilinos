//@HEADER
/*
************************************************************************

              Epetra: Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_DCRS_Matrix.h"
#ifdef PETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif

int main(int argc, char *argv[])
{
  int i;

#ifdef PETRA_MPI
  MPI_Init(&argc,&argv);
#endif

  // get number of processors and the name of this processor
 
#ifdef PETRA_MPI
  Petra_Comm& comm = *new Petra_Comm(MPI_COMM_WORLD);
#else
  Petra_Comm& comm = *new Petra_Comm();
#endif

  int NumProc = comm.getNumProc();
  int MyPID   = comm.getMyPID();
  cout << "Processor " << MyPID << " of " <<  NumProc << " is alive." << endl;

  // Get the number of local equations from the command line
  if (argc!=2)
   {
     if (MyPID==0) cout << "Usage: " << argv[0] << " number_of_equations" << endl;
    exit(1);
   }
  int numGlobalEquations = atoi(argv[1]);

  if (numGlobalEquations < NumProc)
      {
     if (MyPID==0)
       cout << "numGlobalBlocks = " << numGlobalEquations 
	    << " cannot be < number of processors = " << NumProc << endl;
     exit(1);
      }

  // Construct a map that puts approximately the same number of equations on each processor

  Petra_Map& map = *new Petra_Map(numGlobalEquations, comm);
  
  // Get update list and number of local equations from newly created map
  int * UpdateList = map.getUpdateList();
  int numLocalEquations = map.numLocalEquations();

  // Create an integer vector numNz that is used to build the Petra Matrix.
  // numNz[i] is the number of OFF-DIAGONAL term for the ith global equation on this processor

  int * numNz = new int[numLocalEquations];

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (i=0; i<numLocalEquations; i++)
    if (UpdateList[i]==0 || UpdateList[i] == numGlobalEquations-1)
      numNz[i] = 1;
    else
      numNz[i] = 2;

  // Create a Petra_Matrix

  Petra_RDP_DCRS_Matrix& A = *new Petra_RDP_DCRS_Matrix(map);
  
  // Allocate space using numNz
  
  assert(A.allocate(numNz)==0);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal values will always be -1


  double *values = new double[2];
  values[0] = -1.0; values[1] = -1.0;
  int *indices = new int[2];
  double two = 2.0;
  int numEntries;
  
  for (i=0; i<numLocalEquations; i++)
    {
    if (UpdateList[i]==0)
      {
	indices[0] = 1;
	numEntries = 1;
      }
    else if (UpdateList[i] == numGlobalEquations-1)
      {
	indices[0] = numGlobalEquations-2;
	numEntries = 1;
      }
    else
      {
	indices[0] = UpdateList[i]-1;
	indices[1] = UpdateList[i]+1;
	numEntries = 2;
      }
     assert(A.putRow(UpdateList[i], numEntries, values, indices)==0);
     assert(A.putRow(UpdateList[i], 1, &two, UpdateList+i)==0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(A.fillComplete()==0);

  // Create vectors for Power method

  Petra_RDP_Vector& q = *new Petra_RDP_Vector(map);
  Petra_RDP_Vector& z = *new Petra_RDP_Vector(map);
  Petra_RDP_Vector& resid = *new Petra_RDP_Vector(map);

  // Fill z with random numbers
  z.random();

  // variable needed for iteration
  double normz, lambda, residual;

  // Iterate
  int niters = 500*numGlobalEquations;
  double tolerance = 1.0e-10;
  for (int iter = 0; iter < niters; iter++)
    {
      z.norm2(&normz); // Compute 2-norm of z
      q.scaleCopy(z, 1.0/normz);
      A.matvec(q, z); // Compute z = A*q
      q.dotProd(z, &lambda); // Approximate maximum eigenvaluE
      if (iter%100==0 || iter+1==niters)
	{
	  resid.linComb(z, -lambda, q); // Compute A*q - lambda*q
	  resid.norm2(&residual);
	  if (MyPID==0) cout << "Iter = " << iter << "  Lambda = " << lambda 
			     << "  Residual of A*q - lambda*q = " << residual << endl;
	} 
      if (residual < tolerance) break;
    }
  
  // Release all objects

  delete [] numNz;
  delete [] values;
  delete [] indices;

  delete &resid;
  delete &z;
  delete &q;
  delete &A;
  delete &map;
  delete &comm;
				       
#ifdef PETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return 0 ;
}
