
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
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
// ***********************************************************************
// @HEADER

// Trilinos Tutorial
// -----------------
// Example of use of Epetra_Time and Epetra_Flops
//
// (output reported at the end of the file)

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;
  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = Map.MyGlobalElements( );

  // ============= CONSTRUCTION OF THE MATRIX ===========================
  // Create a Epetra_Matrix

  Epetra_CrsMatrix A(Copy,Map,3);

  // Add  rows one-at-a-time

  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  int *Indices = new int[2];
  double two = 2.0;
  int NumEntries;

  for( int i=0 ; i<NumMyElements; ++i ) {
    if (MyGlobalElements[i]==0) {
	Indices[0] = 1;
	NumEntries = 1;
    } else if (MyGlobalElements[i] == NumGlobalElements-1) {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
    } else {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      NumEntries = 2;
    }
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, &two, MyGlobalElements+i)==0);
  }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // ================ CONSTRUCTION OF VECTORS =======================
  
  // build up two distributed vectors q and z, and compute
  // q = A * z
  Epetra_Vector q(A.RowMap());
  Epetra_Vector z(A.RowMap());

  // Fill z with 1's
  z.PutScalar( 1.0 );

  // ================ USE OF TIME AND FLOPS =========================
  
  Epetra_Flops counter;
  A.SetFlopCounter(counter);
  Epetra_Time timer(Comm);

  A.Multiply(false, z, q); // Compute q = A*z

  double elapsed_time = timer.ElapsedTime();
  double total_flops =counter.Flops();
  cout << "Total ops: " << total_flops << endl;
  double MFLOPs = total_flops/elapsed_time/1000000.0;
  cout << "Total MFLOPs  for mat-vec = " << MFLOPs << endl<< endl;

  double dotProduct;
  z.SetFlopCounter(counter);
  timer.ResetStartTime();
  z.Dot( q, &dotProduct );

  total_flops =counter.Flops();
  cout << "Total ops: " << total_flops << endl;

  elapsed_time = timer.ElapsedTime();
  MFLOPs = total_flops/elapsed_time/1000000.0;
  cout << "Total MFLOPs for vec-vec = " << MFLOPs << endl<< endl;

  cout << "q dot z = " << dotProduct << endl;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( 0 );

} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex20.exe

Total ops: 26
Total MFLOPs  for mat-vec = 0.242947

Total ops: 26
Total MFLOPs  for mat-vec = 0.0416668

Total ops: 10
Total MFLOPs for vec-vec = 0.0280859

q dot z = 2
Total ops: 10
Total MFLOPs for vec-vec = 0.00914901

q dot z = 2

*/
