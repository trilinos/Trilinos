
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
// Basic definition of communicator.
// This code should be run with at least two processes
// However, not that the SAME code works even if Epetra
// has been configured without MPI!
//
// (output reported at the end of the file)

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  // define an Epetra communicator
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // get the proc ID of this process
  int MyPID = Comm.MyPID();
  
  // get the total number of processes
  int NumProc = Comm.NumProc();
  
  // output some information to std output
  cout << Comm << endl;
  
  // ======================== //
  // now some basic MPI calls //
  // ------------------------ //
  
  int    ivalue, ivalue2, ivalues[NumProc], ivalues2[NumProc];
  double dvalue, dvalue2, dvalues[NumProc], dvalues2[NumProc];
  int root = 0;
  
  // equivalent to MPI_Barrier
  
  Comm.Barrier();
   
  if( MyPID == root ) dvalue = 12.0;

  // on input, the root processor contains the list of values
  // (in this case, a single value). On exit, all processes will
  // have he same list of values. Note that all values must be allocated
  // vefore the broadcast
  
  // equivalent to  MPI_Broadcast
    
  Comm.Broadcast( &dvalue, 1, root );

  // as before, but with integer values. As C++ can bind to the appropriate
  // interface based on argument typing, the type of data is not required.
  
  Comm.Broadcast( &ivalue, 1, root );

  // equivalent MPI_Allgather

  Comm.GatherAll( dvalues, dvalues2, 1);

  // equivalent to MPI_Allreduce with MPI_SUM

  dvalue = 1.0*MyPID;

  Comm.SumAll( &dvalue, dvalues, 1 );

  // equivalent to MPI_Allreduce with MPI_SUM

  Comm.MaxAll( &dvalue, dvalues, 1 );

  // equiavant to MPI_Scan with MPI_SUM

  dvalue = 1.0*MyPID;
  
  Comm.ScanSum( &dvalue, &dvalue2, 1 );

  cout << "On proc " << MyPID << " dvalue2  = " << dvalue2 << endl;
  
  // ======================= //
  // Finalize MPI and return //
  // ----------------------- //
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );
  
} /* main */

/*

Output of this program (NOTE: the output produced by our code can be
slightly different)

[msala:epetra]> mpirun -np 2 ./ex1.exe
Epetra::MpiComm
  Processor 0 of 2 total processors
Epetra::MpiComm
  Processor 1 of 2 total processors
On proc 0 dvalue2  = 0
On proc 1 dvalue2  = 1

*/

