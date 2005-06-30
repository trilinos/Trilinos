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

// Epetra_BlockMap Test routine

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {
  bool verbose = false;
  // Check if we should print results to standard out
  if (argc > 1) {
    if ((argv[1][0] == '-') && (argv[1][1] == 'v')) {
      verbose = true;
    }
  }

  int i, ierr = 0, returnierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (!verbose) {
    Comm.SetTracebackMode(0); // This should shut down any error traceback reporting
  }
  int MyPID = Comm.MyPID();
  int NumProcs = Comm.NumProc();

  int verbose_int = verbose ? 1 : 0;
  Comm.Broadcast(&verbose_int, 1, 0);
  verbose = verbose_int==1 ? true : false;

  if (verbose && MyPID==0)
    cout << Epetra_Version() << endl << endl;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}
