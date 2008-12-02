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

// Epetra_FEVector Test routine

#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_FEVector.h"
#include "ExecuteTestProblems.h"
#include "../src/Epetra_test_functions.h"
#include "Epetra_Comm.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int ierr = 0;

  //epetra_test::create_comm initializes mpi if EPETRA_MPI is defined
  Epetra_Comm* epetra_comm = epetra_test::create_comm(argc, argv);
  Epetra_Comm& Comm = *epetra_comm;

  int NumProc = Comm.NumProc(); 

  bool verbose =
    epetra_test::global_check_for_flag_on_proc_0("-v",argc,argv, Comm);

  if (verbose) cout << Epetra_Version() << endl << endl;

  int NumVectors = 1;
  int NumMyElements = 4;
  int NumGlobalElements = NumMyElements*NumProc;
  int IndexBase = 0;
  int ElementSize = 1;

  Epetra_BlockMap BlockMap(NumGlobalElements, NumMyElements,
                           ElementSize, IndexBase, Comm);
  BlockMap.SetTracebackMode(0); // This should shut down any error tracing

  EPETRA_TEST_ERR(MultiVectorTests(BlockMap, NumVectors, verbose),ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec0(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec1(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec2(Comm, verbose), ierr);

  Comm.Barrier();
  if (verbose)cout << endl;
  Comm.Barrier();

  EPETRA_TEST_ERR( fevec3(Comm, verbose), ierr);
  delete epetra_comm;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
 
  return ierr;
}

