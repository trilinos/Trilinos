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
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "BuildTestProblems.h"
#include "ExecuteTestProblems.h"
#include "../epetra_test_err.h"
#include "Epetra_Version.h"

int main(int argc, char *argv[]) {

  int ierr = 0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  Comm.SetTracebackMode(0); // This should shut down any error tracing
  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;

  //  char tmp;
  //  if (rank==0) cout << "Press any key to continue..."<< endl;
  //  if (rank==0) cin >> tmp;
  //  Comm.Barrier();

  Comm.SetTracebackMode(0); // This should shut down any error traceback reporting

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc(); 

  int numMyElements = 3;
  int numGlobalElements = NumProc*numMyElements;

  Epetra_BlockMap pll_map(numGlobalElements, numMyElements, 1, 0, Comm);

  Epetra_MultiVector pll_vec(pll_map, 1);

  double* vals = pll_vec.Values();

  for(int i=0; i<numMyElements; ++i) {
    double val = (MyPID*numMyElements+i)*1.0;
    vals[i] = val;
  }

  std::cout << "pll_vec: " <<std::endl<< pll_vec << std::endl;

 int numGlobal = pll_map.NumGlobalElements();

  int my_pid = pll_map.Comm().MyPID();

  //construct a new Epetra_BlockMap where proc 0
  //has all of the elements in pll_map.

  //numLocal == numGlobal if my_pid == 0,
  //numLocal == 0 if my_pid != 0
  int numLocal = my_pid==0 ? numGlobal : 0;

  Epetra_BlockMap my_map(numGlobal, numLocal, 1, 0, pll_map.Comm());

  //construct an import object
  //with target==my_map, source==pll_map

  Epetra_Import importer(my_map, pll_map);

  //construct the multivector that will hold all elements
  //on proc 0
  Epetra_MultiVector my_vec(my_map, pll_vec.NumVectors());

  //now import pll_vec's contents into my_vec
  my_vec.Import(pll_vec, importer, Insert);

  std::cout << "my_vec: " << std::endl << my_vec << std::endl;

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ierr;
}

