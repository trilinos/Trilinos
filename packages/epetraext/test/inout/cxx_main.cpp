/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "Epetra_ConfigDefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Trilinos_Util.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"

#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  cout << comm << endl;

  int MyPID = comm.MyPID();

  bool verbose = true; 
  if (MyPID==0) verbose = true;

  if(argc < 2 && verbose) {
    cerr << "Usage: " << argv[0] 
	 << " HB_filename" << endl;
    return(1);

  }

  // Uncomment the next three lines to debug in mpi mode
  //int tmp;
  //if (MyPID==0) cin >> tmp;
  //Comm.Barrier();

  Epetra_Map * map;
  Epetra_CrsMatrix * A; 
  Epetra_Vector * x; 
  Epetra_Vector * b;
  Epetra_Vector * xexact;

  int nx = 20*comm.NumProc();
  int ny = 30;
  int npoints = 7;
  int xoff[] = {-1,  0,  1, -1,  0,  1,  0};
  int yoff[] = {-1, -1, -1,  0,  0,  0,  1};

   
  // Call routine to read in HB problem
  Trilinos_Util_GenerateCrsProblem(nx, ny, npoints, xoff, yoff, comm, map, A, x, b, xexact);

  Epetra_Vector tmp1(*map);
  A->Multiply(false, *xexact, tmp1);
  tmp1.Update(1.0, *b, -1.0);
  double residual;
  tmp1.Norm2(&residual);
  if (verbose) cout << "Norm of difference between compute Ax and Ax from file            = " << residual << endl;
  comm.Barrier();

  EpetraExt::RowMatrixToMatrixMarketFile("test.mm", *A, "test matrix", "This is a test matrix");
				       
  delete A;
  delete x;
  delete b;
  delete xexact;
  delete map;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
