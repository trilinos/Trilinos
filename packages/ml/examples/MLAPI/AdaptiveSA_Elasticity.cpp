
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
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
// ************************************************************************
//@HEADER

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI.h"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

using namespace Teuchos;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  char filename[200];
  int i = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  int proc  = comm.MyPID();
  int nproc = comm.NumProc();

  // Initialize the workspace and set the output level
  Init();
  
  // create a Epetra_Map from file
  sprintf(filename,"%s/data_update%d.txt",argv[1],nproc);
  Epetra_Map* map = Epetra_ML_readupdatevector(filename,comm);
  if (!map) {
     cout << "**ERR**: could not read map, number of procs ok?\n"; throw -1; }
     
  // read the Epetra_RowMatrix
  sprintf(filename,"%s/data_matrix.txt",argv[1]);
  Epetra_CrsMatrix* Afine = Epetra_ML_readaztecmatrix(filename,*map,comm);
     
     

#ifdef HAVE_MPI
  MPI_Finalize(); 
#endif

  return(0);
}




#endif // #if defined(HAVE_ML_MLAPI)
