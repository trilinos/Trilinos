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

// AztecOO_SetParameters Test routine
#include <AztecOO.h>

#ifdef HAVE_AZTECOO_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#endif

#include "Epetra_SerialComm.h"
#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#endif

int main(int argc, char* argv[]) {
  bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;

  int ierr = 0;
  int returnierr = 0;
  int size = 1;
  int rank = 0;

  if (verbose1) {
  }

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  AztecOO azoo;
#ifdef HAVE_AZTECOO_TEUCHOS
  Teuchos::ParameterList paramlist;
  paramlist.set("AZ_solver", AZ_cg);
  int err = azoo.SetParameters(paramlist);
  if (err != 0) {
    if (verbose1) {
      cerr << "err " << err << " returned from AztecOO::SetParameters"<<endl;
    }
    return(-1);
  }

  const int* options = azoo.GetAllAztecOptions();
  if (options[AZ_solver] != AZ_cg) {
    if (verbose1) {
      cerr << "SetParameters test failed."<<endl;
    }
    return(-1);
  }
#endif

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(returnierr);
}

