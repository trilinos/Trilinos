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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int main(int argc, char* argv[]) {
  bool verbose = false;  // used to set verbose false on non-root processors
  bool verbose1 = false; // user's command-line argument
  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose1 = true;

  int err;
  int returnierr = 0;
  int size = 1;
  int rank = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  if (verbose1) {
  }

#ifdef HAVE_AZTECOO_TEUCHOS
  AztecOO azoo;

  Teuchos::ParameterList paramlist;
  paramlist.set("AZ_solver", AZ_cg);
  paramlist.set("max_ITER", 200);
  paramlist.set("az_Tol", 1.0e-09);
  paramlist.set("precond", AZ_Jacobi);

  paramlist.set("AZ_kspace", 2.5);//bad type

  if (verbose1==true) {
    cout << "Test: parameter 'AZ_kspace' given bad type (double), warning should"
         << "be printed to cerr"<<endl;
    err = azoo.SetParameters(paramlist, verbose1);
  }
  else {
    err = azoo.SetParameters(paramlist);
  }
  if (err != 0) {
    if (verbose1) {
      cerr << "err " << err << " returned from AztecOO::SetParameters"<<endl;
    }
    return(-1);
  }

  const int* options = azoo.GetAllAztecOptions();
  if (options[AZ_solver] != AZ_cg) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_solver."<<endl;
    }
    return(-1);
  }

  if (options[AZ_max_iter] != 200) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_max_iter."<<endl;
    }
    return(-1);
  }

  if (options[AZ_precond] != AZ_Jacobi) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_precond."<<endl;
    }
    return(-1);
  }

  const double* params = azoo.GetAllAztecParams();
  if (params[AZ_tol] != 1.0e-09) {
    if (verbose1) {
      cerr << "SetParameters test failed to correctly set AZ_tol."<<endl;
    }
    return(-1);
  }
#endif

  if (verbose1==true) {
    cout << "********* Test passed **********" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(returnierr);
}

