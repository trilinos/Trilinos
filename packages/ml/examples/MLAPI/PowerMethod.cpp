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

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_RowMatrix.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "ml_include.h"
#include "MLAPI.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  try {

    // must be a square number
    int NumGlobalRows = 16;

    // initialize the MLAPI workspace
    MLAPI::Init();

    // Create the space in which vectors and operators live
    MLAPI::Space MySpace(NumGlobalRows);

    // Define vectors and operator
    MLAPI::MultiVector x(MySpace), y(MySpace);
    MLAPI::Operator A = Gallery("laplace_1d", MySpace);

    // We can now start coding the power method. We want a random vector of
    // unitary 2-norm. First, we set random elements in the vector. Then,
    // we divide each element by the 2-norm of the vector, computed as
    // x * x. For two general vectors, x * y represents the scalar product
    // between x and y (defined on the same space).
    x.Random();
    x = x / sqrt(x * x);

    // loop for power method.
    int MaxIters = 20;
    for (int i = 0 ; i < MaxIters ; ++i) {
      // matrix-vector product
      y = A * x;
      // note that you need parenthesis in the following expression
      // and that x * y is a global operation
      double RQ = (y * x) / (x * x);
      if (GetMyPID() == 0)
        cout << "iter = " << i << ", RQ = " << RQ << endl;
      x = y / sqrt(y * y);
    }

    MLAPI::Finalize();
  } 
  catch (exception& e) {
    cout << e.what() << endl;
  } 
  catch (...) {
    cout << "problems here..." << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);

}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");

  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
