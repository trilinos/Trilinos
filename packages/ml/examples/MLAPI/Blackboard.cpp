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
#include "ml_include.h"
#include "MLAPI.h"

using namespace Teuchos;
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

    // initialize the MLAPI workspace (this has to be done only
    // once, anywhere in the code). To avoid memory leaks, users 
    // may want to call MLAPI::Finalize() before quitting their applications.
    Init();

#if 0
    // Declaration of some empty variables
    Operator A, B, C;
    InverseOperator invA;
#endif
    // define the space in which operators and vectors will live
    int NumGlobalElements = 16;
    Space MySpace(NumGlobalElements);

    DoubleVector x(MySpace), y(MySpace), z(MySpace);

    x = 1.0;
    y = 2.0;
    z = 3.0;
    cout << x + x;
    cout << x;

#if 0
    x = 1.0;
    for (int i = 0 ; i < y.MyLength() ; ++i)
      y(i) = 1.0 * i;

    z = x + y;

    // Create a square matrix corresponding to a 2D Laplacian
    A = Gallery("laplace_2d", MySpace);

    // B will be a diagonal matrix so that B_{i,i} = A_{i,i}
    B = Diagonal(Diagonal(A));

    // verify that the diagonal of A - B is zero
    z = Diagonal(B - A);

    // to print a vector to standard output, simply type
    if (MyPID() == 0)
      cout << z;

    // C is the product between A and B (works for *any* compatible
    // MLAPI::Operators)
    C = A * B;

    // or the sum between A and B
    C = A + B;

    // or the difference
    C = A - B;

    // set to random valules
    x.Random();

    // scale x by the dot product between x and y
    x = x / (x * y);

    // use Amesos to apply the inverse of A using LU factorization
    ParameterList MLAPIList;
    invA.Reshape(A,"Amesos", MLAPIList);

    // verify that x == inv(A) * A * x
    x = invA * (A * x) - x;
    double NormX = sqrt(x * x);
    if (MyPID() == 0)
      cout << "Norm of inv(A) * A * x - x = " << NormX << endl;

#endif

    // Finalize the MLAPI work space before leaving the application
    Finalize();
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
