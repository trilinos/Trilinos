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
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Expressions.h"

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

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 16);
  Epetra_RowMatrix* EpetraA = Gallery.GetMatrix();

  try {

    // initialize the MLAPI workspace (this has to be done only
    // once, anywhere in the code). Users may want to call
    // MLAPI::Finalize() before quitting their applications.
    MLAPI::Init();

    // All MLAPI objects (that is, DoubleVectors and Operators) requires
    // in input one or more Space's. A space is a light-weight object,
    // which defined the local and global number of elements in the space.
    // A Space object also contains an Epetra_Comm object, so that it can 
    // be used for global reduction operations.
    // Space constructor is very simple, and basically requires the
    // local number of rows. Each row *must* be assigned to exactly one
    // process.
    MLAPI::Space MySpace(EpetraA->NumMyRows());

    // Once a space has been defined, we can construct DoubleVector's as
    // follows:
    MLAPI::DoubleVector x(MySpace), y(MySpace), z(MySpace);

    // Assigning a constant to a vector is as simple as
    x = 1.0;

    // each vector element can be set using the [] operator:
    for (int i = 0 ; i < MySpace.NumMyElements() ; ++i) {
      y(i) = 2.0 * i;
    }

    // DoubleVectors can be added in MATLAB-like notation...
    z = x + y;

#if 0
    // ... also with several vectors:
    z = x - y - z;
#endif

    // to print a vector to standard output, simply type
    cout << z;

    // Operator's are defined as any application that maps from 
    // one space (the so-called DomainSpace) to another space
    // (so-called RangeSpace()). The two spaces can coincide.
    // The Operator constructor accepts any ML_Operator and 
    // Epetra_RowMatrix object. The MLAPI::Operator is a
    // very simple, light-weighted wrapper, and it can be defined as
    MLAPI::Operator A(MySpace,MySpace,*EpetraA);
    // (where the first `MySpace refers to the DomainSpace, and the
    // second to the RangeSpace).

    // We can now start coding the power method. We want a random vector of
    // unitary 2-norm. First, we set random elements in the vector. Then,
    // we divide each element by the 2-norm of the vector, computed as
    // x * x. For two general vectors, x * y represents the scalar product
    // between x and y (defined on the same space).
    x.Random();
    x = x / sqrt(x * x);

    // loop for power method.
    int MaxIters = 10;
    for (int i = 0 ; i < MaxIters ; ++i) {
      // matrix-vector product
      y = A * x;
      // note that you need parenthesis in the following expression!
      cout << "iter = " << i << ", RQ = " << (y * x) / (x * x) << endl;
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
