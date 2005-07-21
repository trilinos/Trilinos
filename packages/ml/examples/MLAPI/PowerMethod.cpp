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
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)

#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"

using namespace MLAPI;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  try {

    // initialize the MLAPI workspace
    
    Init();

    // Global dimension of our problem. It must be a square number
    
    int NumGlobalRows = 16;

    // Create the space in which vectors and operators live
    
    Space MySpace(NumGlobalRows);

    // Define vectors and operator, based on `MySpace'. The linear system
    // matrix is created using the MLAPI Gallery, which is only a 
    // convenient wrapper of the Triutils gallery.
    
    MultiVector x(MySpace), y(MySpace);
    Operator A = Gallery("laplace_1d", MySpace);

    // We can now start coding the power method. We want a random vector of
    // unitary 2-norm. First, we set random elements in the vector. Then,
    // we divide each element by the 2-norm of the vector, computed as
    // x * x. For two general vectors, x * y represents the scalar product
    // between x and y (defined on the same space).
    
    x.Random();
    x = x / sqrt(x * x);

    // ==================== //
    // loop of power method //
    // ==================== //
    
    int MaxIters = 10;

    for (int i = 0 ; i < MaxIters ; ++i) {

      y = A * x;  // matrix-vector product

      // note that you need parenthesis in the following expression
      // and that x * y is a global operation (therefore all
      // processes must execute it)
      
      double RQ = (y * x) / (x * x);
      if (GetMyPID() == 0)
        cout << "iter = " << i << ", RQ = " << RQ << endl;

      x = y / sqrt(y * y);

    }

    // finalize the MLAPI workspace

    Finalize();

  } 
  catch (const int e) {
    cout << "Caught integer exception, code = " << e << endl;
  } 
  catch (...) {
    cout << "problems here..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("The ML API requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif // #if defined(HAVE_ML_MLAPI)
