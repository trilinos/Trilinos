
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

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif
#include "ml_config.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), required Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-triutils (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
// includes required by ML

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "MLAPI_Space.h"
#include "MLAPI_Vector.h"
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
  // once, anywhere in the code)
  MLAPI::Init();

  // define the space on this the ML matrix and vectors will live
  MLAPI::Space MySpace(EpetraA->NumMyRows(),Comm);

  // here we define three vectors, and we compute their sum
  // (store in one of the three vectors)
  MLAPI::Vector v1(MySpace);
  v1 = 1.0;
  MLAPI::Vector v2(MySpace);
  v2 = 2.0;
  MLAPI::Vector v3(MySpace);
  v3 = 3.0;

  v3 = v1 + v2 + v3;
  cout << v3;

  // now we wrap (light-weight conversion) the Epetra matrix into
  // an MLAPI::Operator
  MLAPI::Operator A(MySpace,MySpace,*EpetraA);

  MLAPI::Vector x(MySpace), y(MySpace);

  // define a vector of norm 1
  x = 1.0;
  x = x / sqrt(x * x);

  // simple power method applied to A. Note that (x * x) requires
  // parenthesis.
  // x * y == scalar product between x and y
  for (int i = 0 ; i < 10 ; ++i) {
    double alpha;
    y = A * x;
    cout << "iter = " << i << ", RQ = " << (y * x) / (x * x) << endl;
    x = y / sqrt(y * y);
  }

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
