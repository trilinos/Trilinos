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

// To be modified to support Galeri and no longer Triutils
// I am not sure when this test was executed and passed... MS
#ifdef FIXME
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_ZOLTAN)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace ML_Epetra;

// *) define VIZ to visualize the aggregates (does not work for
//    all the aggregation schemes)
//
// \author Marzio Sala, SNL 9214
//
// \date last updated on 26-Mar-05

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int NumNodes = 65536;
  int NumPDEEqns = 2;

  VbrMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumNodes);
  Gallery.Set("output",0);

  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix* A = Gallery.GetVbrMatrix(NumPDEEqns);
  Epetra_LinearProblem* Problem = Gallery.GetVbrLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  double* x_coord = 0;
  double* y_coord = 0;
  double* z_coord = 0; // the problem is 2D, here z_coord will be NULL
  
  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  ParameterList MLList;
  MLList.set("ML output",8);

  MLList.set("max levels",10);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  // *) if a low number, it will use all the available processes
  // *) if a big number, it will use only processor 0 on the next level
  MLList.set("aggregation: next-level aggregates per process", 1);

  MLList.set("aggregation: type (level 0)", "Zoltan");
  MLList.set("aggregation: type (level 1)", "Uncoupled");
  MLList.set("aggregation: type (level 2)", "Zoltan");

  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", z_coord);

  // specify the reduction with respect to the previous level
  // (very small values can break the code)
  int ratio = 16;
  MLList.set("aggregation: global aggregates (level 0)", 
             NumNodes / ratio);
  MLList.set("aggregation: global aggregates (level 1)", 
             NumNodes / (ratio * ratio));
  MLList.set("aggregation: global aggregates (level 2)", 
             NumNodes / (ratio * ratio * ratio));

#ifdef VIZ
  MLList.set("viz: enable", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);
#endif

  // create the preconditioner object and compute hierarchy
  // See comments in "ml_example_epetra_preconditioner.cpp"

  MultiLevelPreconditioner* MLPrec = 
    new MultiLevelPreconditioner(*A, MLList, true);

#ifdef VIZ  
  MLPrec->VisualizeAggregates();
#endif

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidualVbr(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutionsVbr(&diff);
  
  if (Comm.MyPID() == 0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // delete memory for coordinates
  if (x_coord) delete[] x_coord;
  if (y_coord) delete[] y_coord;
  if (z_coord) delete[] z_coord;
  
  if (Comm.MyPID() && residual > 1e-5) {
    cout << "TEST FAILED!!!!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  if (Comm.MyPID() == 0)
    cout << "TEST PASSED" << endl;

  exit(EXIT_SUCCESS);
  
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  /* still need to deal with MPI, some architecture don't like */
  /* an exit(0) without MPI_Finalize() */
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils --with-ml_zoltan");
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_ZOLTAN) */

#endif
