
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

// \author Marzio Sala, SNL 9214
// \data Last modified on 05-Dec-04

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), required Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-triutils (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

// epetra objects
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
// required to build the example matrix
#include "Trilinos_Util_CrsMatrixGallery.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_user.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// ============================================= //
// defines the label for this aggregation scheme //
// ============================================= //

char*UserLabel()
{
  return("Linear");
}

// ================================ //
// defines the linear decomposition //
// ================================ //
//
// Input:
// - Amatrix:           ML_Operator containing the matrix whose
//                      local graph has to be partitioned. This matrix
//                      has already been "amalgamated", so that all
//                      equations of a given node are glued together.
//                      The values of Amatrix are defined as the Frobenius
//                      norm of each block. All matrix entries are
//                      positive. Amatrix contains local rows and
//                      columns, *and* columns corresponding to external
//                      (ghost) nodes, in ML ordering (which may differ
//                      from the user's ordering also for the finest-level
//                      matrix).
//                      [input only]
//
// - bdry_nodes:        array of character, of size Amatrix->invec_leng.
//                      The user must fill this array with 'T' and 'F'
//                      ('T' = is a boundary node, 'F' is not a boundary
//                      node).
//                      [output only]
//
// - epsilon:           dropping parameter
//                      [input only]
//
// - x, y, z            double arrays, containing the coordinates of 
//                      each block node. If the problem is 2D, z is 0.
//                      If coordinates have not been given to ML,
//                      x, y and z are all 0.
//                      [input only]
//
// - partitions         integer array, of size Amatrix->invec_leng.
//                      The user must fill this array with the *local*
//                      aggregate number of each node. For example,
//                      if partitions[10] = 12, this means that the
//                      local node 10 belongs to the local aggregate 12.
//                      [output only]
//
// - Nonzeros           returns the number of local nonzeros. This is
//                      required only to compute the complexity of the
//                      operator. It can be set to 0.
//                      [output only]
//
// return value is the number of local aggregates.

int UserPartitions(ML_Operator* Amatrix, char* bdry_nodes,
                   double epsilon,
                   double* x,double* y,double* z,
                   int* partitions, int* Nonzeros)
{
  int ChunkSize = 10;

  // get the number of local rows in Amatrix
  int NumLocalRows = Amatrix->invec_leng;
  
  // set the number of local aggregates
  int NumAggregates = NumLocalRows / ChunkSize;
  if (NumLocalRows % ChunkSize)
    ++NumAggregates;

  // fill `partitions'
  for (int i = 0 ; i < NumLocalRows ; ++i) {
    partitions[i] = i / ChunkSize;
  }

  // detect the boundary nodes (i.e., rows with no off-diagonal elements)
  // This requires calling ML's getrow(), as detailed in points 
  // 1, 2 and 3.

  // 1.- get the pointer to the getrow() function from the
  //     ML_Operator input object
  int (*getrow)(ML_Operator *,int,int*,int,int*,double*,int*) =
    Amatrix->getrow->func_ptr;

  // 2.- allocate tentative space for getrow(), in this case
  //     128, but any nonzero number (for example, 1) is ok.
  //     Also, allocate space to hold the indices of the nonzeros
  //     in each row, and their values.
  int MaxRowEntries = 128;
  vector<int> Indices(MaxRowEntries);
  vector<double> Values(MaxRowEntries);
  
  // 3.- loop over all local rows
  for (int i = 0 ; i < NumLocalRows ; ++i) {
    int RowEntries;
    // if return value is zero, then MaxRowEntries is too small
    // to copy the nonzeros in the provided array. So, we
    // reallocate them, and recall the function
    while(getrow(Amatrix,1,&i,MaxRowEntries, 
                      &Indices[0],&Values[0],&RowEntries) == 0) {
      MaxRowEntries *= 2;
      Indices.resize(MaxRowEntries);
      Values.resize(MaxRowEntries);
    }

    /* prints out matrix entries in local numbering
    for (int j = 0 ; j < RowEntries ; ++j)
      cout << "A(" << i << "," << Indices[j] << ") = " << Values[j] << endl;
    */

    // if only one nonzero is in this row, then we skip
    // this singleton
    if (RowEntries <= 1) {
      bdry_nodes[i] = 'T';
    } else {
      bdry_nodes[i] = 'F';
      *Nonzeros += RowEntries;
    }
  }
  
  // return the number of local aggregates
  return(NumAggregates);
}

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

  // Create the linear problem using the class `Trilinos_Util::CrsMatrixGallery.'
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 100);
  
  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  double * x_coord = 0;
  double * y_coord = 0;
  double * z_coord = 0; // the problem is 2D, here z_coord will be NULL
  
  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);
  int NumMyRows = A->NumMyRows();

  // As we wish to use AztecOO, we need to construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation.
  ML_Epetra::SetDefaults("SA",MLList);
  
  // use user's defined aggregation scheme to create the aggregates
  // 1.- set "user" as aggregation scheme (for all levels, or for
  //     a specify level only)
  MLList.set("aggregation: type", "user");
  // 2.- set the label (for output)
  ML_SetUserLabel(UserLabel);
  // 3.- set the aggregation scheme (see function above)
  ML_SetUserPartitions(UserPartitions);
  // 4.- set the coordinates.
  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("aggregation: dimensions", 2);

  // also setup some variables to visualize the aggregates
  // (more details are reported in example `ml_viz.cpp'.
  MLList.set("viz: enable", true);

  // now we create the preconditioner
  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  MLPrec->VisualizeAggregates();

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(500, 1e-12);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  if (residual > 1e-5)
    exit(EXIT_FAILURE);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

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
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos");
  puts("--enable-aztecoo --enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) */
