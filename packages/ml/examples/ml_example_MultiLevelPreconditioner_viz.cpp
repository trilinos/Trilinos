
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

// The following header file contains macro definitions for ML. In particular, HAVE_ML_EPETRA,
// HAVE_ML_TEUCHOS, HAVE_ML_TRIUTILS are defines in this file.
#include "ml_config.h"

// the following code cannot be compiled without these Trilinos
// packages. Note that triutils is required in the examples only (to
// generate the linear system), not by the ML library
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS)

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

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Create the linear problem using the class `Trilinos_Util::CrsMatrixGallery.'
  // Various matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  // The matrix here is a VBR matrix, but the code works with any matrix
  // as well. Just note that the coordinates refer to the block rows: if
  // the matrix size if 10000*NumPDEEqns, the coordinate vectors are
  // allocated of size 10000. 
  
  int NumPDEEqns = 5;

  VbrMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 900);

  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetVbrMatrix(NumPDEEqns);
  Epetra_LinearProblem * Problem = Gallery.GetVbrLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== definition of coordinates =================
  
  double * x_coord = 0;
  double * y_coord = 0;
  double * z_coord = 0; // the problem is 2D, here z_coord will be NULL
  
  // use the following triutils matrix gallery function to get the
  // coordinates for a Cartesian grid. Note however that the
  // visualization capabilites of Trilinos accept non-structured grid as
  // well. Visualization and statistics occurs just after the ML
  // preconditioner has been build.

  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // maximum number of levels
  MLList.set("max levels",3);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("smoother: type", "symmetric Gauss-Seidel");

  // aggregation scheme set to METIS. Note that MIS can be
  // visualized only for serial runs, while Uncoupled, METIS and
  // ParMETIS for serial and parallel runs.
  MLList.set("aggregation: type", "METIS");
  MLList.set("aggregation: local aggregates (level 0)", 16);
  MLList.set("aggregation: local aggregates (level 1)", 2);

  // ======================== //
  // visualization parameters //
  // ======================== //
  // 
  // - set "viz: enable" to `false' to disable visualization and
  //   statistics.
  // - set "viz: x-coordinates" to the pointer of x-coor
  // - set "viz: equation to plot" to the number of equation to 
  //   be plotted (for vector problems only). Default is -1 (that is,
  //   plot all the equations)
  // - set "viz: print starting solution" to print on file 
  //   the starting solution vector, that was used for pre-
  //   and post-smoothing, and for the cycle. This may help to
  //   understand whether the smoothed solution is "smooth" 
  //   or not.
  //
  // NOTE: visualization occurs *after* the creation of the ML preconditioner,
  // by calling VisualizeAggregates(), VisualizeSmoothers(), and
  // VisualizeCycle(). However, the user *must* enable visualization 
  // *before* creating the ML object. This is because ML must store some
  // additional information about the aggregates.
  // 
  // NOTE: the options above work only for "viz: output format" == "xyz"
  // (default value). If "viz: output format" == "dx", the user
  // can only plot the aggregates. However, "xyz" works in 2D only.

  MLList.set("viz: enable", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);
  MLList.set("viz: print starting solution", true);

  // =============================== //
  // end of visualization parameters //
  // =============================== //

  // create the preconditioner object and compute hierarchy
  // See comments in "ml_example_epetra_preconditioner.cpp"

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // ============= //
  // visualization //
  // ============= //

  // 1.- print out the shape of the aggregates, plus some
  //     statistics
  // 2.- print out the effect of presmoother and postsmoother
  //     on a random vector. Input integer number represent 
  //     the number of applications of presmoother and postmsoother,
  //     respectively
  // 3.- print out the effect of the ML cycle on a random vector.
  //     The integer parameter represents the number of cycles.

  MLPrec->VisualizeAggregates();
  MLPrec->VisualizeSmoothers(5,1);
  MLPrec->VisualizeCycle(10);

  // ==================== //
  // end of visualization //
  // ==================== //

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  solver.Iterate(500, 1e-12);

  delete MLPrec;
  
  // compute the real residual

  double residual, diff;
  Gallery.ComputeResidualVbr(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutionsVbr(diff);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  // delete memory for coordinates
  if( x_coord ) delete [] x_coord;
  if( y_coord ) delete [] y_coord;
  if( z_coord ) delete [] z_coord;
  
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  return 0 ;
  
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

