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

// Goal of this example is to show the (limited) ML analysis capabilities.
//
// \author Marzio Sala, SNL 9214
// \date Last modified on 17-Nov-04

#include "ml_include.h"

// The C++ interface of ML (more precisely,
// ML_Epetra::MultiLevelPreconditioner), requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// required --enable-triutils (for the definition of the linear systems)

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

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
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

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
  // Several matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  // Most of the examples using the ML_Epetra::MultiLevelPreconditioner
  // class are based on Epetra_CrsMatrix. Example
  // `ml_EpetraVbr.cpp' shows how to define a Epetra_VbrMatrix.
  
  // `laplace_2d' is a symmetric matrix; an example of non-symmetric
  // matrices is `recirc_2d' (advection-diffusion in a box, with
  // recirculating flow). The number of nodes must be a square number

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  int ProblemSize = 256;
  Gallery.Set("problem_size", ProblemSize);
  
  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix* A = Gallery.GetMatrix();
  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  // As we wish to use AztecOO, we need to construct a solver object 
  // for this problem
  AztecOO solver(*Problem);

  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  // use Uncoupled scheme to create the aggregate
  MLList.set("aggregation: type", "Uncoupled");
  
  // fix the smoother
  MLList.set("smoother: type","symmetric Gauss-Seidel");

  // =================================== //
  // V I S U A L I Z A T I O N   P A R T //
  // =================================== //

  // Here we set parameters to visualize the effect of the actual smoothers
  // and the ML cycle on a random vector.
  //
  // First, we get the nodal coordinates. 
  // NOTE: This example can work with VBR matrices as well.
  // NOTE 2: memory for x_coord, y_coord and z_coord is allocated using
  // `new' in GetCartesianCoordinates(). (Actually, z_coord is not allocated,
  // as the problem is 2D.)

  double* x_coord = 0;
  double* y_coord = 0;
  double* z_coord = 0; // the problem is 2D, here z_coord will be 0
  
  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  // set parameters for visualization
  
  MLList.set("viz: enable", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);

  // create the preconditioning object.
  // Note that users need to set "viz: enable" == true in order to
  // visualize!

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // for 2D Cartesian grid, you can print the stencil of your operator
  // using this simple function.
  
  MLPrec->PrintStencil2D(16,16);

  // =========================== //
  // E I G E N - A N A L Y S I S //
  // =========================== //

  if (ProblemSize < 1024 && Comm.NumProc() == 1) {

    // Analysis with "dense" involves dense eigenvalue problems.
    // This can be very expensive, as the sparse matrix on each level
    // is converted to a dense one, then LAPACK routines are used
    // Also, it can be done only for serial computations.
    // `5,5' refers to the number of pre-smoother and post-smoother
    // applications to a random vector.
 
    MLPrec->AnalyzeSmoothersDense(5,5);
    MLPrec->AnalyzeMatrixEigenvaluesDense("A");
    MLPrec->AnalyzeMatrixEigenvaluesDense("P^{-1}A");
    
  }

  // ================================================= //
  // A N A L Y S I S   O F   T H E   H I E R A R C H Y //
  // ================================================= //

  // Method AnalyzeHierarchy() can be used to validate an
  // already built hierarchy.
  // - `true' means perform a "cheap" analysis of each level's matrix
  // - Then, each level's smoothers and the complete cycle are 
  //   applied to solve the problem
  //     A e = 0
  //   with a random initial solution, to get a sense of the effectiveness
  //   of the smoothers and the cycle itself. The parameters are:
  //   * NumPreCycles and NumPostCycles specify the number of post
  //     and pre smoother applications;
  //   * NumMLCycles specifies the number of applications of the 
  //     complete cycle.

  MLPrec->AnalyzeHierarchy(true, NumPreCycles, NumPostCycles, NumMLCycles);

  // ================================================= //
  // A N A L Y S I S   O F   T H E   S M O O T H E R S //
  // ================================================= //

  // Method TestSmoothers() can be used to analyze different smoothers
  // on a given problem. The cycle is built following the parameters
  // specified in MLTestList.
  // Please refer to the user's guide for more details about the following
  // parameters. Not all smoothers are supported by this testing.

  Teuchos::ParameterList MLTestList;
  ML_Epetra::SetDefaults("SA",MLTestList);
  
  MLPrec->TestSmoothers(MLTestList);

  // =============== //
  // end of analysis //
  // =============== //

  delete MLPrec;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
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
