
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
  // `ml_example_epetra_preconditioner_vbr.cpp' shows how to define a
  // Epetra_VbrMatrix.
  
  // `laplace_2d' is a symmetric matrix; an example of non-symmetric
  // matrices is `recirc_2d' (advection-diffusion in a box, with
  // recirculating flow). The number of nodes must be a square number

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  int ProblemSize = 256;
  Gallery.Set("problem_size", ProblemSize);
  
  // The following methods of CrsMatrixGallery are used to get pointers
  // to internally stored Epetra_RowMatrix and Epetra_LinearProblem.

  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  // As we wish to use AztecOO, we need to construct a solver object for this problem
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

  // Here we set parameters to visualiza the effect of the actual smoothers
  // and the ML cycle on a random vector.
  //
  // First, we get the nodal coordinates. 
  // NOTE: This example can work with VBR matrices as well.
  // NOTE 2: memory for x_coord, y_coord and z_coord is allocated using
  // `new' in GetCartesianCoordinates(). (Actually, z_coord is not allocated,
  // as the problem is 2D.)

  double * x_coord = 0;
  double * y_coord = 0;
  double * z_coord = 0; // the problem is 2D, here z_coord will be 0
  
  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  // set parameters for visualization
  MLList.set("viz: enable", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);

  // create the preconditioning object.
  // NOTE that users need to set "viz: enable" == true in order to
  // visualize!

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  MLPrec->PrintStencil2D(16,16);
  exit(0);

  // =============== //
  // A N A L Y S I S //
  // =============== //

  if (ProblemSize < 1024 && Comm.NumProc() == 1) {

    // Analysis with "dense" involves dense eigenvalue problems.
    // This can be very expensive. Also, it can be done only for serial
    // computations.
 
    MLPrec->AnalyzeSmoothersDense(5,5);
    MLPrec->AnalyzeMatrixEigenvaluesDense("A");
    MLPrec->AnalyzeMatrixEigenvaluesDense("P^{-1}A");
    
  }
  else {
    
    // on the other hand, "sparse" analysis can be applied to serial and
    // parallel, of any size, but we cannot get the entire spectrum of the
    // operators.
    // NOTE: the eigencomputation can be expensive!
 
    MLPrec->AnalyzeSmoothersSparse(5,5);
    MLPrec->AnalyzeMatrixEigenvaluesSparse("A");
    MLPrec->AnalyzeMatrixEigenvaluesSparse("P^{-1}A");


  }

  // still to perform:
  // 1.- a "cheap" analysis of the matrix (mainly, whether it is
  //     diagonally domimant and well scaled)
  // 2.- analyze the effect of the ML cycle on a random vector

  MLPrec->AnalyzeMatrixCheap();
  MLPrec->AnalyzeCycle(5);

  // Here we set parameters to compare different smoothers.
  // Please refer to the user's guide for more details about the following
  // parameters.

  Teuchos::ParameterList MLTestList;
  ML_Epetra::SetDefaults("SA",MLTestList);
  
  MLTestList.set("test: Jacobi", true);
  MLTestList.set("test: Gauss-Seidel", true);
  MLTestList.set("test: symmetric Gauss-Seidel", true);
  MLTestList.set("test: Aztec", true);
  MLTestList.set("test: Aztec solver", false);

  MLTestList.set("test: sweeps", 5);

  MLPrec->TestSmoothers(MLTestList);

  // =============== //
  // end of analysis //
  // =============== //

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-12 tolerance  
  // The problem should converge as follows:
  //
  // proc       iterations       condition number
  //   1             14               1.78
  //   2             15               2.39
  //   4             15               2.20

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
