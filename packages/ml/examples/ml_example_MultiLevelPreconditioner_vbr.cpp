
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
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_include.h"

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
  // Various matrix examples are supported; please refer to the
  // Trilinos tutorial for more details.
  // This matrix is a VBR matrix, essentially replicating the
  // point-matrix of the example for each unknown. This example is
  // useful to test the vector capabilities of ML, or to debug 
  // code for vector problems. The number of equations is here
  // hardwired as 5, but any positive number (including 1) can be
  // specified.
  //
  // NOTE: ML cannot work with block matrices with variable block size.
  // The number of equations for each block row MUST be constant (as in
  // this example).
  
  int NumPDEEqns = 5;

  // build up a 9-point Laplacian in 2D. This stencil will lead to
  // "perfect" aggregates, of square shape, using almost all the ML
  // aggregation schemes.
  // The problem size (10000) must be a square number. Otherwise, the user
  // can specify the number of points in the x- and y-direction, and the
  // length of the x- and y-side of the computational domain. Please
  // refer to the Trilinos tutorial for more details.
  //
  // Note also that gallery matrices have no boundary nodes.
  
  VbrMatrixGallery Gallery("laplace_2d_9pt", Comm);
  Gallery.Set("problem_size", 10000);
  
  // retrive pointers for linear system matrix and linear problem
  Epetra_RowMatrix * A = Gallery.GetVbrMatrix(NumPDEEqns);
  Epetra_LinearProblem * Problem = Gallery.GetVbrLinearProblem();

  // Construct a solver object for this problem
  AztecOO solver(*Problem);

  // =========================== begin of ML part ===========================
  
  // create a parameter list for ML options
  ParameterList MLList;

  // set defaults for classic smoothed aggregation
  ML_Epetra::SetDefaults("SA",MLList);
  
  // overwrite some parameters. Please refer to the user's guide
  // for more information
  // some of the parameters do not differ from their default value,
  // and they are here reported for the sake of clarity
  
  // maximum number of levels
  MLList.set("max levels",5);
  MLList.set("increasing or decreasing","increasing");

  // set different aggregation schemes for each level. Depending on the
  // size of your problem, the hierarchy will contain different number
  // of levels. As `Uncoupled' and `METIS' are local aggregation
  // schemes, they should be used only for the finest level. `MIS' and
  // `ParMETIS' are global aggregation schemes (meaning that the
  // aggregates can span across processes), and should be reserved for
  // coarsest levels. 
  // Note also that `Uncoupled' and `MIS' will always try to create
  // aggregates of diameter 3 (in the graph sense), while `METIS' and
  // `ParMETIS' can generate aggregates of any size.

  MLList.set("aggregation: type (level 0)", "Uncoupled");
  MLList.set("aggregation: type (level 1)", "MIS");
  MLList.set("aggregation: type (level 2)", "METIS");
  MLList.set("aggregation: type (level 3)", "ParMETIS");
  
  // this is recognized by `METIS' and `ParMETIS' only
  MLList.set("aggregation: nodes per aggregate", 9);
  
  // smoother is Gauss-Seidel. Example file 
  // ml_example_epetra_preconditioner_2level.cpp shows how to use
  // AZTEC's preconditioners as smoothers
  MLList.set("smoother: type","Gauss-Seidel");

  // use both pre and post smoothing. Non-symmetric problems may take
  // advantage of pre-smoothing or post-smoothing only.
  MLList.set("smoother: pre or post", "both");
  
  // solve with serial direct solver KLU
  MLList.set("coarse: type","Amesos_KLU");
  
  // create the preconditioner object and compute hierarchy
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);

  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);

  // =========================== end of ML part =============================
  
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);

  // solve with 500 iterations and 1e-5 tolerance  
  solver.Iterate(500, 1e-5);

  delete MLPrec;
  
  // compute the real residual. Please refer to the Trilinos tutorial
  // for more details. 

  double residual, diff;
  Gallery.ComputeResidualVbr(residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutionsVbr(diff);
  
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
