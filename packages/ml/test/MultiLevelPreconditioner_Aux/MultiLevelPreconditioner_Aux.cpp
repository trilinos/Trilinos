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

#include "ml_include.h"

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
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  VbrMatrixGallery Gallery("stretched_2d", Comm);
  //VbrMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", Comm.NumProc() * 900);
  Gallery.Set("map_type", "box");
  int NumPDEs = 5;
  Epetra_RowMatrix* A = Gallery.GetVbrMatrix(NumPDEs);
  Epetra_LinearProblem* Problem = Gallery.GetVbrLinearProblem();

  AztecOO solver(*Problem);

  ParameterList MLList;

  ML_Epetra::SetDefaults("DD-ML",MLList);
  
  double* x_coord = 0;
  double* y_coord = 0;
  double* z_coord = 0; 

  Gallery.GetCartesianCoordinates(x_coord, y_coord, z_coord);

  for (int i = 0 ; i < A->NumMyRows() / NumPDEs ; ++i) 
  {
    x_coord[i] *= 100;
  }

  int MaxLevels = 10;

  if (0)
  {
    // this is the old stuff
    MLList.set("x-coordinates", x_coord);
    MLList.set("y-coordinates", y_coord);
    MLList.set("z-coordinates", z_coord);
    MLList.set("aggregation: use auxiliary matrix", true);
    MLList.set("aggregation: threshold", 0.05);
  }
  
  if (1)
  {
    // this is the new stuff
    MLList.set("aggregation: threshold", 0.05);
    MLList.set("aggregation: aux: enable", true);
    MLList.set("aggregation: aux: threshold", 0.05);
    MLList.set("x-coordinates", x_coord);
    MLList.set("y-coordinates", y_coord);
    MLList.set("z-coordinates", z_coord);
    MLList.set("aggregation: aux: max levels", MaxLevels);
  }

  MLList.set("aggregation: damping factor", 0.0);

  MLList.set("output", 10);
  MLList.set("max levels",MaxLevels);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: type","symmetric Gauss-Seidel");
  MLList.set("smoother: pre or post", "both");
  //MLList.set("low memory usage", true);

#ifdef HAVE_ML_AMESOS
  MLList.set("coarse: type","Amesos-KLU");
#else
  MLList.set("smoother: type","Jacobi");
  MLList.set("coarse: type","Jacobi");
#endif

#ifdef IFPACK_SMOOTHER
  if (0)
  {
    MLList.set("smoother: type","IFPACK");
    Teuchos::ParameterList& IFPACKList = MLList.sublist("smoother: ifpack list");
    MLList.set("smoother: ifpack type", "block relaxation stand-alone");
    MLList.set("smoother: ifpack overlap", 0);
    IFPACKList.set("relaxation: zero starting solution",false);
    IFPACKList.set("relaxation: type", "Gauss-Seidel");
    IFPACKList.set("relaxation: sweeps", 1);
    IFPACKList.set("relaxation: damping factor", 0.67);
    IFPACKList.set("partitioner: type", "user");
  }
#endif

#ifdef VIZ_ME
  MLList.set("viz: output format", "xyz");
  MLList.set("viz: enable", true);
  MLList.set("viz: x-coordinates", x_coord);
  MLList.set("viz: y-coordinates", y_coord);
  MLList.set("viz: z-coordinates", z_coord);
  MLList.set("viz: print starting solution", true);
#endif

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // be sure that this is still ok
  MLPrec->ReComputePreconditioner();

#ifdef VIZ_ME
  for (int i = 0 ; i < A->NumMyRows() / NumPDEs ; ++i) {
    x_coord[i] /= 100;
  }
  MLPrec->VisualizeAggregates();
#endif

  solver.GetRHS()->PutScalar(0.0);
  solver.GetLHS()->Random();

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  delete MLPrec;
  
  double norm;
  solver.GetLHS()->Norm2(&norm);

  if (Comm.MyPID() == 0)
    cout << "Error = " << norm << endl;

  if (norm > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
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

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-triutils");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif 
