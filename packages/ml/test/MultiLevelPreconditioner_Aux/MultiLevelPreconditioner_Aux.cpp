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

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Galeri_CrsMatrices.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_Utils.h"

using namespace Teuchos;
using namespace Galeri;

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

  ParameterList GaleriList;
  GaleriList.set("nx", 10);
  GaleriList.set("ny", 10 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1 * Comm.NumProc());
  GaleriList.set("lx", 1.0);
  GaleriList.set("ly", 1.0 * Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Stretched2D", Map, GaleriList);

  int NumPDEs = 5;
  Epetra_VbrMatrix* VbrA = CreateVbrMatrix(CrsA, NumPDEs);

  int NumVectors = 1;
  Epetra_MultiVector LHS(VbrA->Map(), NumVectors);
  Epetra_MultiVector RHS(VbrA->Map(), NumVectors);
  Epetra_LinearProblem Problem(VbrA, &LHS, &RHS);

  Epetra_MultiVector* Coord = CreateCartesianCoordinates("2D", &(CrsA->Map()), 
                                                         GaleriList);
                                                         
  double* x_coord = (*Coord)[0];
  double* y_coord = (*Coord)[1];

  for (int i = 0 ; i < CrsA->NumMyRows() ; ++i) 
  {
    x_coord[i] *= 100;
  }

  AztecOO solver(Problem);

  ParameterList MLList;

  ML_Epetra::SetDefaults("DD-ML",MLList);

  int MaxLevels = 10;

  // this is the new stuff
  MLList.set("aggregation: aux: enable", true);
  MLList.set("aggregation: aux: threshold", 0.05);
  MLList.set("x-coordinates", x_coord);
  MLList.set("y-coordinates", y_coord);
  MLList.set("z-coordinates", (double*)0);
  MLList.set("aggregation: aux: max levels", MaxLevels);

  MLList.set("coarse: max size", 1024);
  MLList.set("aggregation: damping factor", 0.0);

  MLList.set("output", 10);
  MLList.set("max levels",MaxLevels);
  MLList.set("increasing or decreasing","increasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: type","symmetric Gauss-Seidel");
  MLList.set("smoother: pre or post", "both");
  MLList.set("low memory usage", true);

#ifdef HAVE_ML_AMESOS
  MLList.set("coarse: type","Amesos-KLU");
#else
  MLList.set("smoother: type","Gauss-Seidel");
  MLList.set("coarse: type","Gauss-Seidel");
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
  MLList.set("viz: z-coordinates", (double*)0);
  MLList.set("viz: print starting solution", true);
#endif

  ML_Epetra::MultiLevelPreconditioner * MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*VbrA, MLList);

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

  double norm[NumVectors];
  solver.GetLHS()->Norm2(norm);

  if (Comm.MyPID() == 0)
    for (int i=0; i<NumVectors; i++)
      printf("Error[%d] = %e\n",i,norm[i]);

  double maxNorm=0.0;
  for (int i=0; i<NumVectors; i++)
    if (norm[i] > maxNorm) maxNorm = norm[i];
    
  if (maxNorm > 1e-5) {
    if (Comm.MyPID() == 0) printf("Test failed.\n"); 
    exit(EXIT_FAILURE);
  }

  delete Coord;
  delete VbrA;
  delete CrsA;
  delete Map;

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
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
}

#endif 
