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
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_MatrixFreePreconditioner.h"
#include "ml_epetra_utils.h"

using namespace ML_Epetra;
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

  int nx;
  if (argc > 1)
    nx = (int) strtol(argv[1],NULL,10);
  else
    nx = 8;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Recirc2D", Map, GaleriList);
  Epetra_VbrMatrix* VbrA = CreateVbrMatrix(CrsA, 1);

  Epetra_MultiVector NullSpace(VbrA->Map(), 1);
  NullSpace.PutScalar(1.0);
  Teuchos::ParameterList MLList;

  // compute the preconditioner using the matrix-free approach
  MatrixFreePreconditioner* MFP = new
    MatrixFreePreconditioner(*VbrA, VbrA->Graph(), MLList, NullSpace);

  // =========================================== //
  // CHECK 1:                                    //
  //                                             //
  // aggregate and build C using ML on VbrA, and //
  // check that the coarse operator is the same  //
  // =========================================== //

  // wrap VbrA as an ML_Operator
  ML_Operator* VbrA_ML = ML_Operator_Create(MFP->Comm_ML());
  ML_Operator_WrapEpetraMatrix(VbrA, VbrA_ML);
  // then build P, R and C using ML, based on VbrA and default null space
  ML_Aggregate* Aggregates_ML;
  ML_Operator* P_ML,* R_ML,* C_ML;
  MFP->Coarsen(VbrA_ML, &Aggregates_ML, &P_ML, &R_ML, &C_ML);

  Epetra_CrsMatrix* C;
  ML_CHK_ERR(ML_Operator2EpetraCrsMatrix(C_ML, C));

  // compare C * vector using ML and matrix-free
  const Epetra_Map& MapC = C->RowMatrixRowMap();
  Epetra_Vector x(MapC), y(MapC), y_MFP(MapC);
  x.Random();
  C->Apply(x, y);
  MFP->C().Apply(x, y_MFP);
  y.Update(1.0, y_MFP, -1.0);
  double norm;
  y.Norm2(&norm);

  if (Comm.MyPID() == 0)
    cout << "norm = " << norm << endl;

  if (norm > 1e-10) exit(EXIT_FAILURE);

#if 0
  Epetra_LinearProblem Problem(A, &LHS, &RHS);
  AztecOO solver(Problem);

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  // destroy the preconditioner
  delete MLPrec;
  
  // compute the real residual

  double residual;
  LHS.Norm2(&residual);
  
  if( Comm.MyPID()==0 ) {
    cout << "||b-Ax||_2 = " << residual << endl;
  }

  // for testing purposes
  if (residual > 1e-5)
    exit(EXIT_FAILURE);

  delete A;
  delete Map;
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
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
  
  return(EXIT_SUCCESS);
}

#endif
