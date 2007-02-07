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
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_MPI)

#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "AztecOO.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Amesos.h"
#include "ml_Ifpack_ML.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

// ====================================================================== 

template<class T>
int TestAdditiveSchwarz()
{
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  ParameterList GaleriList;
  GaleriList.set("nx", 8);
  GaleriList.set("ny", 8 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorRangeMap());

  LHS.PutScalar(0.0);
  RHS.PutScalar(1.0);

  Epetra_LinearProblem problem(A, &LHS, &RHS);

  AztecOO solver(problem);

  Teuchos::ParameterList List;
  ML_Epetra::SetDefaults("SA", List);
  List.set("ML output", 0);
  List.set("schwarz: combine mode", Add);
  List.set("cycle applications", 10);

  int Overlap = 0;
  Ifpack_AdditiveSchwarz<T> Prec(A, Overlap);
  Prec.SetParameters(List);
  ML_CHK_ERR(Prec.Initialize());
  ML_CHK_ERR(Prec.Compute());

  solver.SetPrecOperator(&Prec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  delete A;
  delete Map;

  return(solver.NumIters());
}

// ====================================================================== 
int TestML(const bool UseSelf) 
{
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  ParameterList GaleriList;
  GaleriList.set("nx", 8);
  GaleriList.set("ny", 8 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);

  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorRangeMap());

  LHS.PutScalar(0.0);
  RHS.PutScalar(1.0);

  Epetra_LinearProblem problem(A, &LHS, &RHS);

  AztecOO solver(problem);

  ParameterList MLList;
  ML_Epetra::SetDefaults("SA",MLList);

  MLList.set("increasing or decreasing", "increasing");
  MLList.set("max levels",5);
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: pre or post","pre");
  MLList.set("ML output", 0);

  // toggle the following to compare with more classical smoothing

  if (UseSelf) 
  {
    MLList.set("smoother: type (level 0)", "Gauss-Seidel");
    MLList.set("smoother: type (level 1)", "self");
    MLList.set("smoother: type (level 2)", "self");
    MLList.set("smoother: type (level 3)", "self");
    MLList.set("smoother: type (level 4)", "self");

    Teuchos::ParameterList& SelfList = MLList.sublist("smoother: self list");
    ML_Epetra::SetDefaults("DD-ML", SelfList);
    SelfList.set("ML output", 0);
    SelfList.set("coarse: max size", 128);
    SelfList.set("cycle applications", 10);
    SelfList.set("aggregation: damping factor", 0.0);
    SelfList.set("smoother: pre or post", "both");
    SelfList.set("max levels", 5);
    SelfList.set("zero starting solution", true);
  } 
  else 
  {
    MLList.set("smoother: type (level 0)", "Gauss-Seidel");
    MLList.set("smoother: type (level 1)", "IFPACK");
    MLList.set("smoother: type (level 2)", "IFPACK");
    MLList.set("smoother: type (level 3)", "IFPACK");
    MLList.set("smoother: type (level 4)", "IFPACK");
    MLList.set("smoother: ifpack type", "Amesos");
    MLList.set("smoother: ifpack overlap", 0);
  }

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  delete MLPrec;

  delete A;
  delete Map;

  return(solver.NumIters());
}

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
  MPI_Init(&argc,&argv);
  int MyPID, NumProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);

  if (NumProcs == 1) {
    cerr << "Please run this example with more than 1 proc" << endl;
    exit(EXIT_SUCCESS);
  }

  int ItersML = TestAdditiveSchwarz<Ifpack_Amesos>();
  int ItersAm = TestAdditiveSchwarz<ML_Epetra::Ifpack_ML>();

  if (MyPID == 0)
    cout << "ML iterations = " << ItersML
         << ", Amesos iterations = " << ItersAm << endl;

  int diff = ItersML - ItersAm;
  if (diff < 0) diff = -diff;
  if (diff > 5)
    exit(EXIT_FAILURE);
    
  ItersML = TestML(true);
  ItersAm = TestML(false);
  
  if (MyPID == 0)
    cout << "ML iterations = " << ItersML
         << ", Amesos iterations = " << ItersAm << endl;

  diff = ItersML - ItersAm;
  if (diff < 0) diff = -diff;
  if (diff > 5)
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
  puts("--enable-galeri");
  puts("--enable-ifpack");
  puts("--enable-mpi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) */
