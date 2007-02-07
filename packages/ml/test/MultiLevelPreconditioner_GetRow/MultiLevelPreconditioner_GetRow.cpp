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

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_AMESOS)

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
#include "Epetra_Time.h"
#include "Galeri_CrsMatrices.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_VbrMatrices.h"
#include "Galeri_Utils.h"
#include "Amesos_TestRowMatrix.h"

using namespace Teuchos;
using namespace Galeri;

int solve(Epetra_RowMatrix&A, const bool UseIFPACK = true)
{
  Epetra_MultiVector X(A.OperatorDomainMap(), 2);
  Epetra_MultiVector B(A.OperatorRangeMap(), 2);

  X.PutScalar(0.0);
  B.PutScalar(1.0);

  Epetra_LinearProblem Problem(&A, &X, &B);

  AztecOO solver(Problem);

  ParameterList MLList;
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("aggregation: damping factor", 0.0);
  MLList.set("ML output", 0);

  if (UseIFPACK)
  {
    // Use IFPACK smoothers, because they can take advantage of the
    // underlying Epetra matrix. The "fast" SGS is used for Epetra_CrsMatrix.
    MLList.set("smoother: type","IFPACK");
    MLList.set("smoother: ifpack type", "point relaxation stand-alone");
    MLList.sublist("smoother: ifpack list").set("relaxation: type", "symmetric Gauss-Seidel");
  }

  Epetra_Time Time(A.Comm());
  Time.ResetStartTime();
  ML_Epetra::MultiLevelPreconditioner MLPrec(A, MLList);
  if (A.Comm().MyPID() == 0)
    cout << "time to build ML = " << Time.ElapsedTime() << endl;

  solver.SetPrecOperator(&MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-5);

  return(solver.NumIters());
}

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
  GaleriList.set("nx", 30);
  GaleriList.set("ny", 30 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1 * Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* CrsA = CreateCrsMatrix("Recirc2D", Map, GaleriList);
  Epetra_VbrMatrix* VbrA = CreateVbrMatrix(CrsA, 1);
  Amesos_TestRowMatrix* RowA = new Amesos_TestRowMatrix(CrsA);

  int RowIters = solve(*RowA);
  int CrsIters = solve(*CrsA);
  int VbrIters = solve(*VbrA);

  if (RowIters != CrsIters || RowIters != VbrIters || CrsIters != VbrIters)
  {
    if (Comm.MyPID() == 0)
      cout << "TEST FAILED!" << endl;
    exit(EXIT_FAILURE);
  }

  if (Comm.NumProc() == 1)
  {
    Epetra_VbrMatrix* VbrA2 = CreateVbrMatrix(CrsA, 5);
    Amesos_TestRowMatrix* RowA2 = new Amesos_TestRowMatrix(VbrA2);

    int VbrIters2 = solve(*VbrA2, false);
    int RowIters2 = solve(*RowA2, false);

    if (RowIters2 != VbrIters2)
    {
      if (Comm.MyPID() == 0)
        cout << "TEST FAILED!" << endl;
      exit(EXIT_FAILURE);
    }

    delete VbrA2;
    delete RowA2;
  }

  delete VbrA;
  delete CrsA;
  delete RowA;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  if (Comm.MyPID() == 0)
    cout << "TEST PASSED!" << endl;
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
  puts("--enable-amesos");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
}

#endif 
