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
#include "Ifpack_ConfigDefs.h"

#if defined(HAVE_IFPACK_TEUCHOS) && defined(HAVE_IFPACK_TRIUTILS) && defined(HAVE_IFPACK_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include "Ifpack.h"

#ifdef HAVE_IFPACK_EPETRAEXT
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "Teuchos_CommandLineProcessor.hpp"
#endif
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
  
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);

  // define USE_EPETRAEXT if you want to read matrices, maps, solution and RHS
  // written by EpetraExt. This can be useful to test the basic performances
  // of IFPACK preconditioners for a given problem. See the EpetraExt
  // documentation for inout classes for more details
  // Otherwise, the matrix is created from the triutils gallery.

#if defined(USE_EPETRAEXT) && defined(HAVE_IFPACK_EPETRAEXT)
  Epetra_Time Time(Comm);
  Teuchos::CommandLineProcessor CLP;
  string MatrixFileName = "not-set";
  string MapFileName = "not-set";
  string LHSFileName = "not-set";
  string RHSFileName = "not-set";

  CLP.setOption("matrix", &MatrixFileName, "Matrix Market file containing matrix");
  CLP.setOption("map", &MapFileName, "Matrix Market file containing map");
  CLP.setOption("lhs", &LHSFileName, "Matrix Market file containing starting solution");
  CLP.setOption("rhs", &RHSFileName, "Matrix Market file containing rhs");

  CLP.throwExceptions(false);
  CLP.parse(argc,argv);

  if (MatrixFileName == "not-set" || MapFileName == "not-set" 
      || LHSFileName == "not-set" || RHSFileName == "not-set") {
    cerr << "You must specify the name of the file" << endl;
    cerr << "containing the matrix, the map, and rhs and the staring solution" << endl;
    cerr << "Run this example with option --help for more details" << endl;
    exit(EXIT_FAILURE);
  }

  Epetra_CrsMatrix* Matrix = 0;
  Epetra_Map* Map = 0;
  Epetra_MultiVector* LHS = 0;
  Epetra_MultiVector* RHS = 0;

  if (verbose)
    cout << "Reading map... " << endl;
  IFPACK_CHK_ERR(EpetraExt::MatrixMarketFileToMap(MapFileName.c_str(), Comm, Map));
  if (verbose)
    cout << "Reading matrix... " << endl;
  IFPACK_CHK_ERR(EpetraExt::MatrixMarketFileToCrsMatrix(MatrixFileName.c_str(), *Map, Matrix));
  if (verbose)
    cout << "Reading RHS... " << endl;
  IFPACK_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector(LHSFileName.c_str(),*Map,LHS));
  if (verbose)
    cout << "Reading LHS... " << endl;
  IFPACK_CHK_ERR(EpetraExt::MatrixMarketFileToMultiVector(RHSFileName.c_str(),*Map,RHS));

  if (verbose == 0)
    cout << "Time for reading data = " << Time.ElapsedTime() << endl;
#else
  // size of the global matrix (must be a square number)
  const int NumPoints = 10000;

  // build the matrix corresponding to a 2D Laplacian on a
  // structured grid.
  CrsMatrixGallery Gallery("recirc_2d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Epetra_RowMatrix*   Matrix = Gallery.GetMatrix();
  Epetra_MultiVector* LHS    = Gallery.GetStartingSolution();
  Epetra_MultiVector* RHS    = Gallery.GetRHS();
#endif

  // ==================== begin of IFPACK TEST part ===========================
  
  // We fix the following components of the preconditioner:
  // - type = ILU
  // - overlap = 0;
  // and we test the following:
  // - level-of-fill = 0, 1, 3, 5
  // - relative threshold = 0.0, 0.01 0.05 0.1
  // - absolute threshold = 0.0, 1e-9, 1e-6
  //
  // NOTE: THIS CAN BE LONG!!
  
  vector<int> lof;
  lof.push_back(0);
  lof.push_back(1);
  lof.push_back(3);
  lof.push_back(5);

  vector<double> rthresh;
  rthresh.push_back(1.0);
  rthresh.push_back(1.01);
  rthresh.push_back(1.05);
  rthresh.push_back(1.1);

  vector<double> athresh;
  athresh.push_back(0.0);
  athresh.push_back(1e-8);
  athresh.push_back(1e-4);
  athresh.push_back(1e-2);

  if (verbose) {
    cout << "Start testing..." << endl;
    cout << endl;
  }

  int count = 0;
  double BestTime = 1e+99;
  int BestTimeCount = -1;
  int BestIterations = 9999;
  int BestIterationsCount = -1;

  for (int i = 0 ; i < (int)lof.size() ; ++i) {
    for (int j = 0 ; j < (int)rthresh.size() ; ++j) {
      for (int k = 0 ; k < (int)athresh.size() ; ++k) {
        
        // parameters for this run
        int this_lof     = lof[i];
        double this_rthresh = rthresh[j];
        double this_athresh = athresh[k];
        if (verbose) {
          printf("[%3d] LOF = %d, rth = %5.2e, ath = %5.2e > ",
                 count, this_lof, this_rthresh, this_athresh);
        }
        // get a copy of the solution vector
        Epetra_MultiVector LHS2(*LHS);
        Epetra_LinearProblem Problem(Matrix,&LHS2,RHS);
        AztecOO solver(Problem);
        // define the IFPACK preconditioner
        Ifpack Factory;
        Ifpack_Preconditioner* IFPACKPrec;

        Epetra_Time Time(Comm);

        IFPACKPrec = Factory.Create("ILUT", Matrix);
        // set the parameters and build the preconditioner
        Teuchos::ParameterList IFPACKList;
        IFPACKList.set("fact: relative threshold", this_rthresh);
        IFPACKList.set("fact: absolute threshold", this_athresh);
        IFPACKList.set("fact: ict level-of-fill", (double)(this_lof+1));
        IFPACKList.set("schwarz: reordering type", "rcm");
        IFPACKList.set("schwarz: filter singletons", true);

        IFPACKPrec->SetParameters(IFPACKList);
        IFPACKPrec->Initialize();
        IFPACKPrec->Compute();
        // set for AztecOO
        solver.SetPrecOperator(IFPACKPrec);
        // other parameters for AztecOO
        solver.SetAztecOption(AZ_solver, AZ_gmres);
        solver.SetAztecOption(AZ_kspace, 500); 
        solver.SetAztecOption(AZ_output, AZ_none);
        // solve with 500 iterations and 1e-12 tolerance  
        solver.Iterate(500, 1e-5);

        int Iterations = solver.NumIters();
        double TotalTime = Time.ElapsedTime();

        if (verbose) {
          printf(" its = %5d      res = %9.2e     time = %f (s)\n",
                 Iterations, solver.TrueResidual(), TotalTime);
        }

        if (TotalTime < BestTime) {
          BestTime = TotalTime;
          BestTimeCount = count;
        }

        if (Iterations < BestIterations) {
          BestIterations = Iterations;
          BestIterationsCount = count;
        }

        ++count;
      }
    }
  }

  if (verbose) {
    cout << endl;
    cout << "Best iteration count was obtain in test # " << BestIterationsCount << endl; 
    cout << "Best CPU time was obtain in test # " << BestTimeCount << endl; 
    cout << endl;
  }

  // ===================== end of IFPACK TEST part =============================

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return(0);
  
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

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */

