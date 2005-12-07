// @HEADER
// ***********************************************************************
// 
//                Amesos: An Interface to Direct Solvers
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Amesos_ConfigDefs.h"
// This example needs Galeri to generate the linear system.
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"

using namespace Teuchos;
using namespace Galeri;

// ===================== //
// M A I N   D R I V E R //
// ===================== //
//
// This example compares all the available Amesos solvers
// for the solution of the same linear system. 
//
// The example can be run in serial and in parallel.
//
// Author: Marzio Sala, SNL 9214
// Last modified: Oct-05.

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool verbose = (Comm.MyPID() == 0);
  double TotalResidual = 0.0;

  // Create the Map, defined as a grid, of size nx x ny x nz,
  // subdivided into mx x my x mz cubes, each assigned to a 
  // different processor.

  ParameterList GaleriList;
  GaleriList.set("nx", 4);
  GaleriList.set("ny", 4);
  GaleriList.set("nz", 4 * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", 1);
  GaleriList.set("mz", Comm.NumProc());
  Epetra_Map* Map = CreateMap("Cartesian3D", Comm, GaleriList);
  
  // Create a matrix, in this case corresponding to a 3D Laplacian
  // discretized using a classical 7-point stencil. Please refer to 
  // the Galeri documentation for an overview of available matrices.
  // 
  // NOTE: matrix must be symmetric if DSCPACK is used.
  Epetra_CrsMatrix* Matrix = CreateCrsMatrix("Laplace3D", Map, GaleriList);

  // build vectors, in this case with 1 vector
  Epetra_MultiVector LHS(*Map, 1); 
  Epetra_MultiVector RHS(*Map, 1); 
    
  // create a linear problem object
  Epetra_LinearProblem Problem(Matrix, &LHS, &RHS);

  // use this list to set up parameters, now it is required
  // to use all the available processes (if supported by the
  // underlying solver). Uncomment the following two lines
  // to let Amesos print out some timing and status information.
  ParameterList List;
  List.set("PrintTiming",true);
  List.set("PrintStatus",true);
  List.set("MaxProcs",Comm.NumProc());

  vector<string> SolverType;
  SolverType.push_back("Amesos_Klu");
  SolverType.push_back("Amesos_Lapack");
  SolverType.push_back("Amesos_Umfpack");
  SolverType.push_back("Amesos_Pardiso");
  SolverType.push_back("Amesos_Taucs");
  SolverType.push_back("Amesos_Superlu");
  SolverType.push_back("Amesos_Superludist");
  SolverType.push_back("Amesos_Mumps");
  SolverType.push_back("Amesos_Dscpack");
  SolverType.push_back("Amesos_Scalapack");

  Epetra_Time Time(Comm);
  
  // this is the Amesos factory object that will create 
  // a specific Amesos solver.
  Amesos Factory;

  // Cycle over all solvers.
  // Only installed solvers will be tested.
  for (unsigned int i = 0 ; i < SolverType.size() ; ++i) 
  {
    // Check whether the solver is available or not
    if (Factory.Query(SolverType[i])) 
    {
      // 1.- set exact solution (constant vector)
      LHS.PutScalar(1.0);
 
      // 2.- create corresponding rhs
      Matrix->Multiply(false, LHS, RHS);
 
      // 3.- randomize solution vector
      LHS.Random();
 
      // 4.- create the amesos solver object
      Amesos_BaseSolver* Solver = Factory.Create(SolverType[i], Problem);
      assert (Solver != 0);

      Solver->SetParameters(List);

      // 5.- factorize and solve
      
      Time.ResetStartTime();
      AMESOS_CHK_ERR(Solver->SymbolicFactorization());
      if (verbose) 
        cout << endl
             << "Solver " << SolverType[i] 
             << ", symbolic factorization time = " 
             << Time.ElapsedTime() << endl;

      AMESOS_CHK_ERR(Solver->NumericFactorization());
      if (verbose) 
        cout << "Solver " << SolverType[i] 
             << ", numeric factorization time = " 
             << Time.ElapsedTime() << endl;

      AMESOS_CHK_ERR(Solver->Solve());
      if (verbose) 
        cout << "Solver " << SolverType[i] 
             << ", solve time = " 
             << Time.ElapsedTime() << endl;
  
      // 6.- compute difference between exact solution and Amesos one
      //     (there are other ways of doing this in Epetra, but let's
      //     keep it simple)
      double d = 0.0, d_tot = 0.0;
      for (int j = 0 ; j< LHS.Map().NumMyElements() ; ++j)
        d += (LHS[0][j] - 1.0) * (LHS[0][j] - 1.0);

      Comm.SumAll(&d,&d_tot,1);
      if (verbose)
        cout << "Solver " << SolverType[i] << ", ||x - x_exact||_2 = " 
             << sqrt(d_tot) << endl;

      // 7.- delete the object
      delete Solver;

      TotalResidual += d_tot;
    }
  }

  delete Matrix;
  delete Map;

  if (TotalResidual > 1e-9) 
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
} // end of main()
