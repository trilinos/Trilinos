
// @HEADER
// ***********************************************************************
// 
//            Amesos: An Interface to Direct Solvers
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
// You must have configured Trilinos with --enable-galeri
// in order to compile this example

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Amesos.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Galeri_ReadHB.h"

// ==================== //
// M A I N  D R I V E R //
// ==================== //
//
// This example will:
// 1.- Read an H/B matrix from file;
// 2.- redistribute the linear system matrix to the
//     available processes
// 3.- set up LHS/RHS if not present
// 4.- create an Amesos_BaseSolver object
// 5.- solve the linear problem.
//
// This example can be run with any number of processors.
//
// Author: Marzio Sala, ETHZ/COLAB
// Last modified: Oct-05

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  std::string matrix_file;
  std::string solver_type;

  if (argc > 1)
    matrix_file = argv[1]; // read it from command line
  else
    matrix_file = "662_bus.rsa"; // file containing the HB matrix.

  if (argc > 2)
    solver_type = argv[2]; // read it form command line
  else
    solver_type = "Klu"; // default

  if (Comm.MyPID() == 0)
    std::cout << "Reading matrix `" << matrix_file << "'";
  
  // ================= //
  // reading HB matrix //
  // ================= //
  
  // HB files are for serial matrices. Hence, only
  // process 0 reads this matrix (and if present
  // solution and RHS). Then, this matrix will be redistributed
  // using epetra capabilities.
  // All variables that begin with "read" refer to the
  // HB matrix read by process 0.
  Epetra_Map* readMap;
  Epetra_CrsMatrix* readA; 
  Epetra_Vector* readx; 
  Epetra_Vector* readb;
  Epetra_Vector* readxexact;
  
  try 
  {
    Galeri::ReadHB(matrix_file.c_str(), Comm, readMap,
                   readA, readx, readb, readxexact);
  }
  catch(...)
  {
    std::cout << "Caught exception, maybe file name is incorrect" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#else
    // not to break test harness
    exit(EXIT_SUCCESS);
#endif
  }

  // Create uniform distributed map.
  // Note that linear map are used for simplicity only!
  // Amesos (through Epetra) can support *any* map.
  Epetra_Map* mapPtr = 0;
  
  if(readMap->GlobalIndicesInt())
    mapPtr = new Epetra_Map((int) readMap->NumGlobalElements(), 0, Comm);
  else if(readMap->GlobalIndicesLongLong())
    mapPtr = new Epetra_Map(readMap->NumGlobalElements(), 0, Comm);
  else
    assert(false);
  
  Epetra_Map& map = *mapPtr;
  
  // Create the distributed matrix, based on Map.
  Epetra_CrsMatrix A(Copy, map, 0);

  const Epetra_Map &OriginalMap = readA->RowMatrixRowMap() ; 
  assert (OriginalMap.SameAs(*readMap)); 
  Epetra_Export exporter(OriginalMap, map);

  Epetra_Vector x(map);          // distributed solution
  Epetra_Vector b(map);          // distributed rhs
  Epetra_Vector xexact(map);     // distributed exact solution

  // Exports from data defined on processor 0 to distributed.
  x.Export(*readx, exporter, Add);
  b.Export(*readb, exporter, Add);
  xexact.Export(*readxexact, exporter, Add);
  A.Export(*readA, exporter, Add);
  A.FillComplete();
    
  // deletes memory
  delete readMap;
  delete readA; 
  delete readx; 
  delete readb;
  delete readxexact;
  delete mapPtr;

  // Creates an epetra linear problem, contaning matrix
  // A, solution x and rhs b.
  Epetra_LinearProblem Problem(&A,&x,&b);
  
  // ======================================================= //
  // B E G I N N I N G   O F   T H E   A M E S O S   P A R T //
  // ======================================================= //

  if (!Comm.MyPID()) 
    std::cout << "Calling Amesos..." << std::endl;

  // For comments on the commands in this section, please
  // see file example_AmesosFactory.cpp.
  
  Amesos_BaseSolver* Solver = 0;
  Amesos Factory;
  
  Solver = Factory.Create(solver_type,Problem);

  // Factory.Create() returns 0 if the requested solver
  // is not available
  if (Solver == 0) {
    std::cerr << "Selected solver (" << solver_type << ") is not available" << std::endl;
    // return ok not to break test harness even if
    // the solver is not available
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return(EXIT_SUCCESS);
  }

  // Calling solve to compute the solution. This calls the symbolic
  // factorization and the numeric factorization.
  Solver->Solve();

  // Print out solver timings and get timings in parameter list.
  Solver->PrintStatus();
  Solver->PrintTiming();

  Teuchos::ParameterList TimingsList;
  Solver->GetTiming( TimingsList );

  // You can find out how much time was spent in ...
  double sfact_time, nfact_time, solve_time;
  double mtx_conv_time, mtx_redist_time, vec_redist_time;

  // 1) The symbolic factorization
  //    (parameter doesn't always exist)
  sfact_time = TimingsList.get( "Total symbolic factorization time", 0.0 );

  // 2) The numeric factorization
  //    (always exists if NumericFactorization() is called)
  nfact_time = Teuchos::getParameter<double>( TimingsList, "Total numeric factorization time" );

  // 3) Solving the linear system
  //    (always exists if Solve() is called)
  solve_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );

  // 4) Converting the matrix to the accepted format for the solver
  //    (always exists if SymbolicFactorization() is called)
  mtx_conv_time = Teuchos::getParameter<double>( TimingsList, "Total solve time" );

  // 5) Redistributing the matrix for each solve to the accepted format for the solver
  mtx_redist_time = TimingsList.get( "Total matrix redistribution time", 0.0 );

  // 6) Redistributing the vector for each solve to the accepted format for the solver
  vec_redist_time = TimingsList.get( "Total vector redistribution time", 0.0 );

  // =========================================== //
  // E N D   O F   T H E   A M E S O S   P A R T //
  // =========================================== //

  // Computes ||Ax - b|| //

  double residual;

  Epetra_Vector Ax(b.Map());
  A.Multiply(false, x, Ax);
  Ax.Update(1.0, b, -1.0);
  Ax.Norm2(&residual);

  if (!Comm.MyPID()) 
    std::cout << "After Amesos solution, ||b - A * x||_2 = " << residual << std::endl;

  // delete Solver. Do this before calling MPI_Finalize() because
  // MPI calls can occur.
  delete Solver;

  if (residual > 1e-5)
    return(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);

} // end of main()
