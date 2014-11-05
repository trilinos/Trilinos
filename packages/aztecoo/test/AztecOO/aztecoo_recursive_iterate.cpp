// Demonstrate growth in memory usage over repeated AztecOO solves
// It appears that calling "recursiveIterate" results in memory
//  associated with the Aztec internal preconditioner not being
//  cleaned up on construction.  This results in memory growth
//  over the course of multiple solves.

// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>

// Trilinos Includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif
#include "Epetra_LocalMap.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

using Teuchos::RCP;
using Teuchos::rcp;


// Build 3D Laplacian on Nx x Ny x Nz grid
RCP<Epetra_CrsMatrix> build_matrix(int Nx, int Ny, int Nz)
{
  int N = Nx*Ny*Nz;
  int entriesPerRow = 7;

  // No parallelism, so we have a serial communicator and a local map
#ifdef HAVE_MPI
  Epetra_MpiComm comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm comm;
#endif
  Epetra_LocalMap map( N, 0, comm );

  Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,entriesPerRow) );

  int err;
  std::vector<int>    ind(entriesPerRow);
  std::vector<double> val(entriesPerRow);
  for( int iz=0; iz<Nz; ++iz )
  {
    for( int iy=0; iy<Ny; ++iy )
    {
      for( int ix=0; ix<Nx; ++ix )
      {
        int row = ix + iy*Nx + iz*Nx*Ny;

        int col_ind = 0;

        ind[col_ind] = row;
        val[col_ind] = 6.0;
        col_ind++;

        if( ix>0 )
        {
          ind[col_ind] = row-1;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        if( ix < Nx-1 )
        {
          ind[col_ind] = row+1;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        if( iy>0 )
        {
          ind[col_ind] = row-Nx;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        if( iy < Ny-1 )
        {
          ind[col_ind] = row+Nx;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        if( iz>0 )
        {
          ind[col_ind] = row-Nx*Ny;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        if( iz < Nz-1 )
        {
          ind[col_ind] = row+Nx*Ny;
          val[col_ind]  = -1.0;
          col_ind++;
        }

        // Set row in matrix
        err = A->InsertGlobalValues( row, col_ind, &val[0], &ind[0] );
        if( err != 0 )
        {
          std::cout << "InsertGlobalValues returned error code "
                    << err << " for row " << row << std::endl;
        }
      }
    }
  }

  // Complete construction
  err = A->FillComplete();
  if( err != 0 )
    std::cout << "FillComplete returned error code " << err << std::endl;

  err = A->OptimizeStorage();
  if( err != 0 )
    std::cout << "OptimizeStorage returned error code " << err << std::endl;

  return A;
}


int main( int argc, char *argv[] )
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  using std::cout;

  using Teuchos::CommandLineProcessor;

  bool all_passed = true;

	try {

    CommandLineProcessor clp(false); // Don't throw exceptions

    const int numSolverTypes = 3;
    const int solverTypeValues[numSolverTypes] = { AZ_gmres, AZ_bicgstab, AZ_cg };
    const char* solverTypeNames[numSolverTypes] = { "gmres", "bicgstab", "cg" };
    int solverType = AZ_gmres;
    clp.setOption( "solver", &solverType,
      numSolverTypes, solverTypeValues, solverTypeNames,
      "The solver type to use." );

    int num_trials = 5;
    clp.setOption( "num-trials", &num_trials,
      "Number of times to set up and solve linear system" );

    int Ni= 5;
    clp.setOption( "Ni", &Ni, "Discretization size used for Nx, Ny, and Nz" );

    double tol = 1e-6;
    clp.setOption( "tol", &tol, "Tolerance for the linear solves" );

    int max_iters = 100;
    clp.setOption( "max-iters", &max_iters, "Maximum number of iterations" );

		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
		if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
			cout << "\nEnd Result: TEST FAILED" << std::endl;
			return parse_return;
		}

    for( int itrial=0; itrial<num_trials; ++itrial )
    {
      int Nx = Ni;
      int Ny = Ni;
      int Nz = Ni;
      int err;

      // Build matrix and vectors for this level
      RCP<Epetra_CrsMatrix> A = build_matrix(Nx,Ny,Nz);

      RCP<Epetra_MultiVector> x = rcp( new Epetra_MultiVector(A->OperatorDomainMap(),1) );
      RCP<Epetra_MultiVector> b = rcp( new Epetra_MultiVector(A->OperatorDomainMap(),1) );
      x->PutScalar(0.0);
      b->PutScalar(1.0);

      RCP<AztecOO> solver = rcp( new AztecOO() );

      int precType = AZ_ilu;
      if (solverType == AZ_gmres) {
        cout << "\nUsing AZ_gmres and AZ_ilu ...\n";
        precType = AZ_ilu;
      }
      else if (solverType == AZ_bicgstab) {
        cout << "\nUsing AZ_bicgstab and AZ_ilu ...\n";
        precType = AZ_ilu;
      }
      else if (solverType == AZ_cg) {
        cout << "\nUsing AZ_cg and AZ_icc ...\n";
        precType = AZ_icc;
      }

      solver->SetAztecOption(AZ_solver, solverType);
      solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
      solver->SetAztecOption(AZ_subdomain_solve, precType);

      // Prepare for solve
      solver->SetLHS( x.get() );
      solver->SetRHS( b.get() );
      err = solver->SetUserMatrix( A.get() );
      if( err != 0 ) {
        cout << "SetUserMatrix returned error code " << err << std::endl;
        all_passed = false;
      }

      // Explicitly build preconditioner
      double condest;
      cout << "\nConstructing a new preconditioner ...\n";
      err = solver->ConstructPreconditioner(condest);
      if( err != 0 ) {
        cout << "ConstructPreconditioner returned error code "
             << err << std::endl;
        all_passed = false;
      }

      solver->SetAztecOption(AZ_keep_info, 1);

      // Solve problem
      err = solver->recursiveIterate(max_iters,tol);
      if( err != 0 ) {
        cout << "recursiveIterate returned error code " << err << std::endl;
        all_passed = false;
      }

      // Destroy preconditioner
      err = solver->DestroyPreconditioner();
      if( err != 0 ) {
        cout << "DestroyPreconditioner returned error code "
             << err << std::endl;
        all_passed = false;
      }

      //cout << "End of solve, press Enter to continue..." << std::endl;
      //std::cin.ignore(100,'\n');

    }

	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, all_passed);


  if (all_passed) {
    cout << "ALL TESTS PASSED\n";
    return 0;
  }
  else {
    cout << "AT LEAST ONE TEST FAILED\n";
    return 1;
  }

}

