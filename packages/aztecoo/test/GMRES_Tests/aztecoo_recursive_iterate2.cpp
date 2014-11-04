
// Demonstrate memory leak in AztecOO using "recursiveIterate"
// The leak appears to be internal to the gmres solver, as switching
// to different solver options eliminates the errors.

// C++ Includes
#include <iostream>
#include <vector>
#include <cmath>

// Trilinos Includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

using Teuchos::RCP;
using Teuchos::rcp;

// Build 1D Laplacian of size Nx
Teuchos::RCP<Epetra_CrsMatrix> build_matrix(int Nx)
{
  int entriesPerRow = 3;

  Epetra_SerialComm comm;
  Epetra_LocalMap map( Nx, 0, comm );

  Teuchos::RCP<Epetra_CrsMatrix> A( new Epetra_CrsMatrix(Copy,map,entriesPerRow) );

  int err;
  std::vector<int>    ind(entriesPerRow);
  std::vector<double> val(entriesPerRow);
  for( int ix=0; ix<Nx; ++ix )
  {
    int col_ind = 0;

    ind[col_ind] = ix;
    val[col_ind] = 2.0;
    col_ind++;

    if( ix>0 )
    {
      ind[col_ind] = ix-1;
      val[col_ind]  = -1.0;
      col_ind++;
    }

    if( ix < Nx-1 )
    {
      ind[col_ind] = ix+1;
      val[col_ind]  = -1.0;
      col_ind++;
    }

    // Set row in matrix
    err = A->InsertGlobalValues( ix, col_ind, &val[0], &ind[0] );
    if( err != 0 )
    {
      std::cout << "InsertGlobalValues returned error code "
                << err << " for row " << ix << std::endl;
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
  using std::endl;

  bool all_passed = true;

	try {

    int Nx = 1000;
    int err;

    double tol = 1e-6;
    int max_iters = 100;

    // Build matrix and vectors for this level
    RCP<Epetra_CrsMatrix> A = build_matrix(Nx);

    RCP<Epetra_MultiVector> x = rcp( new Epetra_MultiVector(A->OperatorDomainMap(),1) );
    RCP<Epetra_MultiVector> b = rcp( new Epetra_MultiVector(A->OperatorDomainMap(),1) );
    x->PutScalar(0.0);
    b->PutScalar(1.0);

    RCP<AztecOO> solver = rcp( new AztecOO() );
    solver->SetAztecOption(AZ_solver,AZ_gmres);
    solver->SetAztecOption(AZ_precond,AZ_dom_decomp);
    solver->SetAztecOption(AZ_subdomain_solve,AZ_ilu);

    // Prepare for solve
    solver->SetLHS( x.get() );
    solver->SetRHS( b.get() );
    err = solver->SetUserMatrix( A.get() );
    if( err != 0 ) {
      cout << "SetUserMatrix returned error code " << err << endl;
      all_passed = false;
    }

    // Explicitly build preconditioner
    double condest;
    err = solver->ConstructPreconditioner(condest);
    if( err != 0 ) {
      cout << "ConstructPreconditioner returned error code "
                << err << endl;
      all_passed = false;
    }

    solver->SetAztecOption(AZ_keep_info, 1);

    // Solve problem
    err = solver->recursiveIterate(max_iters,tol);
    if( err != 0 ) {
      cout << "recursiveIterate returned error code " << err << endl;
      all_passed = false;
    }

    // Destroy preconditioner
    err = solver->DestroyPreconditioner();
    if( err != 0 ) {
      cout << "DestroyPreconditioner returned error code "
                << err << endl;
      all_passed = false;
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
