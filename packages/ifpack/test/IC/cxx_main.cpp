
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
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

#include "Ifpack_ConfigDefs.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "AztecOO.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
//#include "Ifpack_CrsRick.h"
#include "Ifpack.h"
#include "Ifpack_DiagPreconditioner.h"
#include "Teuchos_RefCountPtr.hpp"

// function for fancy output

string toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

string toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

// main driver

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  // The problem is defined on a 2D grid, global size is nx * nx.
  int nx = 30;
  Teuchos::ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  LHS->PutScalar(0.0); RHS->Random();

  // ========================================= //
  // Compare IC preconditioners to no precond. //
  // ----------------------------------------- //

  const double tol = 1e-5;
  const int maxIter = 500;

  // Baseline: No preconditioning
  // Compute number of iterations, to compare to IC later.

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  //solver.SetPrecOperator(&*PrecDiag);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  int Iters = solver.NumIters();
  cout << "No preconditioner iterations: " << Iters << endl;

#if 0 
  // Not sure how to use Ifpack_CrsRick - leave out for now.
  //
  // I wanna test funky values to be sure that they have the same
  // influence on the algorithms, both old and new
  int    LevelFill = 2;
  double DropTol = 0.3333;
  double Condest;
  
  Teuchos::RefCountPtr<Ifpack_CrsRick> IC;
  Ifpack_IlukGraph mygraph (A->Graph(), 0, 0);
  IC = Teuchos::rcp( new Ifpack_CrsRick(*A, mygraph) );
  IC->SetAbsoluteThreshold(0.00123);
  IC->SetRelativeThreshold(0.9876);
  // Init values from A
  IC->InitValues(*A);
  // compute the factors
  IC->Factor();
  // and now estimate the condition number
  IC->Condest(false,Condest);
  
  if( Comm.MyPID() == 0 ) {
    cout << "Condition number estimate (level-of-fill = "
	 << LevelFill <<  ") = " << Condest << endl;
  }

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRick Preconditioner: LevelFill = " + toString(LevelFill) + 
                                                 " Overlap = 0"; 
  IC->SetLabel(label.c_str());
  
  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetPrecOperator(&*IC);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  int RickIters = solver.NumIters();
  cout << "Ifpack_Rick iterations: " << RickIters << endl;

  // Compare to no preconditioning
  if (RickIters > Iters/2)
    IFPACK_CHK_ERR(-1);

#endif

  //////////////////////////////////////////////////////
  // Same test with Ifpack_IC
  // This is Crout threshold Cholesky, so different than IC(0)

  Ifpack Factory;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> PrecIC = Teuchos::rcp( Factory.Create("IC", &*A) );

  Teuchos::ParameterList List;
  //List.get("fact: level-of-fill", 2);
  //List.get("fact: drop tolerance", 0.3333);
  //List.get("fact: absolute threshold", 0.00123);
  //List.get("fact: relative threshold", 0.9876);
  //List.get("fact: relaxation value", 0.0);

  IFPACK_CHK_ERR(PrecIC->SetParameters(List));
  IFPACK_CHK_ERR(PrecIC->Compute());

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  //AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetPrecOperator(&*PrecIC);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  int ICIters = solver.NumIters();
  cout << "Ifpack_IC iterations: " << ICIters << endl;

  // Compare to no preconditioning
  if (ICIters > Iters/2)
    IFPACK_CHK_ERR(-1);

  //////////////////////////////////////////////////////
  // Same test with Ifpack_ICT
  // This is threshold Cholesky, so different than IC(0)

  Teuchos::RefCountPtr<Ifpack_Preconditioner> PrecICT = Teuchos::rcp( Factory.Create("ICT", &*A) );

  //Teuchos::ParameterList List;
  //List.get("fact: level-of-fill", 2);
  //List.get("fact: drop tolerance", 0.3333);
  //List.get("fact: absolute threshold", 0.00123);
  //List.get("fact: relative threshold", 0.9876);
  //List.get("fact: relaxation value", 0.0);

  IFPACK_CHK_ERR(PrecICT->SetParameters(List));
  IFPACK_CHK_ERR(PrecICT->Compute());

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetPrecOperator(&*PrecICT);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  int ICTIters = solver.NumIters();
  cout << "Ifpack_ICT iterations: " << ICTIters << endl;

  // Compare to no preconditioning
  if (ICTIters > Iters/2)
    IFPACK_CHK_ERR(-1);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
