
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
#include "Ifpack_CrsRiluk.h"
#include "Ifpack.h"
#include "Teuchos_RefCountPtr.hpp"

#include "Ifpack_HIPS.h"

extern "C"{
#include "hips.h"
};

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

  int ierr=HIPS_Initialize(1);
  HIPS_ExitOnError(ierr);

  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  Teuchos::ParameterList GaleriList;
  int nx = 100; 
  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
  //  GaleriList.set("ny", nx);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
                 
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Laplace2D", &*Map, GaleriList) );
  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  LHS->PutScalar(0.0); RHS->Random();

  // ============================ //
  // Construct ILU preconditioner //
  // ---------------------------- //

  Teuchos::RefCountPtr<Ifpack_HIPS> RILU;

  RILU = Teuchos::rcp( new Ifpack_HIPS(&*A) );

  Teuchos::ParameterList List;
  List.set("hips: id",0);
  List.set("hips: setup output",2);
  List.set("hips: iteration output",0);
  List.set("hips: drop tolerance",5e-3);
  List.set("hips: graph symmetric",1);
  RILU->SetParameters(List);

  
  RILU->Initialize();
  RILU->Compute();

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  int Niters = 50;

  AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_gmres);
  solver.SetPrecOperator(&*RILU);
  solver.SetAztecOption(AZ_output, 1); 
  solver.Iterate(Niters, 1.0e-8);

  int OldIters = solver.NumIters();


  HIPS_Finalize();
  
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
