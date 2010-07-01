
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
#include "Ifpack_IHSS.h"
#include "Ifpack_SORa.h"
#include "Ifpack.h"
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

  Teuchos::ParameterList GaleriList;
  int nx = 30; 

  GaleriList.set("nx", nx);
  //  GaleriList.set("ny", nx * Comm.NumProc());
  GaleriList.set("ny", nx);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());
  GaleriList.set("alpha", .0);
  GaleriList.set("diff", 1.0);
  GaleriList.set("conv", 100.0);

  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Cartesian2D", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("UniFlow2D", &*Map, GaleriList) );
  Teuchos::RefCountPtr<Epetra_MultiVector> LHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  Teuchos::RefCountPtr<Epetra_MultiVector> RHS = Teuchos::rcp( new Epetra_MultiVector(*Map, 1) );
  LHS->PutScalar(0.0); RHS->Random();
  Ifpack Factory;  
  int Niters = 100;

  // ============================= //
  // Construct IHSS preconditioner //
  // ============================= //
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create("IHSS", &*A,0) );
  Teuchos::ParameterList List;
  List.set("ihss: hermetian type","ILU");
  List.set("ihss: skew hermetian type","ILU");
  List.set("ihss: ratio eigenvalue",100.0);
  // Could set sublist values here to better control the ILU, but this isn't needed for this example.
  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Compute());

  // ============================= //
  // Create solver Object          //
  // ============================= //

  AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_gmres);
  solver.SetPrecOperator(&*Prec);
  solver.SetAztecOption(AZ_output, 1); 
  solver.Iterate(Niters, 1e-8);

  // ============================= //
  // Construct SORa preconditioner //
  // ============================= //
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec2 = Teuchos::rcp( Factory.Create("SORa", &*A,0) );
  Teuchos::ParameterList List2;
  List2.set("sora: sweeps",1);
  // Could set sublist values here to better control the ILU, but this isn't needed for this example.
  IFPACK_CHK_ERR(Prec2->SetParameters(List2));
  IFPACK_CHK_ERR(Prec2->Compute());

  // ============================= //
  // Create solver Object          //
  // ============================= //
  AztecOO solver2;
  LHS->PutScalar(0.0);
  solver2.SetUserMatrix(&*A);
  solver2.SetLHS(&*LHS);
  solver2.SetRHS(&*RHS);
  solver2.SetAztecOption(AZ_solver,AZ_gmres);
  solver2.SetPrecOperator(&*Prec2);
  solver2.SetAztecOption(AZ_output, 1); 
  solver2.Iterate(Niters, 1e-8);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
