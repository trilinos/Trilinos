
/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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

std::string toString(const int& x) {
  char s[100];
  sprintf(s, "%d", x);
  return std::string(s);
}

std::string toString(const double& x) {
  char s[100];
  sprintf(s, "%g", x);
  return std::string(s);
}

// main driver

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

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
  List2.set("sora: use global damping",true);
  List2.set("sora: eigen-analysis: random seed",(unsigned int)24601);
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

  // ============================= //
  // Construct a second SORa preconditioner to check seeds //
  // ============================= //
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec3 = Teuchos::rcp( Factory.Create("SORa", &*A,0) );
  Teuchos::ParameterList List3;
  List3.set("sora: sweeps",1);
  List3.set("sora: use global damping",true);
  List3.set("sora: eigen-analysis: random seed",(unsigned int)24601);
  // Could set sublist values here to better control the ILU, but this isn't needed for this example.
  IFPACK_CHK_ERR(Prec3->SetParameters(List2));
  IFPACK_CHK_ERR(Prec3->Compute());

  Teuchos::RCP<Ifpack_SORa> Prec2_SORa = Teuchos::rcp_dynamic_cast<Ifpack_SORa>(Prec2);
  Teuchos::RCP<Ifpack_SORa> Prec3_SORa = Teuchos::rcp_dynamic_cast<Ifpack_SORa>(Prec3);
  double diff = Prec2_SORa->GetLambdaMax()-Prec3_SORa->GetLambdaMax();
  if(diff > 1e-12) return EXIT_FAILURE;

 
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
