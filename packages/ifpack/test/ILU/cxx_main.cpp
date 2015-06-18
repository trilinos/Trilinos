
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
#include "Ifpack_CrsRiluk.h"
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

  Teuchos::ParameterList GaleriList;
  int nx = 30; 

  GaleriList.set("nx", nx);
  GaleriList.set("ny", nx * Comm.NumProc());
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

  // I wanna test funky values to be sure that they have the same
  // influence on the algorithms, both old and new
  int    LevelFill = 2;
  double DropTol = 0.3333;
  double Athresh = 0.0123;
  double Rthresh = 0.9876;
  double Relax   = 0.1;
  int    Overlap = 2;
  
  Teuchos::RefCountPtr<Ifpack_IlukGraph> Graph;
  Teuchos::RefCountPtr<Ifpack_CrsRiluk> RILU;

  Graph = Teuchos::rcp( new Ifpack_IlukGraph(A->Graph(), LevelFill, Overlap) );
  int ierr;
  ierr = Graph->ConstructFilledGraph();
  IFPACK_CHK_ERR(ierr);

  RILU = Teuchos::rcp( new Ifpack_CrsRiluk(*Graph) );
  RILU->SetAbsoluteThreshold(Athresh);
  RILU->SetRelativeThreshold(Rthresh);
  RILU->SetRelaxValue(Relax);
  int initerr = RILU->InitValues(*A);
  if (initerr!=0) cout << Comm << "*ERR* InitValues = " << initerr;

  RILU->Factor();

  // Define label for printing out during the solve phase
  string label = "Ifpack_CrsRiluk Preconditioner: LevelFill = " + toString(LevelFill) +
                                                 " Overlap = " + toString(Overlap) +
                                                 " Athresh = " + toString(Athresh) +
                                                 " Rthresh = " + toString(Rthresh);
  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  int Niters = 1200;

  AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_gmres);
  solver.SetPrecOperator(&*RILU);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(Niters, 5.0e-5);

  int OldIters = solver.NumIters();

  // now rebuild the same preconditioner using RILU, we expect the same
  // number of iterations

  Ifpack Factory;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create("ILU", &*A, Overlap) );

  Teuchos::ParameterList List;
  List.get("fact: level-of-fill", LevelFill);
  List.get("fact: drop tolerance", DropTol);
  List.get("fact: absolute threshold", Athresh);
  List.get("fact: relative threshold", Rthresh);
  List.get("fact: relax value", Relax);

  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Compute());

  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_gmres);
  solver.SetPrecOperator(&*Prec);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(Niters, 5.0e-5);

  int NewIters = solver.NumIters();

  if (OldIters != NewIters)
    IFPACK_CHK_ERR(-1);


#ifdef HAVE_IFPACK_SUPERLU
  // Now test w/ SuperLU's ILU, if we've got it
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec2 = Teuchos::rcp( Factory.Create("SILU", &*A,0) ); 
  Teuchos::ParameterList SList;
  SList.set("fact: drop tolerance",1e-4);
  SList.set("fact: zero pivot threshold",.1);
  SList.set("fact: maximum fill factor",10.0);
  // NOTE: There is a bug in SuperLU 4.0 which will crash the code if the maximum fill factor is set too low.
  // This bug was reported to Sherry Li on 4/8/10.
  SList.set("fact: silu drop rule",9);

  IFPACK_CHK_ERR(Prec2->SetParameters(SList));
  IFPACK_CHK_ERR(Prec2->Compute());

  LHS->PutScalar(0.0);
  solver.SetPrecOperator(&*Prec2);
  solver.Iterate(Niters, 5.0e-5);
  Prec2->Print(cout);

#endif


#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
