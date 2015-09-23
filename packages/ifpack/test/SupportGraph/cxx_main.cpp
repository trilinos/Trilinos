
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
//#include "Ifpack_CrsRick.h"
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
  
  int MyPID = Comm.MyPID();
  bool verbose = false; 
  if (MyPID==0) verbose = true;

  
  /*int npRows = -1;
  int npCols = -1;
  bool useTwoD = false;
  int randomize = 1;
  std::string matrix = "Laplacian";
  
  Epetra_CrsMatrix *AK = NULL;
    std::string filename = "email.mtx";
  read_matrixmarket_file((char*) filename.c_str(), Comm, AK,
			 useTwoD, npRows, npCols,
			 randomize, false,
			 (matrix.find("Laplacian")!=std::string::npos));
  Teuchos::RCP<Epetra_CrsMatrix> A(AK);
  const Epetra_Map *AMap = &(AK->DomainMap());
  Teuchos::RCP<const Epetra_Map> Map(AMap, false);*/

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

  // ==================================================== //
  // Compare support graph preconditioners to no precond. //
  // ---------------------------------------------------- //

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
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  int Iters = solver.NumIters();


  int SupportIters;
  Ifpack Factory;
  Teuchos::ParameterList List;

#ifdef HAVE_IFPACK_AMESOS
  //////////////////////////////////////////////////////
  // Same test with Ifpack_SupportGraph
  // Factored with Amesos

  
  Teuchos::RefCountPtr<Ifpack_Preconditioner> PrecSupportAmesos = Teuchos::rcp( Factory.Create("MSF Amesos", &*A) );
  List.set("amesos: solver type","Klu");
  List.set("MST: keep diagonal", 1.0);
  List.set("MST: randomize", 1);
  //List.set("fact: absolute threshold", 3.0);
  
  IFPACK_CHK_ERR(PrecSupportAmesos->SetParameters(List));
  IFPACK_CHK_ERR(PrecSupportAmesos->Initialize());
  IFPACK_CHK_ERR(PrecSupportAmesos->Compute());


  // Here we create an AztecOO object
  LHS->PutScalar(0.0);

  //AztecOO solver;
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetPrecOperator(&*PrecSupportAmesos);
  solver.SetAztecOption(AZ_output, 16); 
  solver.Iterate(maxIter, tol);

  SupportIters = solver.NumIters();




  
  // Compare to no preconditioning
  if (SupportIters > 2*Iters)
    IFPACK_CHK_ERR(-1);

#endif

  //////////////////////////////////////////////////////
  // Same test with Ifpack_SupportGraph
  // Factored with IC
  

  
  Teuchos::RefCountPtr<Ifpack_Preconditioner> PrecSupportIC = Teuchos::rcp( Factory.Create("MSF IC", &*A) );

  

  IFPACK_CHK_ERR(PrecSupportIC->SetParameters(List));
  IFPACK_CHK_ERR(PrecSupportIC->Compute());


  // Here we create an AztecOO object                                                                                                                                                                                                        
  LHS->PutScalar(0.0);

  //AztecOO solver;                                                                                                                                                                                                                          
  solver.SetUserMatrix(&*A);
  solver.SetLHS(&*LHS);
  solver.SetRHS(&*RHS);
  solver.SetAztecOption(AZ_solver,AZ_cg);
  solver.SetPrecOperator(&*PrecSupportIC);
  solver.SetAztecOption(AZ_output, 16);
  solver.Iterate(maxIter, tol);

  SupportIters = solver.NumIters();

  // Compare to no preconditioning                                                                                                                                                                                                           
  if (SupportIters > 2*Iters)
    IFPACK_CHK_ERR(-1);
  




#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
