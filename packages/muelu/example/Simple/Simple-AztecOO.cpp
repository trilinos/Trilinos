// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// To compile and run this example, the optional dependency AztecOO must be available.
//

#include <iostream>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

#include <MueLu_TrilinosSmoother.hpp> //TODO: remove

// Header files defining default types for template parameters.
// These headers must be included after other MueLu/Xpetra headers.
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=int, GlobalOrdinal=int
#include <MueLu_UseShortNames.hpp>    // => typedef MueLu::FooClass<Scalar, LocalOrdinal, ...> Foo

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <MueLu_EpetraOperator.hpp>

int main(int argc, char *argv[]) {

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp; //

  //
  // MPI initialization using Teuchos
  //

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  //
  // Parameters
  //

  GlobalOrdinal numGlobalElements = 256; // problem size

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor
  const Epetra_Map map(numGlobalElements, 0, comm);

  // Get update list and number of local equations from newly created map.
  const size_t         numMyElements    = map.NumMyElements();
  const GlobalOrdinal* myGlobalElements = map.MyGlobalElements();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<Epetra_CrsMatrix> A = rcp(new Epetra_CrsMatrix(Copy, map, 3));

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {

      //TODO: should be rewritten in an Epetra style
      A->InsertGlobalValues(myGlobalElements[i], 2, 
                            Teuchos::tuple<Scalar> (2.0, -1.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1).getRawPtr());

    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->InsertGlobalValues(myGlobalElements[i], 2,
                            Teuchos::tuple<Scalar> (-1.0, 2.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]).getRawPtr());
    }
    else {
      A->InsertGlobalValues(myGlobalElements[i], 3,
                            Teuchos::tuple<Scalar> (-1.0, 2.0, -1.0).getRawPtr(),
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1).getRawPtr());
    }
  }

  // Complete the fill, ask that storage be reallocated and optimized
  A->FillComplete();

  //
  // Construct a multigrid preconditioner
  //

  // Turns a Epetra_CrsMatrix into a MueLu::Matrix
  RCP<Xpetra::CrsMatrix<SC, LO, GO, NO, LMO> > mueluA_ = rcp(new Xpetra::EpetraCrsMatrix(A)); //TODO: should not be needed
  RCP<Xpetra::Matrix <SC, LO, GO, NO, LMO> > mueluA  = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO, LMO>(mueluA_));

  // Multigrid Hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy(mueluA));
  H->setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  H->Setup();

  //
  // Define RHS / LHS
  //

  RCP<Epetra_Vector> X = rcp(new Epetra_Vector(map));
  RCP<Epetra_Vector> B = rcp(new Epetra_Vector(map));
  
  X->PutScalar((Scalar) 0.0);
  B->SetSeed(846930886); B->Random();
  //B->PutScalar((Scalar) 1.0);

  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //

  Epetra_LinearProblem epetraProblem(A.get(), X.get(), B.get());

  AztecOO aztecSolver(epetraProblem);
  aztecSolver.SetAztecOption(AZ_solver, AZ_cg);

  MueLu::EpetraOperator aztecPrec(H);
  aztecSolver.SetPrecOperator(&aztecPrec);

  /* example of options
    aztecSolver.SetAztecParam(AZ_rthresh, 1.4);
    aztecSolver.SetAztecParam(AZ_athresh, 10.0);
    aztecSolver.SetAztecParam(AZ_ill_cond_thresh, 1.0e200);
    aztecSolver.SetAztecOption(AZ_kspace, 160);
  */

  int maxIts = 20;
  double tol = 1e-4;

  aztecSolver.Iterate(maxIts, tol);

  //TODO: test solution

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
