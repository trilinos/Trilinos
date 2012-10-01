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
#include <iostream>
//#include <complex>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>

// Define default template types
//typedef std::complex<double> Scalar;
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  // problem size
  GlobalOrdinal numGlobalElements = 256;

  // linear algebra library
  // this example require Tpetra+Ifpack2+Amesos2 or Epetra+Ifpack+Amesos
#if   defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2) && defined(HAVE_MUELU_AMESOS2)
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra; 
#elif defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)  && defined(HAVE_MUELU_AMESOS)
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
#endif

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal>::createUniformContigMap(lib, numGlobalElements, comm);

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal> > A = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal>(map, 3));

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<Scalar> (2.0, -1.0));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i]), 
                            Teuchos::tuple<Scalar> (-1.0, 2.0));
    }
    else {
      A->insertGlobalValues(myGlobalElements[i], 
                            Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1), 
                            Teuchos::tuple<Scalar> (-1.0, 2.0, -1.0));
    }
  }

  // Complete the fill, ask that storage be reallocated and optimized
  A->fillComplete();

  //
  // Construct a multigrid preconditioner
  //

  // Multigrid Hierarchy
  MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal> H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  H.Setup();

  //
  // Solve Ax = b
  //

  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > X = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal>::Build(map);
  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal> > B = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal>::Build(map);
  
  X->putScalar((Scalar) 0.0);
  B->setSeed(846930886); B->randomize();

  // Use AMG directly as an iterative solver (not as a preconditionner)
  int nIts = 9;

  H.Iterate(*B, nIts, *X);

  // Print relative residual norm
  typedef Teuchos::ScalarTraits<Scalar> ST;
  ST::magnitudeType residualNorms = MueLu::Utils<Scalar, LocalOrdinal, GlobalOrdinal>::ResidualNorm(*A, *X, *B)[0];
  if (comm->getRank() == 0)
    std::cout << "||Residual|| = " << residualNorms << std::endl;

  return EXIT_SUCCESS;
}
