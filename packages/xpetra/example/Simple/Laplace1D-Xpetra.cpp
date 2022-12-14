// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MultiVector.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>

#include <TpetraCore_config.h>

#if ((defined(HAVE_TPETRA_INST_OPENMP) || defined(HAVE_TPETRA_INST_SERIAL)) && \
    (defined(HAVE_TPETRA_INST_INT_INT) || defined(HAVE_TPETRA_INST_INT_LONG_LONG)) && \
    defined(HAVE_TPETRA_INST_DOUBLE))

// Choose types Tpetra is instantiated on
typedef double Scalar;
typedef int    LocalOrdinal;
#if defined(HAVE_TPETRA_INST_INT_INT)
typedef int    GlobalOrdinal;
#elif defined(HAVE_TPETRA_INST_INT_LONG_LONG)
typedef long long GlobalOrdinal;
#endif
#if defined(HAVE_TPETRA_INST_OPENMP)
typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;
#elif defined(HAVE_TPETRA_INST_SERIAL)
typedef Kokkos::Compat::KokkosSerialWrapperNode Node;
#endif

int main(int argc, char *argv[]) {

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = false;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    GlobalOrdinal numGlobalElements = 256; // problem size

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::createUniformContigMap(lib, numGlobalElements, comm);

    const size_t numMyElements = map->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

    RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A =  Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map, 3);

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

    A->fillComplete();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
#else
int main(int argc, char *argv[]) { std::cout << "Tpetra is not instantiated on SC=double, GO=int/long long and Node=Serial/OpenMP. Skip example." << std::endl; return EXIT_SUCCESS; }
#endif // Tpetra instantiated on SC=double, GO=int/long long and Node=Serial/OpenMP
