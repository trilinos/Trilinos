// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#ifdef HAVE_XPETRA_TPETRA

#include <TpetraCore_config.h>

#if ((defined(HAVE_TPETRA_INST_OPENMP) || defined(HAVE_TPETRA_INST_SERIAL)) &&         \
     (defined(HAVE_TPETRA_INST_INT_INT) || defined(HAVE_TPETRA_INST_INT_LONG_LONG)) && \
     defined(HAVE_TPETRA_INST_DOUBLE))

// Choose types Tpetra is instantiated on
typedef double Scalar;
typedef int LocalOrdinal;
#if defined(HAVE_TPETRA_INST_INT_INT)
typedef int GlobalOrdinal;
#elif defined(HAVE_TPETRA_INST_INT_LONG_LONG)
typedef long long GlobalOrdinal;
#endif
#if defined(HAVE_TPETRA_INST_OPENMP)
typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode Node;
#elif defined(HAVE_TPETRA_INST_SERIAL)
typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode Node;
#endif

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = false;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    GlobalOrdinal numGlobalElements = 256;  // problem size

    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;

    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::createUniformContigMap(lib, numGlobalElements, comm);

    const size_t numMyElements                               = map->getLocalNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> myGlobalElements = map->getLocalElementList();

    RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A = Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map, 3);

    for (size_t i = 0; i < numMyElements; i++) {
      if (myGlobalElements[i] == 0) {
        A->insertGlobalValues(myGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i], myGlobalElements[i] + 1),
                              Teuchos::tuple<Scalar>(2.0, -1.0));
      } else if (myGlobalElements[i] == numGlobalElements - 1) {
        A->insertGlobalValues(myGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] - 1, myGlobalElements[i]),
                              Teuchos::tuple<Scalar>(-1.0, 2.0));
      } else {
        A->insertGlobalValues(myGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                              Teuchos::tuple<Scalar>(-1.0, 2.0, -1.0));
      }
    }

    A->fillComplete();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
#else
int main(int argc, char *argv[]) {
  std::cout << "Tpetra is not instantiated on SC=double, GO=int/long long and Node=Serial/OpenMP. Skip example." << std::endl;
  return EXIT_SUCCESS;
}
#endif  // Tpetra instantiated on SC=double, GO=int/long long and Node=Serial/OpenMP
#else
int main(int argc, char *argv[]) {
  std::cout << "Xpetra has been compiled without Tpetra support. Skip example." << std::endl;
  return EXIT_SUCCESS;
}
#endif  // HAVE_XPETRA_TPETRA
