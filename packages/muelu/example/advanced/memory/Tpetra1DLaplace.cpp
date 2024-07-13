// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "MueLu_MemoryProfiler.hpp"

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard mpiSession(&argc, &argv);

  bool success = false;
  bool verbose = true;
  try {
    typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
    typedef Tpetra::Map<>::local_ordinal_type LO;
#if defined(HAVE_TPETRA_INST_INT_INT)
    // mfh 07 Aug 2015: Prefer GO = int, for consistency with Epetra,
    // but use the default GO type if GO = int is not enabled.
    typedef int GO;
#else
    typedef Tpetra::Map<LO>::global_ordinal_type GO;
#endif  // HAVE_TPETRA_INT_INT
    typedef Tpetra::Map<LO, GO>::node_type Node;
    typedef Tpetra::Map<LO, GO, Node> Map;
    typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> CrsMatrix;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::tuple;

    RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // const int myRank = comm->getRank();

    // int numGlobalElements = 10000000;
    int numGlobalElements = 100;

    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("numGlobalElements", &numGlobalElements, "Global problem size.");
    if (cmdp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    const GO indexBase = 0;
    RCP<const Map> map =
        rcp(new Map(static_cast<Tpetra::global_size_t>(numGlobalElements),
                    indexBase, comm));
    const size_t numMyElements                    = map->getLocalNumElements();
    Teuchos::ArrayView<const GO> myGlobalElements = map->getLocalElementList();

    MemoryUsageStart("Tpetra");
    PrintMemoryUsage("Initial memory usage", "tpetra-init.heap");

    RCP<CrsMatrix> A = Tpetra::createCrsMatrix<Scalar>(map, 3);

    PrintMemoryUsage("Memory after CrsMatrix constructor", "tpetra-after-ctor.heap");

    for (size_t i = 0; i < numMyElements; i++) {
      if (myGlobalElements[i] == 0) {
        A->insertGlobalValues(myGlobalElements[i],
                              tuple<GO>(myGlobalElements[i], myGlobalElements[i] + 1),
                              tuple<Scalar>(2.0, -1.0));
      } else if (myGlobalElements[i] == numGlobalElements - 1) {
        A->insertGlobalValues(myGlobalElements[i],
                              tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i]),
                              tuple<Scalar>(-1.0, 2.0));
      } else {
        A->insertGlobalValues(myGlobalElements[i],
                              tuple<GO>(myGlobalElements[i] - 1, myGlobalElements[i], myGlobalElements[i] + 1),
                              tuple<Scalar>(-1.0, 2.0, -1.0));
      }
    }

    PrintMemoryUsage("Memory after insertGlobalValues()", "tpetra-after-insert.heap");

    A->fillComplete();  // DoOptimizeStorage by default

    PrintMemoryUsage("Memory after fillComplete()", "tpetra-after-fillcomplete.heap");

    MemoryUsageStop();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
