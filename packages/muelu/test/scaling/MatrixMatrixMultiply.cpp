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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#ifndef _MSC_VER
#include <sys/time.h>
#endif

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"

namespace MueLuTests {

// generate random matrix
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > generateRandomMatrix(LocalOrdinal const &minEntriesPerRow, LocalOrdinal const &maxEntriesPerRow,
                                                                                              Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                                                                                              Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap);

// generate random map
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > generateRandomContiguousMap(size_t const minRowsPerProc, size_t const maxRowsPerProc,
                                                                                                Teuchos::RCP<const Teuchos::Comm<int> > const &comm, Xpetra::UnderlyingLib const lib);

// generate "random" whole number in interval [a,b]
template <class T>
size_t generateRandomNumber(T a, T b);

// generate distinct seeds for each process.  The seeds are different each run, but print to screen in case the run needs reproducibility.
unsigned int generateSeed(Teuchos::Comm<int> const &comm, const double initSeed = -1);

}  // namespace MueLuTests

//- -- --------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::StackedTimer;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using namespace MueLuTests;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

#include <MueLu_UseShortNames.hpp>
  Xpetra::Parameters xpetraParameters(clp);
  int optNmults = 1;
  clp.setOption("nmults", &optNmults, "number of matrix matrix multiplies to perform");

  GO optMaxRowsPerProc = 1000;
  clp.setOption("maxrows", &optMaxRowsPerProc, "maximum number of rows per processor");
  GO optMinRowsPerProc = 500;
  clp.setOption("minrows", &optMinRowsPerProc, "minimum number of rows per processor");
  LO optMaxEntriesPerRow = 50;
  clp.setOption("maxnnz", &optMaxEntriesPerRow, "maximum number of nonzeros per row");
  LO optMinEntriesPerRow = 25;
  clp.setOption("minnnz", &optMinEntriesPerRow, "minimum number of nonzeros per row");

  bool optDumpMatrices = false;
  clp.setOption("dump", "nodump", &optDumpMatrices, "write matrix to file");
  bool optTimings = false;
  clp.setOption("timings", "notimings", &optTimings, "print timings to screen");
  double optSeed = -1;
  clp.setOption("seed", &optSeed, "random number seed (cast to unsigned int)");

  std::ostringstream description;
  description << "This is a parallel-only test of the sparse matrix matrix multiply.  Two\n"
              << "\"random\" matrices are generated. They are most likely rectangular.  Each\n"
              << "has random contiguous domain and range maps. Each row has a random number of\n"
              << "entries with random global column IDs. The matrix values are identically 1,\n"
              << "and success is declared if the multiply finishes.";
  clp.setDocString(description.str().c_str());

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  unsigned int seed = generateSeed(*comm, optSeed);
  Teuchos::ScalarTraits<SC>::seedrandom(seed);

  RCP<StackedTimer> timer = rcp(new StackedTimer("MatrixMatrix Multiply: Total"));
  TimeMonitor::setStackedTimer(timer);

  for (int jj = 0; jj < optNmults; ++jj) {
    RCP<Matrix> A;
    RCP<Matrix> B;
    {
      TimeMonitor tm(*TimeMonitor::getNewTimer("MatrixMatrixMultiplyTest: 1 - Matrix creation"));

      size_t maxRowsPerProc = optMaxRowsPerProc;
      size_t minRowsPerProc = optMinRowsPerProc;
      if (minRowsPerProc > maxRowsPerProc) minRowsPerProc = maxRowsPerProc;

      // Create row map.  This will also be used as the range map.
      RCP<const Map> rowMapForA = generateRandomContiguousMap<LO, GO, NO>(minRowsPerProc, maxRowsPerProc, comm, xpetraParameters.GetLib());
      // Create domain map for A. This will also be the row map for B.
      RCP<const Map> domainMapForA = generateRandomContiguousMap<LO, GO, NO>(minRowsPerProc, maxRowsPerProc, comm, xpetraParameters.GetLib());
      // Create domain map for B.
      RCP<const Map> domainMapForB = generateRandomContiguousMap<LO, GO, NO>(minRowsPerProc, maxRowsPerProc, comm, xpetraParameters.GetLib());

      A = generateRandomMatrix<SC, LO, GO, NO>(optMinEntriesPerRow, optMaxEntriesPerRow, rowMapForA, domainMapForA);
      B = generateRandomMatrix<SC, LO, GO, NO>(optMinEntriesPerRow, optMaxEntriesPerRow, domainMapForA, domainMapForB);

      if (comm->getRank() == 0) {
        std::cout << "case " << jj << " of " << optNmults - 1 << " : "
                  << "seed = " << seed
                  << ", minEntriesPerRow = " << optMinEntriesPerRow
                  << ", maxEntriesPerRow = " << optMaxEntriesPerRow
                  << std::endl
                  << "  A stats:"
                  << " #rows = " << rowMapForA->getGlobalNumElements()
                  << " #cols = " << domainMapForA->getGlobalNumElements()
                  << ", nnz = " << A->getGlobalNumEntries()
                  << std::endl
                  << "  B stats:"
                  << " #rows = " << domainMapForA->getGlobalNumElements()
                  << " #cols = " << domainMapForB->getGlobalNumElements()
                  << ", nnz = " << B->getGlobalNumEntries()
                  << std::endl;
      }

      if (optDumpMatrices) {
        std::string fileName = "checkA.mm";
        Xpetra::IO<SC, LO, GO, Node>::Write(fileName, *A);
        fileName = "checkB.mm";
        Xpetra::IO<SC, LO, GO, Node>::Write(fileName, *B);
      }

    }  // scope for timing matrix creation

    {
      TimeMonitor tm(*TimeMonitor::getNewTimer("MatrixMatrixMultiplyTest: 2 - Multiply"));

      RCP<Matrix> AB;
      RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(rcp(new Teuchos::oblackholestream()));
      AB                             = Xpetra::MatrixMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Multiply(*A, false, *B, false, *fos);
      // if (optDumpMatrices) {
      //   std::string fileName="checkAB.mm";
      //   Utils::Write( fileName,*AB);
      // }

    }  // scope for multiply

    seed = generateSeed(*comm);

  }  // for (int jj=0; jj<optNmults; ++jj)

  timer->stopBaseTimer();
  StackedTimer::OutputOptions options;
  options.print_warnings = false;
  timer->report(std::cout, comm, options);

  auto xmlOut = timer->reportWatchrXML(std::string("MueLu MatrixMatrix Multiply ") + std::to_string(comm->getSize()) + " ranks", comm);
  if (xmlOut.length())
    std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';

  if (comm->getRank() == 0)
    std::cout << "End Result: TEST PASSED";

  return 0;
}
// main_

namespace MueLuTests {

template <class T>
size_t generateRandomNumber(T a, T b) {
  double rv            = (Teuchos::ScalarTraits<double>::random() + 1) * 0.5;  // shift "random" value from interval [-1,1] to [0,1]
  size_t numMyElements = Teuchos::as<size_t>(a + rv * (b - a));
  return numMyElements;
}

//- -- --------------------------------------------------------

unsigned int generateSeed(const Teuchos::Comm<int> &comm, const double initSeed) {
  timeval t1;
  gettimeofday(&t1, NULL);

  unsigned int seed;
  if (initSeed > -1)
    seed = Teuchos::as<unsigned int>(initSeed);
  else
    seed = t1.tv_usec * t1.tv_sec;

  const Teuchos::MpiComm<int> *mpiComm = dynamic_cast<const Teuchos::MpiComm<int> *>(&comm);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiComm == 0, MueLu::Exceptions::RuntimeError, "Cast to MpiComm failed");

  // use variant of proc 0's seed so we can always reproduce the results
  MPI_Bcast((void *)&seed, 1, MPI_UNSIGNED, 0, *(mpiComm->getRawMpiComm()));
  seed = seed * (1 + comm.getRank());

  return seed;
}

//- -- --------------------------------------------------------
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > generateRandomContiguousMap(size_t const minRowsPerProc, size_t const maxRowsPerProc, Teuchos::RCP<const Teuchos::Comm<int> > const &comm, Xpetra::UnderlyingLib const lib) {
  size_t numMyElements                                                    = generateRandomNumber(minRowsPerProc, maxRowsPerProc);
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > map = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::createContigMap(lib,
                                                                                                                                                   Teuchos::OrdinalTraits<size_t>::invalid(),
                                                                                                                                                   numMyElements,
                                                                                                                                                   comm);
  return map;
}

//- -- --------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > generateRandomMatrix(LocalOrdinal const &minEntriesPerRow, LocalOrdinal const &maxEntriesPerRow,
                                                                                              Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap,
                                                                                              Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &domainMap) {
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
#include <MueLu_UseShortNames.hpp>

  size_t numLocalRowsInA = rowMap->getLocalNumElements();
  // Now allocate random number of entries per row for each processor, between minEntriesPerRow and maxEntriesPerRow
  Teuchos::ArrayRCP<size_t> eprData(numLocalRowsInA);
  for (Teuchos::ArrayRCP<size_t>::iterator i = eprData.begin(); i != eprData.end(); ++i) {
    *i = generateRandomNumber(minEntriesPerRow, maxEntriesPerRow);
  }

  // Populate CrsMatrix with random nnz per row and with random column locations.
  ArrayView<const GlobalOrdinal> myGlobalElements = rowMap->getLocalElementList();
  GO numGlobalRows                                = rowMap->getGlobalNumElements();

  LO realMaxEntriesPerRow = maxEntriesPerRow;
  if (maxEntriesPerRow > numGlobalRows)
    realMaxEntriesPerRow = numGlobalRows;
  RCP<Matrix> A = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, eprData));

  Array<Scalar> vals(realMaxEntriesPerRow);
  // stick in ones for values
  for (LO j = 0; j < realMaxEntriesPerRow; ++j) vals[j] = Teuchos::ScalarTraits<SC>::one();
  Array<GO> cols(realMaxEntriesPerRow);
  for (size_t i = 0; i < numLocalRowsInA; ++i) {
    ArrayView<SC> av(&vals[0], eprData[i]);
    ArrayView<GO> iv(&cols[0], eprData[i]);
    for (size_t j = 0; j < eprData[i]; ++j) {
      // we assume there are no gaps in the indices
      cols[j] = Teuchos::as<GO>(generateRandomNumber(Teuchos::as<GO>(0), domainMap->getMaxAllGlobalIndex()));
    }
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete(domainMap, rowMap);

  return A;
}

}  // namespace MueLuTests

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
