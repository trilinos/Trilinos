// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_BinaryIO.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Tpetra::global_size_t;

std::string makeFilename(const std::string& stem,
                         const Teuchos::Comm<int>& comm,
                         const char suffix[]) {
  std::ostringstream os;
  os << stem << "_p" << comm.getSize() << suffix;
  return os.str();
}

void cleanupFile(const std::string& filename,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                 const bool keepFiles) {
  comm->barrier();
  if (!keepFiles && comm->getRank() == 0) {
    std::remove(filename.c_str());
  }
  comm->barrier();
}

template <class LO, class GO, class Node>
RCP<const Tpetra::Map<LO, GO, Node> >
makeCyclicMap(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
              const global_size_t globalNumElts) {
  Array<GO> gids;
  const global_size_t stride = static_cast<global_size_t>(comm->getSize());
  for (global_size_t gid = static_cast<global_size_t>(comm->getRank());
       gid < globalNumElts;
       gid += stride) {
    gids.push_back(static_cast<GO>(gid));
  }
  return rcp(new Tpetra::Map<LO, GO, Node>(globalNumElts,
                                           gids(),
                                           static_cast<GO>(0),
                                           comm));
}

template <class ST, class LO, class GO, class Node>
RCP<Tpetra::MultiVector<ST, LO, GO, Node> >
makeDenseTestMultiVector(const RCP<const Tpetra::Map<LO, GO, Node> >& map,
                         const size_t numVecs) {
  using multivector_type = Tpetra::MultiVector<ST, LO, GO, Node>;

  auto X          = rcp(new multivector_type(map, numVecs));
  const auto gids = map->getLocalElementList();
  for (size_t j = 0; j < numVecs; ++j) {
    auto data = X->getDataNonConst(j);
    for (size_t i = 0; i < static_cast<size_t>(gids.size()); ++i) {
      data[i] = static_cast<ST>(1000 + 100 * j + static_cast<size_t>(gids[i]));
    }
  }
  return X;
}

template <class MV>
void assertSameMultiVector(const MV& X,
                           const MV& Y) {
  TEUCHOS_TEST_FOR_EXCEPTION(!X.getMap()->isSameAs(*Y.getMap()),
                             std::logic_error,
                             "MultiVector maps differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getGlobalLength() != Y.getGlobalLength(),
                             std::logic_error,
                             "MultiVector global lengths differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(),
                             std::logic_error,
                             "MultiVector column counts differ.");

  for (size_t j = 0; j < X.getNumVectors(); ++j) {
    const auto xData = X.getData(j);
    const auto yData = Y.getData(j);
    TEUCHOS_TEST_FOR_EXCEPTION(xData.size() != yData.size(),
                               std::logic_error,
                               "MultiVector local column lengths differ.");
    for (size_t i = 0; i < static_cast<size_t>(xData.size()); ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION(xData[i] != yData[i],
                                 std::logic_error,
                                 "MultiVector entry mismatch at column " << j
                                                                         << ", local row " << i << ".");
    }
  }
}

template <class ST, class LO, class GO, class Node>
RCP<Tpetra::CrsMatrix<ST, LO, GO, Node> >
makeBandedMatrix(const RCP<const Tpetra::Map<LO, GO, Node> >& map,
                 const int nnzPerRow) {
  using matrix_type = Tpetra::CrsMatrix<ST, LO, GO, Node>;

  const global_size_t numGlobalRows = map->getGlobalNumElements();
  const GO maxColumn                = static_cast<GO>(numGlobalRows - 1);
  const GO halfBandwidth            = static_cast<GO>(nnzPerRow / 2);
  Teuchos::Array<size_t> rowLengths(map->getLocalNumElements());
  for (size_t lclRow = 0; lclRow < map->getLocalNumElements(); ++lclRow) {
    const GO gblRow   = map->getGlobalElement(static_cast<LO>(lclRow));
    const GO firstCol = gblRow > halfBandwidth ? static_cast<GO>(gblRow - halfBandwidth) : static_cast<GO>(0);
    const GO lastCol  = std::min(maxColumn, static_cast<GO>(gblRow + halfBandwidth));
    rowLengths[lclRow] = static_cast<size_t>(lastCol - firstCol + static_cast<GO>(1));
  }
  auto A = rcp(new matrix_type(map, rowLengths()));

  for (size_t lclRow = 0; lclRow < map->getLocalNumElements(); ++lclRow) {
    const GO gblRow   = map->getGlobalElement(static_cast<LO>(lclRow));
    const GO firstCol = gblRow > halfBandwidth ? static_cast<GO>(gblRow - halfBandwidth) : static_cast<GO>(0);
    const GO lastCol  = std::min(maxColumn, static_cast<GO>(gblRow + halfBandwidth));

    Array<GO> cols;
    Array<ST> vals;
    cols.reserve(nnzPerRow);
    vals.reserve(nnzPerRow);
    for (GO col = firstCol; col <= lastCol; ++col) {
      cols.push_back(col);
      vals.push_back(col == gblRow ? static_cast<ST>(2) : static_cast<ST>(-1));
    }
    A->insertGlobalValues(gblRow, cols(), vals());
  }
  A->fillComplete(map, map);
  return A;
}

template <class MatrixType>
void assertSameMatrix(const MatrixType& A,
                      const MatrixType& B) {
  using GO = typename MatrixType::global_ordinal_type;
  using ST = typename MatrixType::scalar_type;
  using LO = typename MatrixType::local_ordinal_type;

  TEUCHOS_TEST_FOR_EXCEPTION(!A.getRowMap()->isSameAs(*B.getRowMap()),
                             std::logic_error,
                             "Matrix row maps differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.getDomainMap()->isSameAs(*B.getDomainMap()),
                             std::logic_error,
                             "Matrix domain maps differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(!A.getRangeMap()->isSameAs(*B.getRangeMap()),
                             std::logic_error,
                             "Matrix range maps differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getGlobalNumRows() != B.getGlobalNumRows(),
                             std::logic_error,
                             "Matrix global row counts differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getGlobalNumCols() != B.getGlobalNumCols(),
                             std::logic_error,
                             "Matrix global column counts differ.");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getGlobalNumEntries() != B.getGlobalNumEntries(),
                             std::logic_error,
                             "Matrix global entry counts differ.");

  const auto rowMap = A.getRowMap();
  for (size_t lclRow = 0; lclRow < rowMap->getLocalNumElements(); ++lclRow) {
    const GO gblRow   = rowMap->getGlobalElement(static_cast<LO>(lclRow));
    const size_t aNum = A.getNumEntriesInGlobalRow(gblRow);
    const size_t bNum = B.getNumEntriesInGlobalRow(gblRow);
    TEUCHOS_TEST_FOR_EXCEPTION(aNum != bNum,
                               std::logic_error,
                               "Matrix row lengths differ for global row " << gblRow << ".");

    typename MatrixType::nonconst_global_inds_host_view_type aInds("aInds", aNum);
    typename MatrixType::nonconst_values_host_view_type aVals("aVals", aNum);
    typename MatrixType::nonconst_global_inds_host_view_type bInds("bInds", bNum);
    typename MatrixType::nonconst_values_host_view_type bVals("bVals", bNum);
    size_t aRead = 0;
    size_t bRead = 0;
    A.getGlobalRowCopy(gblRow, aInds, aVals, aRead);
    B.getGlobalRowCopy(gblRow, bInds, bVals, bRead);
    TEUCHOS_TEST_FOR_EXCEPTION(aRead != bRead,
                               std::logic_error,
                               "Matrix copied row lengths differ for global row " << gblRow << ".");

    std::vector<std::pair<GO, ST> > aEntries;
    std::vector<std::pair<GO, ST> > bEntries;
    aEntries.reserve(aRead);
    bEntries.reserve(bRead);
    for (size_t k = 0; k < aRead; ++k) {
      aEntries.push_back(std::make_pair(aInds(k), aVals(k)));
      bEntries.push_back(std::make_pair(bInds(k), bVals(k)));
    }
    auto byColumn = [](const std::pair<GO, ST>& lhs,
                       const std::pair<GO, ST>& rhs) {
      return lhs.first < rhs.first;
    };
    std::sort(aEntries.begin(), aEntries.end(), byColumn);
    std::sort(bEntries.begin(), bEntries.end(), byColumn);
    for (size_t k = 0; k < aEntries.size(); ++k) {
      TEUCHOS_TEST_FOR_EXCEPTION(aEntries[k].first != bEntries[k].first ||
                                     aEntries[k].second != bEntries[k].second,
                                 std::logic_error,
                                 "Matrix entry mismatch in global row " << gblRow
                                                                        << ", position " << k << ".");
    }
  }
}

struct Config {
  long long globalNumRows;
  int numVectors;
  int numRepeats;
  int nnzPerRow;
  bool runDense;
  bool runSparse;
  bool customMap;
  bool keepFiles;
  std::string filenameStem;
};

void validateConfig(const Config& config,
                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  using map_type = Tpetra::Map<>;
  using GO       = map_type::global_ordinal_type;

  TEUCHOS_TEST_FOR_EXCEPTION(config.globalNumRows <= 0,
                             std::invalid_argument,
                             "--globalNumRows must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(config.numVectors <= 0,
                             std::invalid_argument,
                             "--numVectors must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(config.numRepeats <= 0,
                             std::invalid_argument,
                             "--numRepeats must be positive.");
  TEUCHOS_TEST_FOR_EXCEPTION(config.nnzPerRow <= 0 || config.nnzPerRow % 2 == 0,
                             std::invalid_argument,
                             "--nnzPerRow must be a positive odd integer.");
  TEUCHOS_TEST_FOR_EXCEPTION(!config.runDense && !config.runSparse,
                             std::invalid_argument,
                             "At least one of --dense or --sparse must be enabled.");
  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<unsigned long long>(config.globalNumRows) >
                                 static_cast<unsigned long long>(Teuchos::OrdinalTraits<GO>::max()),
                             std::invalid_argument,
                             "--globalNumRows exceeds the maximum global ordinal for this Tpetra instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<unsigned long long>(config.globalNumRows) >
                                 static_cast<unsigned long long>(std::numeric_limits<global_size_t>::max()),
                             std::invalid_argument,
                             "--globalNumRows exceeds Tpetra::global_size_t.");

  if (comm->getRank() == 0) {
    std::cout << "BinaryIO perf configuration:" << std::endl
              << "  global rows : " << config.globalNumRows << std::endl
              << "  num vectors : " << config.numVectors << std::endl
              << "  repeats     : " << config.numRepeats << std::endl
              << "  nnz/row     : " << config.nnzPerRow << std::endl
              << "  dense       : " << (config.runDense ? "on" : "off") << std::endl
              << "  sparse      : " << (config.runSparse ? "on" : "off") << std::endl
              << "  custom map  : " << (config.customMap ? "on" : "off") << std::endl
              << std::endl;
  }
}

void runDenseCase(const Config& config,
                  const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  using ST               = Tpetra::Details::DefaultTypes::scalar_type;
  using LO               = Tpetra::Map<>::local_ordinal_type;
  using GO               = Tpetra::Map<>::global_ordinal_type;
  using Node             = Tpetra::Map<>::node_type;
  using map_type         = Tpetra::Map<LO, GO, Node>;
  using multivector_type = Tpetra::MultiVector<ST, LO, GO, Node>;
  using import_type      = Tpetra::Import<LO, GO, Node>;
  using binary_io_type   = Tpetra::BinaryIO<ST, LO, GO, Node>;

  const global_size_t globalNumRows = static_cast<global_size_t>(config.globalNumRows);
  auto fileMap                      = rcp(new map_type(globalNumRows,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  RCP<const map_type> targetMap;
  if (config.customMap) {
    targetMap = makeCyclicMap<LO, GO, Node>(comm, globalNumRows);
  } else {
    targetMap = fileMap;
  }
  auto X = makeDenseTestMultiVector<ST, LO, GO, Node>(fileMap,
                                                      static_cast<size_t>(config.numVectors));
  const std::string filename = makeFilename(config.filenameStem, *comm, "_dense.bin");

  for (int repeat = 0; repeat < config.numRepeats; ++repeat) {
    comm->barrier();
    {
      auto timer = Teuchos::TimeMonitor::getNewTimer("BinaryIO dense write");
      Teuchos::TimeMonitor timeGuard(*timer);
      binary_io_type::writeDenseFile(filename, *X);
    }

    comm->barrier();
    RCP<multivector_type> Y;
    {
      auto timer = Teuchos::TimeMonitor::getNewTimer(config.customMap ? "BinaryIO dense read custom map" : "BinaryIO dense read default map");
      Teuchos::TimeMonitor timeGuard(*timer);
      Y = config.customMap ? binary_io_type::readDenseFile(filename, targetMap) : binary_io_type::readDenseFile(filename, comm);
    }

    if (config.customMap) {
      auto expected = rcp(new multivector_type(targetMap, X->getNumVectors()));
      import_type importer(fileMap, targetMap);
      expected->doImport(*X, importer, Tpetra::INSERT);
      assertSameMultiVector(*expected, *Y);
    } else {
      assertSameMultiVector(*X, *Y);
    }

    cleanupFile(filename, comm, config.keepFiles);
  }
}

void runSparseCase(const Config& config,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  using ST             = Tpetra::Details::DefaultTypes::scalar_type;
  using LO             = Tpetra::Map<>::local_ordinal_type;
  using GO             = Tpetra::Map<>::global_ordinal_type;
  using Node           = Tpetra::Map<>::node_type;
  using map_type       = Tpetra::Map<LO, GO, Node>;
  using matrix_type    = Tpetra::CrsMatrix<ST, LO, GO, Node>;
  using import_type    = Tpetra::Import<LO, GO, Node>;
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, Node>;

  const global_size_t globalNumRows = static_cast<global_size_t>(config.globalNumRows);
  auto fileMap                      = rcp(new map_type(globalNumRows,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  RCP<const map_type> targetMap;
  if (config.customMap) {
    targetMap = makeCyclicMap<LO, GO, Node>(comm, globalNumRows);
  } else {
    targetMap = fileMap;
  }
  auto A = makeBandedMatrix<ST, LO, GO, Node>(fileMap, config.nnzPerRow);
  const std::string filename = makeFilename(config.filenameStem, *comm, "_sparse.bin");

  for (int repeat = 0; repeat < config.numRepeats; ++repeat) {
    comm->barrier();
    {
      auto timer = Teuchos::TimeMonitor::getNewTimer("BinaryIO sparse write");
      Teuchos::TimeMonitor timeGuard(*timer);
      binary_io_type::writeSparseFile(filename, *A);
    }

    comm->barrier();
    RCP<matrix_type> B;
    {
      auto timer = Teuchos::TimeMonitor::getNewTimer(config.customMap ? "BinaryIO sparse read custom map" : "BinaryIO sparse read default map");
      Teuchos::TimeMonitor timeGuard(*timer);
      B = config.customMap ? binary_io_type::readSparseFile(filename, targetMap, targetMap, targetMap, true) : binary_io_type::readSparseFile(filename, comm);
    }

    if (config.customMap) {
      auto expected = rcp(new matrix_type(targetMap, static_cast<size_t>(config.nnzPerRow)));
      import_type importer(fileMap, targetMap);
      expected->doImport(*A, importer, Tpetra::INSERT);
      expected->fillComplete(targetMap, targetMap);
      assertSameMatrix(*expected, *B);
    } else {
      assertSameMatrix(*A, *B);
    }

    cleanupFile(filename, comm, config.keepFiles);
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scope(&argc, &argv);
  const auto comm = Tpetra::getDefaultComm();

  Config config;
  config.globalNumRows = 0;
  config.numVectors    = 4;
  config.numRepeats    = 3;
  config.nnzPerRow     = 5;
  config.runDense      = true;
  config.runSparse     = true;
  config.customMap     = false;
  config.keepFiles     = false;
  config.filenameStem  = "Tpetra_BinaryIO_PerfTest";

  Teuchos::CommandLineProcessor cmdp(false, true);
  cmdp.setOption("globalNumRows", &config.globalNumRows,
                 "Global number of rows to generate.  If zero, choose a default based on MPI size.");
  cmdp.setOption("numVectors", &config.numVectors,
                 "Number of dense vectors to write and read.");
  cmdp.setOption("numRepeats", &config.numRepeats,
                 "Number of write/read repetitions for each enabled object type.");
  cmdp.setOption("nnzPerRow", &config.nnzPerRow,
                 "Odd number of sparse entries per row in the generated banded matrix.");
  cmdp.setOption("dense", "no-dense", &config.runDense,
                 "Enable or disable dense BinaryIO timing.");
  cmdp.setOption("sparse", "no-sparse", &config.runSparse,
                 "Enable or disable sparse BinaryIO timing.");
  cmdp.setOption("custom-map", "default-map", &config.customMap,
                 "Read back into a cyclic target map instead of the file distribution.");
  cmdp.setOption("keep-files", "delete-files", &config.keepFiles,
                 "Keep generated binary files instead of deleting them after each repetition.");
  cmdp.setOption("filename-stem", &config.filenameStem,
                 "Base filename stem for generated binary files.");

  const Teuchos::CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse(argc, argv);
  if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    if (comm->getRank() == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    return EXIT_SUCCESS;
  }
  if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  if (config.globalNumRows == 0) {
    config.globalNumRows = 10000ll * static_cast<long long>(comm->getSize());
  }

  int ierr = 0;
  Teuchos::RCP<Teuchos::StackedTimer> stackedTimer =
      rcp(new Teuchos::StackedTimer("BinaryIO perf"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

  try {
    validateConfig(config, comm);

    if (config.runDense) {
      runDenseCase(config, comm);
    }
    if (config.runSparse) {
      runSparseCase(config, comm);
    }
  } catch (std::exception& e) {
    ierr = 1;
    if (comm->getRank() == 0) {
      std::cout << e.what() << std::endl;
    }
  }

  stackedTimer->stopBaseTimer();

  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction  = true;
  options.output_histogram = true;
  options.output_minmax    = true;

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);
  stackedTimer->report(out, comm, options);

  if (comm->getRank() == 0) {
    if (ierr == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  return ierr;
}
