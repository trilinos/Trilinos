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
#include "Tpetra_Import.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#include <algorithm>
#include <cstdio>
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
                         const Teuchos::Comm<int>& comm) {
  std::ostringstream os;
  os << stem << "_p" << comm.getSize() << ".bin";
  return os.str();
}

void cleanupFile(const std::string& filename,
                 const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  comm->barrier();
  if (comm->getRank() == 0) {
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
makeTridiagonalMatrix(const RCP<const Tpetra::Map<LO, GO, Node> >& map) {
  using matrix_type = Tpetra::CrsMatrix<ST, LO, GO, Node>;

  const global_size_t numGlobalRows = map->getGlobalNumElements();
  Teuchos::Array<size_t> rowLengths(map->getLocalNumElements(), static_cast<size_t>(3));
  const auto minAllGlobalIndex = map->getMinAllGlobalIndex();
  const auto maxAllGlobalIndex = map->getMaxAllGlobalIndex();
  if (map->getMinGlobalIndex() == minAllGlobalIndex) {
    rowLengths[0] = static_cast<size_t>(2);
  }
  if (map->getMaxGlobalIndex() == maxAllGlobalIndex) {
    rowLengths[map->getLocalNumElements() - 1] = static_cast<size_t>(2);
  }

  auto A = rcp(new matrix_type(map, rowLengths()));
  for (size_t lclRow = 0; lclRow < map->getLocalNumElements(); ++lclRow) {
    const GO gblRow = map->getGlobalElement(static_cast<LO>(lclRow));
    Array<GO> cols;
    Array<ST> vals;
    if (gblRow > static_cast<GO>(0)) {
      cols.push_back(static_cast<GO>(gblRow - 1));
      vals.push_back(static_cast<ST>(-1));
    }
    cols.push_back(gblRow);
    vals.push_back(static_cast<ST>(2 + (static_cast<size_t>(gblRow) % 7)));
    if (static_cast<global_size_t>(gblRow + 1) < numGlobalRows) {
      cols.push_back(static_cast<GO>(gblRow + 1));
      vals.push_back(static_cast<ST>(-1));
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
    const GO gblRow   = rowMap->getGlobalElement(static_cast<typename MatrixType::local_ordinal_type>(lclRow));
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
    TEUCHOS_TEST_FOR_EXCEPTION(aEntries.size() != bEntries.size(),
                               std::logic_error,
                               "Matrix sorted row lengths differ for global row " << gblRow << ".");
    for (size_t k = 0; k < aEntries.size(); ++k) {
      TEUCHOS_TEST_FOR_EXCEPTION(aEntries[k].first != bEntries[k].first ||
                                     aEntries[k].second != bEntries[k].second,
                                 std::logic_error,
                                 "Matrix entry mismatch in global row " << gblRow
                                                                        << ", position " << k << ".");
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, MapRoundTrip,
                                  ST, LO, GO, NODE) {
  using map_type       = Tpetra::Map<LO, GO, NODE>;
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, NODE>;

  auto comm                         = Tpetra::getDefaultComm();
  const global_size_t globalNumElts = 11 * static_cast<global_size_t>(comm->getSize());
  auto map                          = rcp(new map_type(globalNumElts,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  const std::string filename        = makeFilename("Tpetra_BinaryIO_MapRoundTrip", *comm);

  binary_io_type::writeMapFile(filename, *map);
  auto inMap = binary_io_type::readMapFile(filename, comm);

  TEST_ASSERT(map->isSameAs(*inMap));
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, DenseRoundTripDefaultMap,
                                  ST, LO, GO, NODE) {
  using map_type         = Tpetra::Map<LO, GO, NODE>;
  using multivector_type = Tpetra::MultiVector<ST, LO, GO, NODE>;
  using binary_io_type   = Tpetra::BinaryIO<ST, LO, GO, NODE>;

  auto comm                         = Tpetra::getDefaultComm();
  const global_size_t globalNumElts = 13 * static_cast<global_size_t>(comm->getSize());
  auto map                          = rcp(new map_type(globalNumElts,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  auto X                            = makeDenseTestMultiVector<ST, LO, GO, NODE>(map, 3);
  const std::string filename        = makeFilename("Tpetra_BinaryIO_DenseRoundTripDefaultMap", *comm);

  binary_io_type::writeDenseFile(filename, *X);
  auto Y = binary_io_type::readDenseFile(filename, comm);

  assertSameMultiVector(*X, *Y);
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, DenseRoundTripCustomMap,
                                  ST, LO, GO, NODE) {
  using map_type         = Tpetra::Map<LO, GO, NODE>;
  using multivector_type = Tpetra::MultiVector<ST, LO, GO, NODE>;
  using binary_io_type   = Tpetra::BinaryIO<ST, LO, GO, NODE>;
  using import_type      = Tpetra::Import<LO, GO, NODE>;

  auto comm                         = Tpetra::getDefaultComm();
  const global_size_t globalNumElts = 13 * static_cast<global_size_t>(comm->getSize());
  auto fileMap                      = rcp(new map_type(globalNumElts,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  auto targetMap                    = makeCyclicMap<LO, GO, NODE>(comm, globalNumElts);
  auto X                            = makeDenseTestMultiVector<ST, LO, GO, NODE>(fileMap, 2);
  const std::string filename        = makeFilename("Tpetra_BinaryIO_DenseRoundTripCustomMap", *comm);

  binary_io_type::writeDenseFile(filename, *X);
  auto Y = binary_io_type::readDenseFile(filename, targetMap);

  auto expected = rcp(new multivector_type(targetMap, X->getNumVectors()));
  import_type importer(fileMap, targetMap);
  expected->doImport(*X, importer, Tpetra::INSERT);

  assertSameMultiVector(*expected, *Y);
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, SparseRoundTripDefaultMap,
                                  ST, LO, GO, NODE) {
  using map_type       = Tpetra::Map<LO, GO, NODE>;
  using matrix_type    = Tpetra::CrsMatrix<ST, LO, GO, NODE>;
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, NODE>;

  auto comm                         = Tpetra::getDefaultComm();
  const global_size_t globalNumRows = 9 * static_cast<global_size_t>(comm->getSize()) + 1;
  auto map                          = rcp(new map_type(globalNumRows,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  auto A                            = makeTridiagonalMatrix<ST, LO, GO, NODE>(map);
  const std::string filename        = makeFilename("Tpetra_BinaryIO_SparseRoundTripDefaultMap", *comm);

  binary_io_type::writeSparseFile(filename, *A);
  auto B = binary_io_type::readSparseFile(filename, comm);

  assertSameMatrix(*A, *B);
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, SparseRoundTripCustomRowMap,
                                  ST, LO, GO, NODE) {
  using map_type       = Tpetra::Map<LO, GO, NODE>;
  using matrix_type    = Tpetra::CrsMatrix<ST, LO, GO, NODE>;
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, NODE>;
  using import_type    = Tpetra::Import<LO, GO, NODE>;

  auto comm                         = Tpetra::getDefaultComm();
  const global_size_t globalNumRows = 9 * static_cast<global_size_t>(comm->getSize()) + 1;
  auto fileMap                      = rcp(new map_type(globalNumRows,
                                                       static_cast<GO>(0),
                                                       comm,
                                                       Tpetra::GloballyDistributed));
  auto targetMap                    = makeCyclicMap<LO, GO, NODE>(comm, globalNumRows);
  auto A                            = makeTridiagonalMatrix<ST, LO, GO, NODE>(fileMap);
  const std::string filename        = makeFilename("Tpetra_BinaryIO_SparseRoundTripCustomRowMap", *comm);

  binary_io_type::writeSparseFile(filename, *A);
  auto B = binary_io_type::readSparseFile(filename, targetMap, targetMap, targetMap, true);

  auto expected = rcp(new matrix_type(targetMap, 3));
  import_type importer(fileMap, targetMap);
  expected->doImport(*A, importer, Tpetra::INSERT);
  expected->fillComplete(targetMap, targetMap);

  assertSameMatrix(*expected, *B);
  cleanupFile(filename, comm);
}

#if defined(HAVE_TPETRA_INST_DOUBLE)
#define UNIT_TEST_GROUP(LO, GO, NODE)                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, MapRoundTrip, double, LO, GO, NODE)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, DenseRoundTripDefaultMap, double, LO, GO, NODE)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, DenseRoundTripCustomMap, double, LO, GO, NODE)   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, SparseRoundTripDefaultMap, double, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, SparseRoundTripCustomRowMap, double, LO, GO, NODE)
#else
#define UNIT_TEST_GROUP(LO, GO, NODE)
#endif

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

}  // namespace
