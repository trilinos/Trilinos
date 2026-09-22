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
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Core.hpp"
#include "TpetraCore_ETIHelperMacros.h"

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

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
  using map_type    = Tpetra::Map<LO, GO, Node>;
  using device_type = typename map_type::device_type;

  const global_size_t stride = static_cast<global_size_t>(comm->getSize());
  const global_size_t rank   = static_cast<global_size_t>(comm->getRank());
  const size_t localNumElts  = rank < globalNumElts ? static_cast<size_t>((globalNumElts - rank + stride - 1) / stride) : 0;
  Kokkos::View<GO*, device_type> gids("Tpetra_BinaryIO_InOutTest::cyclicGids", localNumElts);
  Kokkos::parallel_for(
      "Tpetra_BinaryIO_InOutTest::fillCyclicGids",
      Kokkos::RangePolicy<typename device_type::execution_space>(0, localNumElts),
      KOKKOS_LAMBDA(const size_t i) {
        gids(i) = static_cast<GO>(rank + static_cast<global_size_t>(i) * stride);
      });
  return rcp(new map_type(globalNumElts,
                          gids,
                          static_cast<GO>(0),
                          comm));
}

template <class LO, class GO, class Node>
RCP<const Tpetra::Map<LO, GO, Node> >
makeImbalancedContiguousMap(const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {
  const size_t numLocalElts         = (comm->getRank() == 0) ? static_cast<size_t>(comm->getSize() + 1) : static_cast<size_t>(1);
  const global_size_t globalNumElts = static_cast<global_size_t>(2 * comm->getSize());
  return rcp(new Tpetra::Map<LO, GO, Node>(globalNumElts,
                                           numLocalElts,
                                           static_cast<GO>(0),
                                           comm));
}

template <class ST, class LO, class GO, class Node>
RCP<Tpetra::MultiVector<ST, LO, GO, Node> >
makeDenseTestMultiVector(const RCP<const Tpetra::Map<LO, GO, Node> >& map,
                         const size_t numVecs) {
  using multivector_type = Tpetra::MultiVector<ST, LO, GO, Node>;
  using impl_scalar_type = typename multivector_type::impl_scalar_type;
  using device_type      = typename multivector_type::device_type;
  using execution_space  = typename device_type::execution_space;

  auto X            = rcp(new multivector_type(map, numVecs));
  auto localX       = X->getLocalViewDevice(Tpetra::Access::OverwriteAll);
  auto localMap     = map->getLocalMap();
  const size_t nvec = numVecs;
  Kokkos::parallel_for(
      "Tpetra_BinaryIO_InOutTest::fillDense",
      Kokkos::MDRangePolicy<execution_space, Kokkos::Rank<2>>({0, 0}, {localX.extent(0), nvec}),
      KOKKOS_LAMBDA(const size_t i, const size_t j) {
        const auto gid = localMap.getGlobalElement(static_cast<LO>(i));
        localX(i, j)  = static_cast<impl_scalar_type>(1000 + 100 * j + static_cast<size_t>(gid));
      });
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
  using matrix_type          = Tpetra::CrsMatrix<ST, LO, GO, Node>;
  using local_graph_type     = typename matrix_type::local_graph_device_type;
  using rowptr_type          = typename local_graph_type::row_map_type::non_const_type;
  using rowptr_value_type    = typename rowptr_type::non_const_value_type;
  using colidx_type          = typename local_graph_type::entries_type::non_const_type;
  using values_type          = typename matrix_type::local_matrix_device_type::values_type::non_const_type;
  using impl_scalar_type     = typename matrix_type::impl_scalar_type;
  using device_type          = typename matrix_type::device_type;
  using memory_space         = typename device_type::memory_space;
  using execution_space      = typename device_type::execution_space;

  const size_t localNumRows    = map->getLocalNumElements();
  const auto minAllGlobalIndex = map->getMinAllGlobalIndex();
  const auto maxAllGlobalIndex = map->getMaxAllGlobalIndex();
  const auto localRowMap       = map->getLocalMap();

  rowptr_type rowPtr("Tpetra_BinaryIO_InOutTest::rowPtr", localNumRows + 1);
  Kokkos::deep_copy(rowPtr, static_cast<rowptr_value_type>(0));
  Kokkos::parallel_scan(
      "Tpetra_BinaryIO_InOutTest::countRows",
      Kokkos::RangePolicy<execution_space>(0, localNumRows),
      KOKKOS_LAMBDA(const size_t lclRow, rowptr_value_type& offset, const bool finalPass) {
        const GO gblRow = localRowMap.getGlobalElement(static_cast<LO>(lclRow));
        if (finalPass) {
          rowPtr(lclRow) = offset;
        }
        offset += static_cast<rowptr_value_type>(1 + (gblRow > minAllGlobalIndex ? 1 : 0) + (gblRow < maxAllGlobalIndex ? 1 : 0));
        if (finalPass && lclRow + 1 == localNumRows) {
          rowPtr(localNumRows) = offset;
        }
      });
  auto rowPtrHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowPtr);
  const size_t localNnz = rowPtrHost(localNumRows);

  Kokkos::View<GO*, memory_space> globalColumns("Tpetra_BinaryIO_InOutTest::globalColumns", localNnz);
  Kokkos::parallel_for(
      "Tpetra_BinaryIO_InOutTest::fillGlobalColumns",
      Kokkos::RangePolicy<execution_space>(0, localNumRows),
      KOKKOS_LAMBDA(const size_t lclRow) {
        const GO gblRow = localRowMap.getGlobalElement(static_cast<LO>(lclRow));
        size_t offset   = rowPtr(lclRow);
        if (gblRow > minAllGlobalIndex) {
          globalColumns(offset++) = static_cast<GO>(gblRow - 1);
        }
        globalColumns(offset++) = gblRow;
        if (gblRow < maxAllGlobalIndex) {
          globalColumns(offset) = static_cast<GO>(gblRow + 1);
        }
      });

  RCP<const Tpetra::Map<LO, GO, Node> > colMap;
  std::ostringstream errStrm;
  const int err = Tpetra::Details::makeColMap<LO, GO, Node>(colMap,
                                                            map,
                                                            globalColumns,
                                                            &errStrm);
  TEUCHOS_TEST_FOR_EXCEPTION(err != 0 || colMap.is_null(),
                             std::runtime_error,
                             "Failed to construct test column map. " << errStrm.str());
  const auto localColMap = colMap->getLocalMap();

  colidx_type colInd("Tpetra_BinaryIO_InOutTest::colInd", localNnz);
  values_type values("Tpetra_BinaryIO_InOutTest::values", localNnz);
  Kokkos::parallel_for(
      "Tpetra_BinaryIO_InOutTest::fillLocalMatrix",
      Kokkos::RangePolicy<execution_space>(0, localNumRows),
      KOKKOS_LAMBDA(const size_t lclRow) {
        const GO gblRow = localRowMap.getGlobalElement(static_cast<LO>(lclRow));
        size_t offset   = rowPtr(lclRow);
        if (gblRow > minAllGlobalIndex) {
          colInd(offset) = localColMap.getLocalElement(static_cast<GO>(gblRow - 1));
          values(offset) = static_cast<impl_scalar_type>(-1);
          ++offset;
        }
        colInd(offset) = localColMap.getLocalElement(gblRow);
        values(offset) = static_cast<impl_scalar_type>(2 + (static_cast<size_t>(gblRow) % 7));
        ++offset;
        if (gblRow < maxAllGlobalIndex) {
          colInd(offset) = localColMap.getLocalElement(static_cast<GO>(gblRow + 1));
          values(offset) = static_cast<impl_scalar_type>(-1);
        }
      });
  Tpetra::Import_Util::sortCrsEntries(rowPtr, colInd, values);

  auto A = rcp(new matrix_type(map, colMap, rowPtr, colInd, values));
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
  auto inMap = Tpetra::readBinaryMapFile<LO, GO, NODE>(filename, comm);

  TEST_ASSERT(map->isSameAs(*inMap));
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, MapRoundTripImbalancedContiguous,
                                  ST, LO, GO, NODE) {
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, NODE>;

  auto comm                  = Tpetra::getDefaultComm();
  auto map                   = makeImbalancedContiguousMap<LO, GO, NODE>(comm);
  const std::string filename = makeFilename("Tpetra_BinaryIO_MapRoundTripImbalancedContiguous", *comm);

  binary_io_type::writeMapFile(filename, *map);
  auto inMap = binary_io_type::readMapFile(filename, comm);

  TEST_ASSERT(map->isSameAs(*inMap));
  cleanupFile(filename, comm);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BinaryIO, DenseRoundTripDefaultMap,
                                  ST, LO, GO, NODE) {
  using map_type       = Tpetra::Map<LO, GO, NODE>;
  using binary_io_type = Tpetra::BinaryIO<ST, LO, GO, NODE>;

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

  import_type importer(fileMap, targetMap);
  auto expected = Tpetra::importAndFillCompleteCrsMatrix<matrix_type>(A, importer, targetMap, targetMap);

  assertSameMatrix(*expected, *B);
  cleanupFile(filename, comm);
}

#if defined(HAVE_TPETRA_INST_DOUBLE)
#define UNIT_TEST_GROUP(LO, GO, NODE)                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, MapRoundTrip, double, LO, GO, NODE)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, MapRoundTripImbalancedContiguous, double, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, DenseRoundTripDefaultMap, double, LO, GO, NODE)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, DenseRoundTripCustomMap, double, LO, GO, NODE)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, SparseRoundTripDefaultMap, double, LO, GO, NODE)        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BinaryIO, SparseRoundTripCustomRowMap, double, LO, GO, NODE)
#else
#define UNIT_TEST_GROUP(LO, GO, NODE)
#endif

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

}  // namespace
