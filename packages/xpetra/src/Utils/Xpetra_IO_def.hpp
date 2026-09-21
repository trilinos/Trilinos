// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_IO_decl.hpp"
#include <fstream>
#include "Xpetra_ConfigDefs.hpp"

#include <Teuchos_CommHelpers.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_BinaryIO.hpp>
#include <Tpetra_Details_OrdinalTraits.hpp>
#include <Tpetra_Details_makeColMap.hpp>
#include <Tpetra_Import_Util2.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraBlockCrsMatrix.hpp>
#include "Tpetra_Util.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixMatrix.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_BlockedCrsMatrix.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MatrixFactory.hpp"

#include <Teuchos_MatrixMarket_Raw_Writer.hpp>
#include <Teuchos_MatrixMarket_Raw_Reader.hpp>
#include <algorithm>
#include <complex>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

namespace Xpetra {

namespace Details {

template <class Scalar>
struct binaryIOAvailableForScalar : std::false_type {};

#if defined(HAVE_TPETRA_INST_FLOAT)
template <>
struct binaryIOAvailableForScalar<float> : std::true_type {};
#endif

#if defined(HAVE_TPETRA_INST_DOUBLE)
template <>
struct binaryIOAvailableForScalar<double> : std::true_type {};
#endif

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
template <>
struct binaryIOAvailableForScalar<std::complex<float>> : std::true_type {};
#endif

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE)
template <>
struct binaryIOAvailableForScalar<std::complex<double>> : std::true_type {};
#endif

#if defined(HAVE_TPETRA_INST_FLOAT128)
template <>
struct binaryIOAvailableForScalar<__float128> : std::true_type {};
#endif

template <class LocalOrdinal>
struct binaryIOAvailableForLocalOrdinalScalar : std::false_type {};

#if !defined(HAVE_TPETRA_REDUCED_ETI) && !defined(HAVE_TPETRA_INST_INT_INT)
template <>
struct binaryIOAvailableForLocalOrdinalScalar<int> : std::true_type {};
#endif

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> readBinaryMap(const std::string& fileName,
                                                                        const RCP<const Teuchos::Comm<int>>& comm) {
  return Tpetra::readBinaryMapFile<LocalOrdinal, GlobalOrdinal, Node>(fileName, comm);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const RCP<const Teuchos::Comm<int>>& comm);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string&,
                     const RCP<const Teuchos::Comm<int>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: binary Tpetra matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
                     const bool callFillComplete);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string&,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const bool) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: binary Tpetra matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

inline bool hasTpetraBinaryHeader(const std::string& fileName,
                                  const RCP<const Teuchos::Comm<int>>& comm) {
  int openOk    = 1;
  int hasHeader = 0;
  if (comm->getRank() == 0) {
    std::ifstream in(fileName.c_str(), std::ios::binary);
    openOk = in.good() ? 1 : 0;
    if (openOk) {
      char magic[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      in.read(magic, sizeof(magic));
      const char expectedMagic[8] = {'T', 'P', 'B', 'I', 'O', '0', '0', '1'};
      hasHeader                   = in.good() && std::memcmp(magic, expectedMagic, sizeof(expectedMagic)) == 0 ? 1 : 0;
    }
  }
  Teuchos::broadcast(*comm, 0, 1, &openOk);
  Teuchos::broadcast(*comm, 0, 1, &hasHeader);
  TEUCHOS_TEST_FOR_EXCEPTION(openOk == 0,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Failed to open binary file '" << fileName << "'.");
  return hasHeader != 0;
}

struct LegacyBinaryHeader {
  int numRows;
  int numCols;
  int numEntries;
};

inline unsigned long long legacyBinaryFileSize(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::binary | std::ios::ate);
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Failed to open legacy binary file '" << fileName << "'.");
  const std::streamoff size = in.tellg();
  TEUCHOS_TEST_FOR_EXCEPTION(size < 0,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Failed to determine the size of legacy binary file '" << fileName << "'.");
  return static_cast<unsigned long long>(size);
}

inline LegacyBinaryHeader readAndValidateLegacyBinaryHeader(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::binary);
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Failed to open legacy binary file '" << fileName << "'.");

  LegacyBinaryHeader header = {0, 0, 0};
  in.read(reinterpret_cast<char*>(&header.numRows), sizeof(header.numRows));
  in.read(reinterpret_cast<char*>(&header.numCols), sizeof(header.numCols));
  in.read(reinterpret_cast<char*>(&header.numEntries), sizeof(header.numEntries));
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             Exceptions::RuntimeError,
                             "Xpetra::IO: File '" << fileName << "' is not a valid legacy Xpetra binary matrix file.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.numRows < 0 || header.numCols < 0 || header.numEntries < 0,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Legacy binary file '" << fileName << "' has negative dimensions or entry count.");

  const unsigned long long fileSize = legacyBinaryFileSize(fileName);
  const unsigned long long minSize  = static_cast<unsigned long long>(3 * sizeof(int)) +
                                     static_cast<unsigned long long>(header.numRows) * static_cast<unsigned long long>(2 * sizeof(int)) +
                                     static_cast<unsigned long long>(header.numEntries) * static_cast<unsigned long long>(sizeof(int) + sizeof(double));
  TEUCHOS_TEST_FOR_EXCEPTION(fileSize != minSize,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: File '" << fileName << "' is not a valid legacy Xpetra binary matrix file."
                                                  << " Expected " << minSize << " bytes from its header, but file has " << fileSize << " bytes.");
  return header;
}

template <class T>
void readLegacyBinaryValue(std::ifstream& in, T& value, const std::string& fileName, const char label[]) {
  in.read(reinterpret_cast<char*>(&value), sizeof(T));
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Failed to read " << label << " from legacy binary file '" << fileName << "'.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readLegacyBinarySparseFile(const std::string& oldFileName,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMapInput,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& colMapInput,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMapInput,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMapInput,
                           const RCP<const Teuchos::Comm<int>>& comm,
                           const bool callFillComplete) {
  using map_type             = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using matrix_type          = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using import_type          = Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  using local_graph_type     = typename matrix_type::local_graph_device_type;
  using rowptr_type          = typename local_graph_type::row_map_type::non_const_type;
  using colidx_type          = typename local_graph_type::entries_type::non_const_type;
  using values_type          = typename matrix_type::local_matrix_device_type::values_type::non_const_type;
  using impl_scalar_type     = typename matrix_type::impl_scalar_type;
  using device_type          = typename matrix_type::device_type;
  using execution_space      = typename device_type::execution_space;
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;

  const LegacyBinaryHeader legacyHeader = readAndValidateLegacyBinaryHeader(oldFileName);
  const int m                           = legacyHeader.numRows;
  const int n                           = legacyHeader.numCols;
  const int nnz                         = legacyHeader.numEntries;
  const int myRank                      = comm->getRank();

  const size_t fileLocalNumRows = (myRank == 0) ? static_cast<size_t>(m) : static_cast<size_t>(0);
  const size_t fileLocalNumCols = (myRank == 0) ? static_cast<size_t>(n) : static_cast<size_t>(0);
  const auto fileRowMap         = Teuchos::rcp(new map_type(static_cast<Tpetra::global_size_t>(m), fileLocalNumRows, static_cast<GlobalOrdinal>(0), comm));
  const auto fileColMap         = Teuchos::rcp(new map_type(static_cast<Tpetra::global_size_t>(n), fileLocalNumCols, static_cast<GlobalOrdinal>(0), comm));
  const auto fileDomainMap      = fileColMap;
  const auto fileRangeMap       = fileRowMap;

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> rowLengths("Xpetra::IO::legacyRowLengths", fileLocalNumRows);
  Kokkos::deep_copy(rowLengths, 0ull);
  unsigned long long entriesRead = 0;

  std::ifstream in;
  if (myRank == 0) {
    in.open(oldFileName.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               Exceptions::RuntimeError,
                               "Xpetra::IO: Failed to open legacy binary file '" << oldFileName << "'.");
    in.seekg(static_cast<std::streamoff>(3 * sizeof(int)), std::ios::beg);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               Exceptions::RuntimeError,
                               "Xpetra::IO: Failed to seek past the header of legacy binary file '" << oldFileName << "'.");

    for (int record = 0; record < m; ++record) {
      int row    = 0;
      int rownnz = 0;
      readLegacyBinaryValue(in, row, oldFileName, "row index");
      readLegacyBinaryValue(in, rownnz, oldFileName, "row entry count");
      TEUCHOS_TEST_FOR_EXCEPTION(row < 0 || row >= m || rownnz < 0,
                                 Exceptions::RuntimeError,
                                 "Xpetra::IO: Legacy binary file '" << oldFileName << "' has an invalid row record.");
      entriesRead += static_cast<unsigned long long>(rownnz);
      rowLengths(static_cast<size_t>(row)) += static_cast<unsigned long long>(rownnz);
      for (int j = 0; j < rownnz; ++j) {
        int col = 0;
        readLegacyBinaryValue(in, col, oldFileName, "column index");
        TEUCHOS_TEST_FOR_EXCEPTION(col < 0 || col >= n,
                                   Exceptions::RuntimeError,
                                   "Xpetra::IO: Legacy binary file '" << oldFileName << "' has a column index outside [0, n).");
      }
      for (int j = 0; j < rownnz; ++j) {
        double value = 0.0;
        readLegacyBinaryValue(in, value, oldFileName, "matrix value");
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(entriesRead != static_cast<unsigned long long>(nnz),
                               Exceptions::RuntimeError,
                               "Xpetra::IO: Legacy binary file '" << oldFileName << "' header entry count does not match row records.");
    TEUCHOS_TEST_FOR_EXCEPTION(entriesRead > static_cast<unsigned long long>(std::numeric_limits<size_t>::max()),
                               Exceptions::RuntimeError,
                               "Xpetra::IO: Legacy binary file '" << oldFileName << "' has too many entries for this platform.");
  }

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> rowPtrHost("Xpetra::IO::legacyRowPtr", fileLocalNumRows + 1);
  unsigned long long localNnz64 = 0;
  Kokkos::parallel_scan(
      "Xpetra::IO::legacyRowPtrScan",
      Kokkos::RangePolicy<host_execution_space>(0, fileLocalNumRows),
      KOKKOS_LAMBDA(const size_t row, unsigned long long& count, const bool finalPass) {
        if (finalPass) {
          rowPtrHost(row) = count;
        }
        count += rowLengths(row);
        if (finalPass && row + 1 == fileLocalNumRows) {
          rowPtrHost(row + 1) = count;
        }
      },
      localNnz64);
  if (fileLocalNumRows == 0) {
    rowPtrHost(0) = 0;
  }
  const size_t localNnz = static_cast<size_t>(localNnz64);
  Kokkos::View<unsigned long long*, Kokkos::HostSpace> nextPtr("Xpetra::IO::legacyNextPtr", fileLocalNumRows);
  Kokkos::parallel_for(
      "Xpetra::IO::legacyInitNextPtr",
      Kokkos::RangePolicy<host_execution_space>(0, fileLocalNumRows),
      KOKKOS_LAMBDA(const size_t row) {
        nextPtr(row) = rowPtrHost(row);
      });

  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> globalColumnsHost("Xpetra::IO::legacyGlobalColumns", localNnz);
  Kokkos::View<Scalar*, Kokkos::HostSpace> valuesHost("Xpetra::IO::legacyValues", localNnz);

  if (myRank == 0) {
    in.clear();
    in.seekg(static_cast<std::streamoff>(3 * sizeof(int)), std::ios::beg);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               Exceptions::RuntimeError,
                               "Xpetra::IO: Failed to rewind legacy binary file '" << oldFileName << "'.");

    for (int record = 0; record < m; ++record) {
      int row    = 0;
      int rownnz = 0;
      readLegacyBinaryValue(in, row, oldFileName, "row index");
      readLegacyBinaryValue(in, rownnz, oldFileName, "row entry count");

      Kokkos::View<int*, Kokkos::HostSpace> columns("Xpetra::IO::legacyRowColumns", static_cast<size_t>(rownnz));
      for (int j = 0; j < rownnz; ++j) {
        readLegacyBinaryValue(in, columns(static_cast<size_t>(j)), oldFileName, "column index");
      }

      size_t offset = static_cast<size_t>(nextPtr(static_cast<size_t>(row)));
      for (int j = 0; j < rownnz; ++j) {
        const int col                                      = columns(static_cast<size_t>(j));
        globalColumnsHost(offset + static_cast<size_t>(j)) = static_cast<GlobalOrdinal>(col);
      }
      for (int j = 0; j < rownnz; ++j) {
        double value = 0.0;
        readLegacyBinaryValue(in, value, oldFileName, "matrix value");
        valuesHost(offset + static_cast<size_t>(j)) = static_cast<Scalar>(value);
      }
      nextPtr(static_cast<size_t>(row)) += static_cast<unsigned long long>(rownnz);
    }
  }

  rowptr_type rowPtrDevice("Xpetra::IO::legacyRowPtrDevice", fileLocalNumRows + 1);
  auto rowPtrDeviceHost = Kokkos::create_mirror_view(rowPtrDevice);
  Kokkos::parallel_for(
      "Xpetra::IO::legacyCopyRowPtrToDeviceHost",
      Kokkos::RangePolicy<host_execution_space>(0, fileLocalNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        rowPtrDeviceHost(i) = static_cast<typename rowptr_type::non_const_value_type>(rowPtrHost(i));
      });
  Kokkos::deep_copy(rowPtrDevice, rowPtrDeviceHost);

  Kokkos::View<GlobalOrdinal*, device_type> globalColumnsDevice("Xpetra::IO::legacyGlobalColumnsDevice", localNnz);
  Kokkos::deep_copy(globalColumnsDevice, globalColumnsHost);

  colidx_type localColumnsDevice("Xpetra::IO::legacyLocalColumns", localNnz);
  const auto localFileColMap             = fileColMap->getLocalMap();
  const LocalOrdinal invalidLocalOrdinal = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
  unsigned long long invalidCount        = 0;
  Kokkos::parallel_reduce(
      "Xpetra::IO::legacyReindexColumns",
      Kokkos::RangePolicy<execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i, unsigned long long& lclInvalidCount) {
        const LocalOrdinal lclCol = localFileColMap.getLocalElement(globalColumnsDevice(i));
        localColumnsDevice(i)     = lclCol;
        if (lclCol == invalidLocalOrdinal) {
          ++lclInvalidCount;
        }
      },
      invalidCount);
  TEUCHOS_TEST_FOR_EXCEPTION(invalidCount != 0,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: Column map construction missed " << invalidCount
                                                                           << " column GID(s) while reading legacy binary file '"
                                                                           << oldFileName << "'.");

  values_type valuesDevice("Xpetra::IO::legacyValuesDevice", localNnz);
  auto valuesDeviceHost = Kokkos::create_mirror_view(valuesDevice);
  Kokkos::parallel_for(
      "Xpetra::IO::legacyCopyValuesToDeviceHost",
      Kokkos::RangePolicy<host_execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i) {
        valuesDeviceHost(i) = static_cast<impl_scalar_type>(valuesHost(i));
      });
  Kokkos::deep_copy(valuesDevice, valuesDeviceHost);

  Tpetra::Import_Util::sortCrsEntries(rowPtrDevice, localColumnsDevice, valuesDevice);
  auto fileMatrix = Teuchos::rcp(new matrix_type(fileRowMap, fileColMap, rowPtrDevice, localColumnsDevice, valuesDevice));
  fileMatrix->fillComplete(fileDomainMap, fileRangeMap);

  if (rowMapInput.is_null()) {
    return fileMatrix;
  }

  RCP<const map_type> requestedDomainMap = domainMapInput.is_null() ? rowMapInput : domainMapInput;
  RCP<const map_type> requestedRangeMap  = rangeMapInput.is_null() ? rowMapInput : rangeMapInput;
  import_type importer(fileRowMap, rowMapInput);
  RCP<matrix_type> matrix;
  if (colMapInput.is_null()) {
    matrix = Teuchos::rcp(new matrix_type(rowMapInput,
                                          static_cast<size_t>(fileMatrix->getGlobalMaxNumRowEntries())));
  } else {
    matrix = Teuchos::rcp(new matrix_type(rowMapInput,
                                          colMapInput,
                                          static_cast<size_t>(fileMatrix->getGlobalMaxNumRowEntries())));
  }
  matrix->doImport(*fileMatrix, importer, Tpetra::INSERT);
  if (callFillComplete) {
    matrix->fillComplete(requestedDomainMap, requestedRangeMap);
  }
  return matrix;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readLegacyBinarySparseFile(const std::string&,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const RCP<const Teuchos::Comm<int>>&,
                           const bool) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: legacy binary matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value, void>::type
convertLegacyBinaryToBinary(const std::string& oldFileName,
                            const std::string& newFileName,
                            const RCP<const Teuchos::Comm<int>>& comm) {
  using binary_io_type = Tpetra::BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  auto matrix          = readLegacyBinarySparseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(oldFileName,
                                                                                      Teuchos::null,
                                                                                      Teuchos::null,
                                                                                      Teuchos::null,
                                                                                      Teuchos::null,
                                                                                      comm,
                                                                                      true);
  binary_io_type::writeSparseFile(newFileName, *matrix);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const RCP<const Teuchos::Comm<int>>& comm) {
  using binary_reader_type = Tpetra::BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  if (hasTpetraBinaryHeader(fileName, comm)) {
    return binary_reader_type::readSparseFile(fileName, comm);
  }

  return readLegacyBinarySparseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fileName,
                                                                               Teuchos::null,
                                                                               Teuchos::null,
                                                                               Teuchos::null,
                                                                               Teuchos::null,
                                                                               comm,
                                                                               true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
                     const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
                     const bool callFillComplete) {
  using binary_reader_type = Tpetra::BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  if (hasTpetraBinaryHeader(fileName, rowMap->getComm())) {
    return binary_reader_type::readSparseFile(fileName, rowMap, colMap, domainMap, rangeMap, callFillComplete);
  }

  return readLegacyBinarySparseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fileName,
                                                                               rowMap,
                                                                               colMap,
                                                                               domainMap,
                                                                               rangeMap,
                                                                               rowMap->getComm(),
                                                                               callFillComplete);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value, void>::type
convertLegacyBinaryToBinary(const std::string&, const std::string&, const RCP<const Teuchos::Comm<int>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: legacy binary conversion is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFile(const std::string& fileName,
                    const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using binary_reader_type = Tpetra::BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  return binary_reader_type::readDenseFile(fileName, map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFile(const std::string&,
                    const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: binary Tpetra multivector input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForLocalOrdinalScalar<LocalOrdinal>::value,
                        RCP<Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFileLocalOrdinal(const std::string& fileName,
                                const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using binary_reader_type = Tpetra::BinaryIO<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  return binary_reader_type::readDenseFile(fileName, map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForLocalOrdinalScalar<LocalOrdinal>::value,
                        RCP<Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFileLocalOrdinal(const std::string&,
                                const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             Exceptions::RuntimeError,
                             "Xpetra::IO: binary Tpetra local-ordinal multivector input is not available for this build. "
                                 << "That path requires a matching Tpetra::BinaryIO<LocalOrdinal,...> explicit instantiation.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

}  // namespace Details

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map2TpetraMap(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) {
  const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>& tmp_TMap = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(rcpFromRef(map));
  if (tmp_TMap == Teuchos::null)
    throw Exceptions::BadCast("Utils::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
  return tmp_TMap->getTpetra_Map();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& M) {
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tmp_Map = rcpFromRef(M);

  const RCP<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>& tmp_TMap =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node>>(tmp_Map);
  if (tmp_TMap != Teuchos::null) {
    RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> TMap = tmp_TMap->getTpetra_Map();
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeMapFile(fileName, *TMap);
    return;
  }

  throw Exceptions::BadCast("Could not cast to EpetraMap or TpetraMap in map writing");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  std::string mapfile = "map_" + fileName;
  Write(mapfile, *(vec.getMap()));

  RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);

  const RCP<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
  if (tmp_TVec != Teuchos::null) {
    RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
    return;
  }

  throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteLOMV(const std::string& fileName, const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  std::string mapfile = "map_" + fileName;
  Write(mapfile, *(vec.getMap()));

  RCP<const Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);
  const RCP<const Xpetra::TpetraMultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
  if (tmp_TVec != Teuchos::null) {
    RCP<const Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
    return;
  } else {
    throw Exceptions::RuntimeError("Xpetra cannot write MV<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> when the underlying library is Epetra.");
  }

  throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteGOMV(const std::string& fileName, const Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>& vec) {
  std::string mapfile = "map_" + fileName;
  Write(mapfile, *(vec.getMap()));

  RCP<const Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> tmp_Vec = Teuchos::rcpFromRef(vec);
  const RCP<const Xpetra::TpetraMultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TVec =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraMultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>(tmp_Vec);
  if (tmp_TVec != Teuchos::null) {
    RCP<const Tpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> TVec = tmp_TVec->getTpetra_MultiVector();
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>::writeDenseFile(fileName, TVec);
    return;
  } else {
    throw Exceptions::RuntimeError("Xpetra cannot write MV<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> when the underlying library is Epetra.");
  }

  throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in multivector writing");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Write(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps) {
  Write("rowmap_" + fileName, *(Op.getRowMap()));
  if (!Op.getDomainMap()->isSameAs(*(Op.getRowMap())) || writeAllMaps)
    Write("domainmap_" + fileName, *(Op.getDomainMap()));
  if (!Op.getRangeMap()->isSameAs(*(Op.getRowMap())) || writeAllMaps)
    Write("rangemap_" + fileName, *(Op.getRangeMap()));
  if (!Op.getColMap()->isSameAs(*(Op.getDomainMap())) || writeAllMaps)
    Write("colmap_" + fileName, *(Op.getColMap()));

  const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp =
      dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_CrsMtx = crsOp.getCrsMatrix();

  const RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_TCrsMtx =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_CrsMtx);
  if (tmp_TCrsMtx != Teuchos::null) {
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A = tmp_TCrsMtx->getTpetra_CrsMatrix();
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>::writeSparseFile(fileName, A);
    return;
  }
  const RCP<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& tmp_BlockCrs =
      Teuchos::rcp_dynamic_cast<const Xpetra::TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmp_CrsMtx);
  if (tmp_BlockCrs != Teuchos::null) {
    std::ofstream outstream(fileName, std::ofstream::out);
    Teuchos::FancyOStream ofs(Teuchos::rcpFromRef(outstream));
    tmp_BlockCrs->getTpetra_BlockCrsMatrix()->describe(ofs, Teuchos::VERB_EXTREME);
    return;
  }

  throw Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteLocal(const std::string& fileName, const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op) {
  const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>& crsOp =
      dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(Op);
  RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmp_CrsMtx = crsOp.getCrsMatrix();

  ArrayRCP<const size_t> rowptr_RCP;
  ArrayRCP<LocalOrdinal> rowptr2_RCP;
  ArrayRCP<const LocalOrdinal> colind_RCP;
  ArrayRCP<const Scalar> vals_RCP;
  tmp_CrsMtx->getAllValues(rowptr_RCP, colind_RCP, vals_RCP);

  ArrayView<const size_t> rowptr       = rowptr_RCP();
  ArrayView<const LocalOrdinal> colind = colind_RCP();
  ArrayView<const Scalar> vals         = vals_RCP();

  rowptr2_RCP.resize(rowptr.size());
  ArrayView<LocalOrdinal> rowptr2 = rowptr2_RCP();
  for (LocalOrdinal j = 0; j < rowptr.size(); j++)
    rowptr2[j] = rowptr[j];

  Teuchos::MatrixMarket::Raw::Writer<Scalar, LocalOrdinal> writer;
  writer.writeFile(fileName + "." + std::to_string(Op.getRowMap()->getComm()->getSize()) + "." + std::to_string(Op.getRowMap()->getComm()->getRank()),
                   rowptr2, colind, vals,
                   rowptr.size() - 1, Op.getColMap()->getLocalNumElements());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::WriteBlockedCrsMatrix(const std::string& fileName, const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Op, const bool& writeAllMaps) {
  // Write all matrix blocks with their maps
  for (size_t row = 0; row < Op.Rows(); ++row) {
    for (size_t col = 0; col < Op.Cols(); ++col) {
      RCP<const Matrix> m = Op.getMatrix(row, col);
      if (m != Teuchos::null) {  // skip empty blocks
        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(m) == Teuchos::null, Xpetra::Exceptions::BadCast,
                                   "Sub block matrix (" << row << "," << col << ") is not of type CrsMatrixWrap.");
        Write(fileName + toString(row) + toString(col) + ".m", *m, writeAllMaps);
      }
    }
  }

  // write map information of map extractors
  RCP<const MapExtractor> rangeMapExtractor  = Op.getRangeMapExtractor();
  RCP<const MapExtractor> domainMapExtractor = Op.getDomainMapExtractor();

  for (size_t row = 0; row < rangeMapExtractor->NumMaps(); ++row) {
    RCP<const Map> map = rangeMapExtractor->getMap(row);
    Write("subRangeMap_" + fileName + toString(row) + ".m", *map);
  }
  Write("fullRangeMap_" + fileName + ".m", *(rangeMapExtractor->getFullMap()));

  for (size_t col = 0; col < domainMapExtractor->NumMaps(); ++col) {
    RCP<const Map> map = domainMapExtractor->getMap(col);
    Write("subDomainMap_" + fileName + toString(col) + ".m", *map);
  }
  Write("fullDomainMap_" + fileName + ".m", *(domainMapExtractor->getFullMap()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ConvertLegacyBinaryToBinary(const std::string& oldFileName,
                                                                                const std::string& newFileName,
                                                                                Xpetra::UnderlyingLib lib,
                                                                                const RCP<const Teuchos::Comm<int>>& comm) {
  if (lib == Xpetra::UseTpetra) {
    Details::convertLegacyBinaryToBinary<Scalar, LocalOrdinal, GlobalOrdinal, Node>(oldFileName, newFileName, comm);
    return;
  }

  throw Exceptions::RuntimeError("Utils::ConvertLegacyBinaryToBinary : binary mode is only implemented for Xpetra::UseTpetra.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm, bool binary) {
  if (binary == false) {
    // Matrix Market file format (ASCII)
    if (lib == Xpetra::UseTpetra) {
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;

      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

      bool callFillComplete = true;

      RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, comm, callFillComplete);

      if (tA.is_null())
        throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

      return A;
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }
  } else {
    if (lib == Xpetra::UseTpetra) {
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;

      RCP<sparse_matrix_type> tA = Details::readBinarySparseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fileName, comm);

      if (tA.is_null())
        throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from BinaryIO::readSparseFile() is null.");

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

      return A;
    } else {
      throw Exceptions::RuntimeError("Utils::Read : binary mode is only implemented for Xpetra::UseTpetra.");
    }
  }  // if (binary == false) ... else

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read(const std::string& filename,
                                                    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                                                    RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,
                                                    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap,
                                                    const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap,
                                                    const bool callFillComplete,
                                                    const bool binary,
                                                    const bool tolerant,
                                                    const bool debug) {
  TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::Read() : rowMap cannot be null");

  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domain = (domainMap.is_null() ? rowMap : domainMap);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> range  = (rangeMap.is_null() ? rowMap : rangeMap);

  const Xpetra::UnderlyingLib lib = rowMap->lib();
  if (binary == false) {
    if (lib == Xpetra::UseTpetra) {
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

      const RCP<const map_type> tpetraRowMap    = Map2TpetraMap(*rowMap);
      RCP<const map_type> tpetraColMap          = (colMap.is_null() ? Teuchos::null : Map2TpetraMap(*colMap));
      const RCP<const map_type> tpetraRangeMap  = (rangeMap.is_null() ? tpetraRowMap : Map2TpetraMap(*rangeMap));
      const RCP<const map_type> tpetraDomainMap = (domainMap.is_null() ? tpetraRowMap : Map2TpetraMap(*domainMap));

      RCP<sparse_matrix_type> tA = reader_type::readSparseFile(filename, tpetraRowMap, tpetraColMap, tpetraDomainMap, tpetraRangeMap,
                                                               callFillComplete, tolerant, debug);
      if (tA.is_null())
        throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

      return A;
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }
  } else {
    if (lib == Xpetra::UseTpetra) {
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

      const RCP<const map_type> tpetraRowMap    = Map2TpetraMap(*rowMap);
      const RCP<const map_type> tpetraColMap    = (colMap.is_null() ? Teuchos::null : Map2TpetraMap(*colMap));
      const RCP<const map_type> tpetraDomainMap = (domainMap.is_null() ? tpetraRowMap : Map2TpetraMap(*domainMap));
      const RCP<const map_type> tpetraRangeMap  = (rangeMap.is_null() ? tpetraRowMap : Map2TpetraMap(*rangeMap));

      RCP<sparse_matrix_type> tA = Details::readBinarySparseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(filename, tpetraRowMap, tpetraColMap, tpetraDomainMap, tpetraRangeMap, callFillComplete);
      if (tA.is_null())
        throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from BinaryIO::readSparseFile() is null.");

      RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA1 = rcp(new Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tA));
      RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tmpA2       = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tmpA1);
      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A              = rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tmpA2));

      return A;
    } else {
      throw Exceptions::RuntimeError("Utils::Read : binary mode is only implemented for Xpetra::UseTpetra.");
    }
  }

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadLocal(const std::string& filename,
                                                                                                                                 const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap,
                                                                                                                                 RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap,
                                                                                                                                 const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domainMap,
                                                                                                                                 const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rangeMap,
                                                                                                                                 const bool callFillComplete,
                                                                                                                                 const bool binary,
                                                                                                                                 const bool tolerant,
                                                                                                                                 const bool debug) {
  TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::ReadLocal() : rowMap cannot be null");
  TEUCHOS_TEST_FOR_EXCEPTION(colMap.is_null(), Exceptions::RuntimeError, "Utils::ReadLocal() : colMap cannot be null");

  using matrix_type   = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_wrap_type = Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using crs_type      = Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> domain = (domainMap.is_null() ? rowMap : domainMap);
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> range  = (rangeMap.is_null() ? rowMap : rangeMap);

  std::string rankFilename = filename + "." + std::to_string(rowMap->getComm()->getSize()) + "." + std::to_string(rowMap->getComm()->getRank());
  RCP<matrix_type> A       = rcp(new crs_wrap_type(rowMap, colMap, 0));

  if (binary == false) {
    RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
    params->set("Parse tolerantly", tolerant);
    params->set("Debug mode", debug);

    LocalOrdinal numRows = rowMap->getLocalNumElements();
    LocalOrdinal numCols = colMap->getLocalNumElements();

    ArrayRCP<LocalOrdinal> rowptr2_RCP;
    ArrayRCP<LocalOrdinal> colind2_RCP;
    ArrayRCP<Scalar> vals2_RCP;

    Teuchos::MatrixMarket::Raw::Reader<Scalar, LocalOrdinal> reader;
    reader.readFile(rowptr2_RCP, colind2_RCP, vals2_RCP,
                    numRows, numCols,
                    rankFilename);

    RCP<crs_type> ACrs = Teuchos::rcp_dynamic_cast<crs_wrap_type>(A)->getCrsMatrix();

    ArrayRCP<size_t> rowptr_RCP;
    ArrayRCP<LocalOrdinal> colind_RCP;
    ArrayRCP<Scalar> vals_RCP;
    ACrs->allocateAllValues(colind2_RCP.size(), rowptr_RCP, colind_RCP, vals_RCP);

    rowptr_RCP.assign(rowptr2_RCP.begin(), rowptr2_RCP.end());
    colind_RCP = colind2_RCP;
    vals_RCP   = vals2_RCP;

    ACrs->setAllValues(rowptr_RCP, colind_RCP, vals_RCP);
  } else {
    // Custom file format (binary)
    std::ifstream ifs = std::ifstream(rankFilename.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << filename << "\"");

    int m, n, nnz;
    ifs.read(reinterpret_cast<char*>(&m), sizeof(m));
    ifs.read(reinterpret_cast<char*>(&n), sizeof(n));
    ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

    TEUCHOS_ASSERT_EQUALITY(Teuchos::as<int>(rowMap->getLocalNumElements()), m);

    Teuchos::ArrayRCP<size_t> rowptrRCP;
    Teuchos::ArrayRCP<LocalOrdinal> indicesRCP;
    Teuchos::ArrayRCP<Scalar> valuesRCP;

    RCP<crs_type> ACrs = Teuchos::rcp_dynamic_cast<crs_wrap_type>(A)->getCrsMatrix();

    ACrs->allocateAllValues(nnz, rowptrRCP, indicesRCP, valuesRCP);

    Teuchos::ArrayView<size_t> rowptr        = rowptrRCP();
    Teuchos::ArrayView<LocalOrdinal> indices = indicesRCP();
    Teuchos::ArrayView<Scalar> values        = valuesRCP();

    bool sorted = true;

    // Read in rowptr
    for (int i = 0; i < m; i++) {
      int row, rownnz;
      ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
      ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));

      rowptr[row + 1] += rownnz;
      ifs.seekg(sizeof(int) * rownnz + sizeof(double) * rownnz, ifs.cur);
    }
    for (int i = 0; i < m; i++)
      rowptr[i + 1] += rowptr[i];
    TEUCHOS_ASSERT(Teuchos::as<int>(rowptr[m]) == nnz);

    // reset to where the data starts
    ifs.seekg(sizeof(int) * 3, ifs.beg);

    // read in entries
    for (int i = 0; i < m; i++) {
      int row, rownnz;
      ifs.read(reinterpret_cast<char*>(&row), sizeof(row));
      ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
      size_t ptr = rowptr[row];
      for (int j = 0; j < rownnz; j++) {
        int index;
        ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
        indices[ptr] = Teuchos::as<LocalOrdinal>(index);
        if (j > 0)
          sorted = sorted & (indices[ptr - 1] < indices[ptr]);
        ++ptr;
      }
      ptr = rowptr[row];
      for (int j = 0; j < rownnz; j++) {
        double value;
        ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
        values[ptr] = Teuchos::as<Scalar>(value);
        ++ptr;
      }
      rowptr[row] += rownnz;
    }
    for (int i = m; i > 0; i--)
      rowptr[i] = rowptr[i - 1];
    rowptr[0] = 0;

    if (!sorted) {
      for (LocalOrdinal lclRow = 0; lclRow < m; lclRow++) {
        size_t rowBegin = rowptr[lclRow];
        size_t rowEnd   = rowptr[lclRow + 1];
        Tpetra::sort2(&indices[rowBegin], &indices[rowEnd], &values[rowBegin]);
      }
    }

    ACrs->setAllValues(rowptrRCP, indicesRCP, valuesRCP);
  }

  if (callFillComplete)
    A->fillComplete(domainMap, rangeMap);
  return A;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(const std::string& fileName,
                                                                                                                                   const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                                                                                                                                   const bool binary) {
  Xpetra::UnderlyingLib lib = map->lib();

  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> mm_reader_type;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> multivector_type;

    RCP<const map_type> temp = toTpetra(map);
    RCP<multivector_type> TMV;
    if (binary) {
      TMV = Details::readBinaryDenseFile<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fileName, temp);
    } else {
      TMV = mm_reader_type::readDenseFile(fileName, map->getComm(), temp, false, false, false);
    }
    RCP<MultiVector> rmv = Xpetra::toXpetra(TMV);
    return rmv;
  } else {
    throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
  }

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVectorLO(const std::string& fileName,
                                                                                                                                           const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                                                                                                                                           const bool binary) {
  Xpetra::UnderlyingLib lib = map->lib();

  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrix<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> mm_reader_type;
    typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;
    typedef Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> multivector_type;

    RCP<const map_type> temp = toTpetra(map);
    RCP<multivector_type> TMV;
    if (binary) {
      TMV = Details::readBinaryDenseFileLocalOrdinal<LocalOrdinal, GlobalOrdinal, Node>(fileName, temp);
    } else {
      TMV = mm_reader_type::readDenseFile(fileName, map->getComm(), temp, false, false, false);
    }
    RCP<Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>> rmv = Xpetra::toXpetra(TMV);
    return rmv;
  } else {
    throw Exceptions::RuntimeError("Utils::ReadMultiVectorLO : only implemented for Tpetra");
  }

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMap(const std::string& fileName,
                                                                                                                 Xpetra::UnderlyingLib lib,
                                                                                                                 const RCP<const Teuchos::Comm<int>>& comm,
                                                                                                                 const bool binary) {
  if (lib == Xpetra::UseTpetra) {
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
    typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> mm_reader_type;

    RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> tMap;
    if (binary) {
      tMap = Details::readBinaryMap<LocalOrdinal, GlobalOrdinal, Node>(fileName, comm);
    } else {
      tMap = mm_reader_type::readMapFile(fileName, comm, false, false, false);
    }
    if (tMap.is_null())
      throw Exceptions::RuntimeError("The Tpetra::Map returned from readMapFile() is null.");

    return Xpetra::toXpetra(tMap);
  } else {
    throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
  }

  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadBlockedCrsMatrix(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int>>& comm) {
  size_t numBlocks = 2;  // TODO user parameter?

  std::vector<RCP<const Map>> rangeMapVec;
  for (size_t row = 0; row < numBlocks; ++row) {
    RCP<const Map> map = ReadMap("subRangeMap_" + fileName + toString(row) + ".m", lib, comm);
    rangeMapVec.push_back(map);
  }
  RCP<const Map> fullRangeMap = ReadMap("fullRangeMap_" + fileName + ".m", lib, comm);

  std::vector<RCP<const Map>> domainMapVec;
  for (size_t col = 0; col < numBlocks; ++col) {
    RCP<const Map> map = ReadMap("subDomainMap_" + fileName + toString(col) + ".m", lib, comm);
    domainMapVec.push_back(map);
  }
  RCP<const Map> fullDomainMap = ReadMap("fullDomainMap_" + fileName + ".m", lib, comm);

  /*std::vector<RCP<const XpMap> > testRgMapVec;
  for(size_t r = 0; r < numBlocks; ++r) {
    RCP<const XpMap> map = ReadMap("rangemap_" + fileName + toString<size_t>(r) + "0.m", lib, comm);
    testRgMapVec.push_back(map);
  }
  std::vector<RCP<const XpMap> > testDoMapVec;
  for(size_t c = 0; c < numBlocks; ++c) {
    RCP<const XpMap> map = ReadMap("domainmap_" + fileName + "0" + toString<size_t>(c) + ".m", lib, comm);
    testDoMapVec.push_back(map);
  }*/

  // create map extractors

  // range map extractor
  bool bRangeUseThyraStyleNumbering = false;
  /*GlobalOrdinal gMinGids = 0;
  for(size_t v = 0; v < testRgMapVec.size(); ++v) {
    gMinGids += testRgMapVec[v]->getMinAllGlobalIndex();
  }
  if ( gMinGids==0 && testRgMapVec.size() > 1 ) bRangeUseThyraStyleNumbering = true;
  */
  RCP<const MapExtractor> rangeMapExtractor =
      Teuchos::rcp(new MapExtractor(fullRangeMap, rangeMapVec, bRangeUseThyraStyleNumbering));

  // domain map extractor
  bool bDomainUseThyraStyleNumbering = false;
  /*gMinGids = 0;
  for(size_t v = 0; v < testDoMapVec.size(); ++v) {
    gMinGids += testDoMapVec[v]->getMinAllGlobalIndex();
  }
  if ( gMinGids==0 && testDoMapVec.size() > 1 ) bDomainUseThyraStyleNumbering = true;
  */
  RCP<const MapExtractor> domainMapExtractor =
      Teuchos::rcp(new MapExtractor(fullDomainMap, domainMapVec, bDomainUseThyraStyleNumbering));

  RCP<BlockedCrsMatrix> bOp = Teuchos::rcp(new BlockedCrsMatrix(rangeMapExtractor, domainMapExtractor, 33));

  // Read all sub-matrices with their maps and place into blocked operator
  for (size_t row = 0; row < numBlocks; ++row) {
    for (size_t col = 0; col < numBlocks; ++col) {
      RCP<const Map> rowSubMap = ReadMap("rowmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
      RCP<const Map> colSubMap = ReadMap("colmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
      RCP<const Map> domSubMap = ReadMap("domainmap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
      RCP<const Map> ranSubMap = ReadMap("rangemap_" + fileName + toString(row) + toString(col) + ".m", lib, comm);
      RCP<Matrix> mat          = Read(fileName + toString(row) + toString(col) + ".m", rowSubMap, colSubMap, domSubMap, ranSubMap);
      bOp->setMatrix(row, col, mat);
    }
  }

  bOp->fillComplete();

  return bOp;
}  // ReadBlockedCrsMatrix

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
std::string IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toString(const T& what) {
  std::ostringstream buf;
  buf << what;
  return buf.str();
}
}  // namespace Xpetra
