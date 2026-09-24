// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BINARYIO_LEGACY_HPP
#define TPETRA_BINARYIO_LEGACY_HPP

/// \file Tpetra_BinaryIO_Legacy.hpp
/// \brief Format-detecting binary matrix input, legacy binary format support,
///   and scalar-availability traits for Tpetra's binary I/O.
///
/// This header holds functionality that used to live in Xpetra::IO but is
/// really Tpetra-level.  It is header-only and SFINAE-gated so that it compiles
/// for every Scalar type (including Sacado/Stokhos ensemble scalars), even
/// though the underlying \c Tpetra::BinaryIO reader is only explicitly
/// instantiated for native Tpetra scalar types.
///
/// The entry points (all in \c Tpetra::Details) are:
///   - \c binaryIOAvailableForScalar / \c binaryIOAvailableForLocalOrdinalScalar
///     — traits describing which Scalar types have a usable BinaryIO reader.
///   - \c hasTpetraBinaryHeader — sniff whether a file is in the current Tpetra
///     binary format or the legacy Xpetra binary matrix format.
///   - \c readBinarySparseFile — read a sparse matrix, auto-detecting the format
///     and falling back to the legacy reader when needed.
///   - \c readLegacyBinarySparseFile — read the legacy Xpetra binary matrix
///     format directly.
///   - \c convertLegacyBinaryToBinary — convert a legacy file to the current
///     format.
///   - \c readBinaryDenseFile / \c readBinaryDenseFileLocalOrdinal — read a
///     dense (multi)vector in the current Tpetra binary format.

#include "TpetraCore_config.h"

#include "Tpetra_BinaryIO.hpp"
#include "Tpetra_ReadBinaryMapFile.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include "Kokkos_Core.hpp"

#include <complex>
#include <cstring>
#include <fstream>
#include <limits>
#include <string>
#include <type_traits>

namespace Tpetra {
namespace Details {

/// \brief Whether \c Tpetra::BinaryIO is explicitly instantiated for \c Scalar.
///
/// Tpetra::BinaryIO is only ETI'd for native Tpetra scalar types, not for
/// Sacado/Stokhos ensemble scalars.  Reading binary files for other scalar
/// types therefore throws rather than failing to link.
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

/// \brief Whether \c Tpetra::BinaryIO is instantiated with Scalar = LocalOrdinal.
///
/// This gates the "read a MultiVector whose Scalar is LocalOrdinal" path.
template <class LocalOrdinal>
struct binaryIOAvailableForLocalOrdinalScalar : std::false_type {};

#if !defined(HAVE_TPETRA_REDUCED_ETI) && !defined(HAVE_TPETRA_INST_INT_INT)
template <>
struct binaryIOAvailableForLocalOrdinalScalar<int> : std::true_type {};
#endif

/// \brief Return true if the file begins with the current Tpetra binary I/O
///   header magic, false if it looks like a legacy Xpetra binary matrix file.
inline bool hasTpetraBinaryHeader(const std::string& fileName,
                                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
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
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to open binary file '" << fileName << "'.");
  return hasHeader != 0;
}

/// \brief On-disk header of a legacy Xpetra binary matrix file.
struct LegacyBinaryHeader {
  int numRows;
  int numCols;
  int numEntries;
};

inline unsigned long long legacyBinaryFileSize(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::binary | std::ios::ate);
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to open legacy binary file '" << fileName << "'.");
  const std::streamoff size = in.tellg();
  TEUCHOS_TEST_FOR_EXCEPTION(size < 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to determine the size of legacy binary file '" << fileName << "'.");
  return static_cast<unsigned long long>(size);
}

inline LegacyBinaryHeader readAndValidateLegacyBinaryHeader(const std::string& fileName) {
  std::ifstream in(fileName.c_str(), std::ios::binary);
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to open legacy binary file '" << fileName << "'.");

  LegacyBinaryHeader header = {0, 0, 0};
  in.read(reinterpret_cast<char*>(&header.numRows), sizeof(header.numRows));
  in.read(reinterpret_cast<char*>(&header.numCols), sizeof(header.numCols));
  in.read(reinterpret_cast<char*>(&header.numEntries), sizeof(header.numEntries));
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File '" << fileName << "' is not a valid legacy Xpetra binary matrix file.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.numRows < 0 || header.numCols < 0 || header.numEntries < 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Legacy binary file '" << fileName << "' has negative dimensions or entry count.");

  const unsigned long long fileSize = legacyBinaryFileSize(fileName);
  const unsigned long long minSize  = static_cast<unsigned long long>(3 * sizeof(int)) +
                                     static_cast<unsigned long long>(header.numRows) * static_cast<unsigned long long>(2 * sizeof(int)) +
                                     static_cast<unsigned long long>(header.numEntries) * static_cast<unsigned long long>(sizeof(int) + sizeof(double));
  TEUCHOS_TEST_FOR_EXCEPTION(fileSize != minSize,
                             std::runtime_error,
                             "Tpetra::BinaryIO: File '" << fileName << "' is not a valid legacy Xpetra binary matrix file."
                                                        << " Expected " << minSize << " bytes from its header, but file has " << fileSize << " bytes.");
  return header;
}

template <class T>
void readLegacyBinaryValue(std::ifstream& in, T& value, const std::string& fileName, const char label[]) {
  in.read(reinterpret_cast<char*>(&value), sizeof(T));
  TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to read " << label << " from legacy binary file '" << fileName << "'.");
}

/// \brief Read a matrix stored in the legacy Xpetra binary matrix format.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readLegacyBinarySparseFile(const std::string& oldFileName,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMapInput,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& colMapInput,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMapInput,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMapInput,
                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
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

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> rowLengths("Tpetra::BinaryIO::legacyRowLengths", fileLocalNumRows);
  Kokkos::deep_copy(rowLengths, 0ull);
  unsigned long long entriesRead = 0;

  std::ifstream in;
  if (myRank == 0) {
    in.open(oldFileName.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open legacy binary file '" << oldFileName << "'.");
    in.seekg(static_cast<std::streamoff>(3 * sizeof(int)), std::ios::beg);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Failed to seek past the header of legacy binary file '" << oldFileName << "'.");

    for (int record = 0; record < m; ++record) {
      int row    = 0;
      int rownnz = 0;
      readLegacyBinaryValue(in, row, oldFileName, "row index");
      readLegacyBinaryValue(in, rownnz, oldFileName, "row entry count");
      TEUCHOS_TEST_FOR_EXCEPTION(row < 0 || row >= m || rownnz < 0,
                                 std::runtime_error,
                                 "Tpetra::BinaryIO: Legacy binary file '" << oldFileName << "' has an invalid row record.");
      entriesRead += static_cast<unsigned long long>(rownnz);
      rowLengths(static_cast<size_t>(row)) += static_cast<unsigned long long>(rownnz);
      for (int j = 0; j < rownnz; ++j) {
        int col = 0;
        readLegacyBinaryValue(in, col, oldFileName, "column index");
        TEUCHOS_TEST_FOR_EXCEPTION(col < 0 || col >= n,
                                   std::runtime_error,
                                   "Tpetra::BinaryIO: Legacy binary file '" << oldFileName << "' has a column index outside [0, n).");
      }
      for (int j = 0; j < rownnz; ++j) {
        double value = 0.0;
        readLegacyBinaryValue(in, value, oldFileName, "matrix value");
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(entriesRead != static_cast<unsigned long long>(nnz),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Legacy binary file '" << oldFileName << "' header entry count does not match row records.");
    TEUCHOS_TEST_FOR_EXCEPTION(entriesRead > static_cast<unsigned long long>(std::numeric_limits<size_t>::max()),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Legacy binary file '" << oldFileName << "' has too many entries for this platform.");
  }

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> rowPtrHost("Tpetra::BinaryIO::legacyRowPtr", fileLocalNumRows + 1);
  unsigned long long localNnz64 = 0;
  Kokkos::parallel_scan(
      "Tpetra::BinaryIO::legacyRowPtrScan",
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
  Kokkos::View<unsigned long long*, Kokkos::HostSpace> nextPtr("Tpetra::BinaryIO::legacyNextPtr", fileLocalNumRows);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::legacyInitNextPtr",
      Kokkos::RangePolicy<host_execution_space>(0, fileLocalNumRows),
      KOKKOS_LAMBDA(const size_t row) {
        nextPtr(row) = rowPtrHost(row);
      });

  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> globalColumnsHost("Tpetra::BinaryIO::legacyGlobalColumns", localNnz);
  Kokkos::View<Scalar*, Kokkos::HostSpace> valuesHost("Tpetra::BinaryIO::legacyValues", localNnz);

  if (myRank == 0) {
    in.clear();
    in.seekg(static_cast<std::streamoff>(3 * sizeof(int)), std::ios::beg);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Failed to rewind legacy binary file '" << oldFileName << "'.");

    for (int record = 0; record < m; ++record) {
      int row    = 0;
      int rownnz = 0;
      readLegacyBinaryValue(in, row, oldFileName, "row index");
      readLegacyBinaryValue(in, rownnz, oldFileName, "row entry count");

      Kokkos::View<int*, Kokkos::HostSpace> columns("Tpetra::BinaryIO::legacyRowColumns", static_cast<size_t>(rownnz));
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

  rowptr_type rowPtrDevice("Tpetra::BinaryIO::legacyRowPtrDevice", fileLocalNumRows + 1);
  auto rowPtrDeviceHost = Kokkos::create_mirror_view(rowPtrDevice);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::legacyCopyRowPtrToDeviceHost",
      Kokkos::RangePolicy<host_execution_space>(0, fileLocalNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        rowPtrDeviceHost(i) = static_cast<typename rowptr_type::non_const_value_type>(rowPtrHost(i));
      });
  Kokkos::deep_copy(rowPtrDevice, rowPtrDeviceHost);

  Kokkos::View<GlobalOrdinal*, device_type> globalColumnsDevice("Tpetra::BinaryIO::legacyGlobalColumnsDevice", localNnz);
  Kokkos::deep_copy(globalColumnsDevice, globalColumnsHost);

  colidx_type localColumnsDevice("Tpetra::BinaryIO::legacyLocalColumns", localNnz);
  const auto localFileColMap             = fileColMap->getLocalMap();
  const LocalOrdinal invalidLocalOrdinal = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
  unsigned long long invalidCount        = 0;
  Kokkos::parallel_reduce(
      "Tpetra::BinaryIO::legacyReindexColumns",
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
                             std::runtime_error,
                             "Tpetra::BinaryIO: Column map construction missed " << invalidCount
                                                                                 << " column GID(s) while reading legacy binary file '"
                                                                                 << oldFileName << "'.");

  values_type valuesDevice("Tpetra::BinaryIO::legacyValuesDevice", localNnz);
  auto valuesDeviceHost = Kokkos::create_mirror_view(valuesDevice);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::legacyCopyValuesToDeviceHost",
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

  Teuchos::RCP<const map_type> requestedDomainMap = domainMapInput.is_null() ? rowMapInput : domainMapInput;
  Teuchos::RCP<const map_type> requestedRangeMap  = rangeMapInput.is_null() ? rowMapInput : rangeMapInput;
  import_type importer(fileRowMap, rowMapInput);
  Teuchos::RCP<matrix_type> matrix;
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
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readLegacyBinarySparseFile(const std::string&,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                           const Teuchos::RCP<const Teuchos::Comm<int>>&,
                           const bool) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: legacy binary matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

/// \brief Convert a legacy Xpetra binary matrix file to the current Tpetra
///   binary format.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value, void>::type
convertLegacyBinaryToBinary(const std::string& oldFileName,
                            const std::string& newFileName,
                            const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
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
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value, void>::type
convertLegacyBinaryToBinary(const std::string&, const std::string&, const Teuchos::RCP<const Teuchos::Comm<int>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: legacy binary conversion is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
}

/// \brief Read a sparse matrix from a binary file, auto-detecting whether it is
///   in the current Tpetra binary format or the legacy Xpetra binary format.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
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
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string&,
                     const Teuchos::RCP<const Teuchos::Comm<int>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: binary Tpetra matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string& fileName,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& colMap,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rangeMap,
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
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinarySparseFile(const std::string&,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&,
                     const bool) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: binary Tpetra matrix input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

/// \brief Read a dense (multi)vector from a file in the current Tpetra binary
///   format.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFile(const std::string& fileName,
                    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using binary_reader_type = Tpetra::BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  return binary_reader_type::readDenseFile(fileName, map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForScalar<Scalar>::value,
                        Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFile(const std::string&,
                    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: binary Tpetra multivector input is not available for this Scalar type. "
                                 << "Tpetra::BinaryIO is only explicitly instantiated for native Tpetra scalar types, "
                                 << "not for Sacado/Stokhos ensemble scalars.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

/// \brief Read a dense (multi)vector whose Scalar is LocalOrdinal from a file
///   in the current Tpetra binary format.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<binaryIOAvailableForLocalOrdinalScalar<LocalOrdinal>::value,
                        Teuchos::RCP<Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFileLocalOrdinal(const std::string& fileName,
                                const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map) {
  using binary_reader_type = Tpetra::BinaryIO<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  return binary_reader_type::readDenseFile(fileName, map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename std::enable_if<!binaryIOAvailableForLocalOrdinalScalar<LocalOrdinal>::value,
                        Teuchos::RCP<Tpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>>>::type
readBinaryDenseFileLocalOrdinal(const std::string&,
                                const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>&) {
  TEUCHOS_TEST_FOR_EXCEPTION(true,
                             std::runtime_error,
                             "Tpetra::BinaryIO: binary Tpetra local-ordinal multivector input is not available for this build. "
                                 << "That path requires a matching Tpetra::BinaryIO<LocalOrdinal,...> explicit instantiation.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_BINARYIO_LEGACY_HPP
