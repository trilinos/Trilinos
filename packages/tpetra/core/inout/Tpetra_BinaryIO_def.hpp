// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BINARYIO_DEF_HPP
#define TPETRA_BINARYIO_DEF_HPP

#include "Tpetra_BinaryIO_decl.hpp"

#include "Tpetra_BinaryIO_Helpers.hpp"
#include "Tpetra_ReadBinaryMapFile_decl.hpp"
#include "Tpetra_Details_MpiTypeTraits.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#include "Tpetra_Details_makeColMap.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Import_Util2.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_Sort.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace Tpetra {

namespace Details {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
makeColumnMapFromGlobalColumns(const Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace>& globalColumns,
                               const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap) {
  using map_type        = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using memory_space    = typename Node::memory_space;
  const size_t localNnz = globalColumns.extent(0);

  int hasLocalEntries     = localNnz == 0 ? 0 : 1;
  int allHaveLocalEntries = 0;
  Teuchos::reduceAll(*domainMap->getComm(), Teuchos::REDUCE_MIN, 1, &hasLocalEntries, &allHaveLocalEntries);
  if (allHaveLocalEntries != 0) {
    Kokkos::View<GlobalOrdinal*, memory_space> globalColumnsDefaultMemory("Tpetra::BinaryIO::globalColumnsForColMap", localNnz);
    Kokkos::deep_copy(globalColumnsDefaultMemory, globalColumns);

    Teuchos::RCP<const map_type> colMap;
    std::ostringstream errStrm;
    const int err = Tpetra::Details::makeColMap<LocalOrdinal, GlobalOrdinal, Node>(colMap,
                                                                                   domainMap,
                                                                                   globalColumnsDefaultMemory,
                                                                                   &errStrm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0 || colMap.is_null(),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Failed to construct column map from sparse matrix file. " << errStrm.str());
    return colMap;
  }

  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> colGids("Tpetra::BinaryIO::colGids", localNnz);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::copyColumnGidsForMap",
      Kokkos::RangePolicy<host_execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i) {
        colGids(i) = globalColumns(i);
      });
  if (localNnz > 1) {
    Kokkos::sort(colGids);
  }

  size_t uniqueCount = 0;
  Kokkos::parallel_scan(
      "Tpetra::BinaryIO::uniqueColumnGidsForMap",
      Kokkos::RangePolicy<host_execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i, size_t& count, const bool finalPass) {
        const bool isUnique = i == 0 || colGids(i) != colGids(i - 1);
        if (isUnique) {
          if (finalPass) {
            colGids(count) = colGids(i);
          }
          ++count;
        }
      },
      uniqueCount);

  const LocalOrdinal localCount = binaryIOCheckedLocalOrdinalCount<LocalOrdinal>(uniqueCount, "column-map local entry count");
  using device_type             = typename map_type::device_type;
  Kokkos::View<GlobalOrdinal*, device_type> colGidsDevice("Tpetra::BinaryIO::colGidsDevice", static_cast<size_t>(localCount));
  Kokkos::deep_copy(colGidsDevice, Kokkos::subview(colGids, Kokkos::make_pair(size_t(0), static_cast<size_t>(localCount))));
  return Teuchos::rcp(new map_type(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                                   colGidsDevice,
                                   domainMap->getIndexBase(),
                                   domainMap->getComm()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
buildSparseMatrixFromLocalCrsViews(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& rowMap,
                                   const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& inputColMap,
                                   const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& domainMap,
                                   const Kokkos::View<unsigned long long*, Kokkos::HostSpace>& rowPtrHost,
                                   const Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace>& globalColumnsHost,
                                   const Kokkos::View<Scalar*, Kokkos::HostSpace>& valuesHost) {
  using matrix_type          = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using local_graph_type     = typename matrix_type::local_graph_device_type;
  using rowptr_type          = typename local_graph_type::row_map_type::non_const_type;
  using colidx_type          = typename local_graph_type::entries_type::non_const_type;
  using values_type          = typename matrix_type::local_matrix_device_type::values_type::non_const_type;
  using impl_scalar_type     = typename matrix_type::impl_scalar_type;
  using device_type          = typename matrix_type::device_type;
  using execution_space      = typename device_type::execution_space;
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;

  const size_t localNumRows = rowMap->getLocalNumElements();
  const size_t localNnz     = globalColumnsHost.extent(0);

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> colMap = inputColMap;
  if (colMap.is_null()) {
    colMap = makeColumnMapFromGlobalColumns<LocalOrdinal, GlobalOrdinal, Node>(globalColumnsHost, domainMap);
  }

  rowptr_type rowPtrDevice("Tpetra::BinaryIO::rowPtr", localNumRows + 1);
  auto rowPtrDeviceHost                  = Kokkos::create_mirror_view(rowPtrDevice);
  using rowptr_value_type                = typename rowptr_type::non_const_value_type;
  unsigned long long rowPtrOverflowCount = 0;
  Kokkos::parallel_reduce(
      "Tpetra::BinaryIO::checkRowPtrRange",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows + 1),
      KOKKOS_LAMBDA(const size_t i, unsigned long long& overflowCount) {
        if (rowPtrHost(i) > static_cast<unsigned long long>(std::numeric_limits<rowptr_value_type>::max())) {
          ++overflowCount;
        }
      },
      rowPtrOverflowCount);
  TEUCHOS_TEST_FOR_EXCEPTION(rowPtrOverflowCount != 0,
                             std::overflow_error,
                             "Tpetra::BinaryIO: local sparse row offset does not fit in the local row pointer type.");
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::copyRowPtrToDeviceHost",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        rowPtrDeviceHost(i) = static_cast<rowptr_value_type>(rowPtrHost(i));
      });
  Kokkos::deep_copy(rowPtrDevice, rowPtrDeviceHost);

  Kokkos::View<GlobalOrdinal*, device_type> globalColumnsDevice("Tpetra::BinaryIO::globalColumns", localNnz);
  Kokkos::deep_copy(globalColumnsDevice, globalColumnsHost);

  colidx_type localColumnsDevice("Tpetra::BinaryIO::localColumns", localNnz);
  const auto localColMap                       = colMap->getLocalMap();
  const LocalOrdinal invalidLocalOrdinal       = Tpetra::Details::OrdinalTraits<LocalOrdinal>::invalid();
  unsigned long long invalidLocalColumnIdCount = 0;
  Kokkos::parallel_reduce(
      "Tpetra::BinaryIO::reindexColumns",
      Kokkos::RangePolicy<execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i, unsigned long long& invalidCount) {
        const LocalOrdinal lclCol = localColMap.getLocalElement(globalColumnsDevice(i));
        localColumnsDevice(i)     = lclCol;
        if (lclCol == invalidLocalOrdinal) {
          ++invalidCount;
        }
      },
      invalidLocalColumnIdCount);
  TEUCHOS_TEST_FOR_EXCEPTION(invalidLocalColumnIdCount != 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Column map is missing " << invalidLocalColumnIdCount
                                                                        << " column GID(s) stored in the sparse matrix file.");

  values_type valuesDevice("Tpetra::BinaryIO::values", localNnz);
  auto valuesDeviceHost = Kokkos::create_mirror_view(valuesDevice);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::copyValuesToDeviceHost",
      Kokkos::RangePolicy<host_execution_space>(0, localNnz),
      KOKKOS_LAMBDA(const size_t i) {
        valuesDeviceHost(i) = static_cast<impl_scalar_type>(valuesHost(i));
      });
  Kokkos::deep_copy(valuesDevice, valuesDeviceHost);

  Tpetra::Import_Util::sortCrsEntries(rowPtrDevice, localColumnsDevice, valuesDevice);
  return Teuchos::rcp(new matrix_type(rowMap, colMap, rowPtrDevice, localColumnsDevice, valuesDevice));
}

}  // namespace Details

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FileHeader
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::makeBaseHeader(const unsigned long long objectKind) {
  FileHeader header;
  std::memset(&header, 0, sizeof(FileHeader));
  const char magic[8] = {'T', 'P', 'B', 'I', 'O', '0', '0', '1'};
  std::memcpy(header.magic, magic, sizeof(magic));
  header.version            = fileVersion;
  header.byteOrderMarker    = byteOrderMarker;
  header.objectKind         = objectKind;
  header.scalarSize         = sizeof(Scalar);
  header.scalarFlags        = typeFlags<Scalar>();
  header.localOrdinalSize   = sizeof(LocalOrdinal);
  header.localOrdinalFlags  = typeFlags<LocalOrdinal>();
  header.globalOrdinalSize  = sizeof(GlobalOrdinal);
  header.globalOrdinalFlags = typeFlags<GlobalOrdinal>();
  return header;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::validateHeader(const FileHeader& header,
                                                                         const unsigned long long expectedObjectKind) {
  const char expectedMagic[8] = {'T', 'P', 'B', 'I', 'O', '0', '0', '1'};
  TEUCHOS_TEST_FOR_EXCEPTION(std::memcmp(header.magic, expectedMagic, sizeof(expectedMagic)) != 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: File does not have a recognized Tpetra binary I/O header. "
                                 << "Legacy Xpetra binary matrix files must be converted with "
                                 << "Xpetra::IO::ConvertLegacyBinaryToBinary before reading with Tpetra::BinaryIO.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.version != fileVersion,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Unsupported file version " << header.version << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(header.byteOrderMarker != byteOrderMarker,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Unsupported byte order marker.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.objectKind != expectedObjectKind,
                             std::runtime_error,
                             "Tpetra::BinaryIO: File object kind " << header.objectKind
                                                                   << " does not match expected kind " << expectedObjectKind << ".");
  if (expectedObjectKind != MAP_OBJECT) {
    TEUCHOS_TEST_FOR_EXCEPTION(header.scalarSize != sizeof(Scalar) || header.scalarFlags != typeFlags<Scalar>(),
                               std::runtime_error,
                               "Tpetra::BinaryIO: File scalar type does not match this BinaryIO instantiation.");
  }
  TEUCHOS_TEST_FOR_EXCEPTION(header.localOrdinalSize != sizeof(LocalOrdinal) || header.localOrdinalFlags != typeFlags<LocalOrdinal>(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File local ordinal type does not match this BinaryIO instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.globalOrdinalSize != sizeof(GlobalOrdinal) || header.globalOrdinalFlags != typeFlags<GlobalOrdinal>(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File global ordinal type does not match this BinaryIO instantiation.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::typeFlags() {
  unsigned long long flags = 0;
  if (Teuchos::ScalarTraits<T>::isOrdinal) {
    flags |= TYPE_IS_ORDINAL;
  }
  if (Teuchos::ScalarTraits<T>::isComplex) {
    flags |= TYPE_IS_COMPLEX;
  }
  if (std::is_floating_point<T>::value || Teuchos::ScalarTraits<T>::isComplex) {
    flags |= TYPE_IS_FLOATING;
  }
  if (std::is_signed<T>::value || Teuchos::ScalarTraits<T>::isComplex || std::is_floating_point<T>::value) {
    flags |= TYPE_IS_SIGNED;
  }
  return flags;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapFlags(const map_type& map) {
  unsigned long long flags = 0;
  if (map.isDistributed()) {
    flags |= 1ull << 0;
  }
  if (map.isContiguous()) {
    flags |= 1ull << 1;
  }
  if (map.isOneToOne()) {
    flags |= 1ull << 2;
  }
  return flags;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkedByteCount(const unsigned long long count,
                                                                                         const size_t elementSize) {
  TEUCHOS_TEST_FOR_EXCEPTION(count > 0 &&
                                 elementSize > static_cast<size_t>(std::numeric_limits<unsigned long long>::max() / count),
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte count overflow while computing transfer size.");
  return count * static_cast<unsigned long long>(elementSize);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkedByteOffset(const unsigned long long dataOffset,
                                                                                          const unsigned long long globalOffset,
                                                                                          const size_t elementSize) {
  const unsigned long long payloadOffset = checkedByteCount(globalOffset, elementSize);
  TEUCHOS_TEST_FOR_EXCEPTION(dataOffset > std::numeric_limits<unsigned long long>::max() - payloadOffset,
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte offset overflow while computing file offset.");
  return dataOffset + payloadOffset;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::broadcastBytesFromRoot(char data[],
                                                                                 const unsigned long long byteCount,
                                                                                 const trcp_tcomm_t& comm) {
  unsigned long long remaining      = byteCount;
  char* current                     = data;
  const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<int>::max());
  while (remaining > 0) {
    const unsigned long long chunk = std::min(remaining, maxChunk);
    Teuchos::broadcast(*comm, 0, static_cast<int>(chunk), current);
    current += static_cast<std::ptrdiff_t>(chunk);
    remaining -= chunk;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::checkFileOpen(const bool success,
                                                                        const std::string& filename,
                                                                        const trcp_tcomm_t& comm,
                                                                        const char mode[]) {
  int opened = success ? 1 : 0;
  Teuchos::broadcast(*comm, 0, 1, &opened);
  TEUCHOS_TEST_FOR_EXCEPTION(opened == 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to open file '" << filename << "' for " << mode << ".");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeHeaderToNewFile(const std::string& filename,
                                                                               const FileHeader& header,
                                                                               const trcp_tcomm_t& comm) {
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()),
                                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                                  MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed for '" << filename << "'.");
    const int truncateErr = MPI_File_set_size(file, static_cast<MPI_Offset>(0));
    TEUCHOS_TEST_FOR_EXCEPTION(truncateErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_set_size failed while truncating '" << filename << "'.");
    if (comm->getRank() == 0) {
      MPI_Status status;
      const int writeErr = MPI_File_write_at(file, 0, const_cast<FileHeader*>(&header),
                                             sizeof(FileHeader), MPI_BYTE, &status);
      TEUCHOS_TEST_FOR_EXCEPTION(writeErr != MPI_SUCCESS, std::runtime_error,
                                 "Tpetra::BinaryIO: MPI_File_write_at failed while writing header.");
    }
    MPI_File_close(&file);
    comm->barrier();
    return;
  }
#endif
  bool success = true;
  if (comm->getRank() == 0) {
    std::ofstream out(filename.c_str(), std::ios::binary | std::ios::trunc);
    success = static_cast<bool>(out);
    if (success) {
      out.write(reinterpret_cast<const char*>(&header), sizeof(FileHeader));
      success = static_cast<bool>(out);
    }
  }
  checkFileOpen(success, filename, comm, "writing");
  comm->barrier();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FileHeader
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readHeaderFromFile(const std::string& filename,
                                                                        const trcp_tcomm_t& comm) {
  FileHeader header;
  std::memset(&header, 0, sizeof(FileHeader));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()),
                                  MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed for '" << filename << "'.");
    MPI_Status status;
    const int readErr = MPI_File_read_at_all(file, 0, &header, sizeof(FileHeader), MPI_BYTE, &status);
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(readErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_read_at_all failed while reading header.");
    return header;
  }
#endif
  bool success = true;
  if (comm->getRank() == 0) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    success = static_cast<bool>(in);
    if (success) {
      in.read(reinterpret_cast<char*>(&header), sizeof(FileHeader));
      success = static_cast<bool>(in);
    }
  }
  checkFileOpen(success, filename, comm, "reading");
  Teuchos::broadcast(*comm, 0, static_cast<int>(sizeof(FileHeader)), reinterpret_cast<char*>(&header));
  return header;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MapSectionHeader
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::makeMapSectionHeader(const map_type& map) {
  MapSectionHeader header;
  header.numGlobalElements = static_cast<unsigned long long>(map.getGlobalNumElements());
  header.indexBase         = static_cast<long long>(map.getIndexBase());
  header.mapFlags          = mapFlags(map);
  header.numRanks          = static_cast<unsigned long long>(map.getComm()->getSize());
  return header;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapSectionSize(const map_type& map) {
  return static_cast<unsigned long long>(sizeof(MapSectionHeader)) +
         static_cast<unsigned long long>(map.getComm()->getSize()) * sizeof(unsigned long long);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapSectionCountsOffset(const unsigned long long mapSectionOffset) {
  return mapSectionOffset + static_cast<unsigned long long>(sizeof(MapSectionHeader));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapSectionPayloadOffset(const unsigned long long mapSectionOffset,
                                                                                                const unsigned long long numRanks) {
  return mapSectionCountsOffset(mapSectionOffset) + numRanks * sizeof(unsigned long long);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeMapSectionHeader(const std::string& filename,
                                                                                const unsigned long long mapSectionOffset,
                                                                                const MapSectionHeader& header,
                                                                                const trcp_tcomm_t& comm) {
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while writing map section header.");
    if (comm->getRank() == 0) {
      MPI_Status status;
      const int writeErr = MPI_File_write_at(file, static_cast<MPI_Offset>(mapSectionOffset),
                                             const_cast<MapSectionHeader*>(&header), sizeof(MapSectionHeader), MPI_BYTE, &status);
      TEUCHOS_TEST_FOR_EXCEPTION(writeErr != MPI_SUCCESS, std::runtime_error,
                                 "Tpetra::BinaryIO: MPI_File_write_at failed while writing map section header.");
    }
    MPI_File_close(&file);
    comm->barrier();
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::fstream out(filename.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to reopen file '" << filename << "' while writing map section header.");
    out.seekp(static_cast<std::streamoff>(mapSectionOffset), std::ios::beg);
    out.write(reinterpret_cast<const char*>(&header), sizeof(MapSectionHeader));
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to write map section header to file '" << filename << "'.");
  }
  comm->barrier();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MapSectionHeader
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readMapSectionHeader(const std::string& filename,
                                                                          const unsigned long long mapSectionOffset,
                                                                          const trcp_tcomm_t& comm) {
  MapSectionHeader header;
  std::memset(&header, 0, sizeof(MapSectionHeader));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading map section header.");
    MPI_Status status;
    const int readErr = MPI_File_read_at_all(file, static_cast<MPI_Offset>(mapSectionOffset),
                                             &header, sizeof(MapSectionHeader), MPI_BYTE, &status);
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(readErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_read_at_all failed while reading map section header.");
    return header;
  }
#endif
  if (comm->getRank() == 0) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while reading map section header.");
    in.seekg(static_cast<std::streamoff>(mapSectionOffset), std::ios::beg);
    in.read(reinterpret_cast<char*>(&header), sizeof(MapSectionHeader));
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to read map section header from file '" << filename << "'.");
  }
  Teuchos::broadcast(*comm, 0, static_cast<int>(sizeof(MapSectionHeader)), reinterpret_cast<char*>(&header));
  return header;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::exclusiveScanUnsignedLongLong(const unsigned long long localValue,
                                                                                                      const trcp_tcomm_t& comm) {
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm          = Details::extractMpiCommFromTeuchos(*comm);
    unsigned long long result = 0;
    const int err             = MPI_Exscan(const_cast<unsigned long long*>(&localValue), &result, 1,
                                           MPI_UNSIGNED_LONG_LONG, MPI_SUM, rawComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_Exscan failed.");
    if (comm->getRank() == 0) {
      result = 0;
    }
    return result;
  }
#endif
  unsigned long long inclusive = 0;
  Teuchos::scan<int, unsigned long long>(*comm, Teuchos::REDUCE_SUM, 1,
                                         const_cast<unsigned long long*>(&localValue), &inclusive);
  return inclusive - localValue;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeArrayFromRoot(const std::string& filename,
                                                                             const unsigned long long dataOffset,
                                                                             const T* data,
                                                                             const unsigned long long count,
                                                                             const trcp_tcomm_t& comm) {
  const unsigned long long byteCount = checkedByteCount(count, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while writing root-owned array data.");
    if (comm->getRank() == 0) {
      MPI_Status status;
      const int writeErr = MPI_File_write_at(file,
                                             static_cast<MPI_Offset>(dataOffset),
                                             const_cast<T*>(data),
                                             static_cast<int>(byteCount),
                                             MPI_BYTE,
                                             &status);
      TEUCHOS_TEST_FOR_EXCEPTION(writeErr != MPI_SUCCESS, std::runtime_error,
                                 "Tpetra::BinaryIO: MPI_File_write_at failed while writing root-owned array data.");
    }
    MPI_File_close(&file);
    comm->barrier();
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::fstream out(filename.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while writing root-owned array data.");
    out.seekp(static_cast<std::streamoff>(dataOffset), std::ios::beg);
    out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(byteCount));
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to write root-owned array data to file '" << filename << "'.");
  }
  comm->barrier();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readArrayFromRoot(const std::string& filename,
                                                                            const unsigned long long dataOffset,
                                                                            T* data,
                                                                            const unsigned long long count,
                                                                            const trcp_tcomm_t& comm) {
  const unsigned long long byteCount = checkedByteCount(count, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading root-owned array data.");
    if (comm->getRank() == 0) {
      MPI_Status status;
      const int readErr = MPI_File_read_at(file,
                                           static_cast<MPI_Offset>(dataOffset),
                                           data,
                                           static_cast<int>(byteCount),
                                           MPI_BYTE,
                                           &status);
      TEUCHOS_TEST_FOR_EXCEPTION(readErr != MPI_SUCCESS, std::runtime_error,
                                 "Tpetra::BinaryIO: MPI_File_read_at failed while reading root-owned array data.");
    }
    MPI_File_close(&file);
    if (comm->getSize() > 1) {
      broadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
    }
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while reading root-owned array data.");
    in.seekg(static_cast<std::streamoff>(dataOffset), std::ios::beg);
    in.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(byteCount));
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to read root-owned array data from file '" << filename << "'.");
  }
  if (comm->getSize() > 1) {
    broadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeArrayCollective(const std::string& filename,
                                                                               const unsigned long long dataOffset,
                                                                               const T* data,
                                                                               const unsigned long long count,
                                                                               const unsigned long long globalOffset,
                                                                               const trcp_tcomm_t& comm) {
  const unsigned long long byteCount  = checkedByteCount(count, sizeof(T));
  const unsigned long long byteOffset = checkedByteOffset(dataOffset, globalOffset, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while writing array data.");

    const char* current                     = reinterpret_cast<const char*>(data);
    unsigned long long remaining            = byteCount;
    unsigned long long currentOffset        = byteOffset;
    const unsigned long long maxChunk       = static_cast<unsigned long long>(std::numeric_limits<int>::max());
    const unsigned long long localNumChunks = (byteCount + maxChunk - 1ull) / maxChunk;
    unsigned long long numChunks            = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localNumChunks, &numChunks);
    int writeErr = MPI_SUCCESS;
    for (unsigned long long chunkIndex = 0; chunkIndex < numChunks && writeErr == MPI_SUCCESS; ++chunkIndex) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      MPI_Status status;
      writeErr = MPI_File_write_at_all(file,
                                       static_cast<MPI_Offset>(currentOffset),
                                       const_cast<char*>(current),
                                       static_cast<int>(chunk),
                                       MPI_BYTE,
                                       &status);
      current += static_cast<std::ptrdiff_t>(chunk);
      currentOffset += chunk;
      remaining -= chunk;
    }
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(writeErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_write_at_all failed while writing array data.");
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::fstream out(filename.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while writing array data.");
    out.seekp(static_cast<std::streamoff>(byteOffset), std::ios::beg);

    const char* current               = reinterpret_cast<const char*>(data);
    unsigned long long remaining      = byteCount;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<std::streamsize>::max());
    while (remaining > 0) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      out.write(current, static_cast<std::streamsize>(chunk));
      TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                                 "Tpetra::BinaryIO: Failed to write array data to file '" << filename << "'.");
      current += static_cast<std::ptrdiff_t>(chunk);
      remaining -= chunk;
    }
  }
  comm->barrier();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
template <class T>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readArrayCollective(const std::string& filename,
                                                                              const unsigned long long dataOffset,
                                                                              T* data,
                                                                              const unsigned long long count,
                                                                              const unsigned long long globalOffset,
                                                                              const trcp_tcomm_t& comm) {
  const unsigned long long byteCount  = checkedByteCount(count, sizeof(T));
  const unsigned long long byteOffset = checkedByteOffset(dataOffset, globalOffset, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (Details::teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = Details::extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading array data.");

    char* current                           = reinterpret_cast<char*>(data);
    unsigned long long remaining            = byteCount;
    unsigned long long currentOffset        = byteOffset;
    const unsigned long long maxChunk       = static_cast<unsigned long long>(std::numeric_limits<int>::max());
    const unsigned long long localNumChunks = (byteCount + maxChunk - 1ull) / maxChunk;
    unsigned long long numChunks            = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localNumChunks, &numChunks);
    int readErr = MPI_SUCCESS;
    for (unsigned long long chunkIndex = 0; chunkIndex < numChunks && readErr == MPI_SUCCESS; ++chunkIndex) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      MPI_Status status;
      readErr = MPI_File_read_at_all(file,
                                     static_cast<MPI_Offset>(currentOffset),
                                     current,
                                     static_cast<int>(chunk),
                                     MPI_BYTE,
                                     &status);
      current += static_cast<std::ptrdiff_t>(chunk);
      currentOffset += chunk;
      remaining -= chunk;
    }
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(readErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_read_at_all failed while reading array data.");
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::ifstream in(filename.c_str(), std::ios::binary);
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while reading array data.");
    in.seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);

    char* current                     = reinterpret_cast<char*>(data);
    unsigned long long remaining      = byteCount;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<std::streamsize>::max());
    while (remaining > 0) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      in.read(current, static_cast<std::streamsize>(chunk));
      TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                                 "Tpetra::BinaryIO: Failed to read array data from file '" << filename << "'.");
      current += static_cast<std::ptrdiff_t>(chunk);
      remaining -= chunk;
    }
  }
  if (comm->getSize() > 1) {
    broadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeMapSection(const std::string& filename,
                                                                          const unsigned long long mapSectionOffset,
                                                                          const map_type& map,
                                                                          const trcp_tcomm_t& comm) {
  const MapSectionHeader sectionHeader = makeMapSectionHeader(map);
  writeMapSectionHeader(filename, mapSectionOffset, sectionHeader, comm);

  const auto gidsDevice               = map.getMyGlobalIndicesDevice();
  const unsigned long long localCount = static_cast<unsigned long long>(gidsDevice.extent(0));
  Kokkos::View<unsigned long long*, Kokkos::HostSpace> localCounts("Tpetra::BinaryIO::mapLocalCounts", comm->getSize());
  Teuchos::gatherAll(*comm, 1, &localCount, comm->getSize(), localCounts.data());
  writeArrayFromRoot(filename,
                     mapSectionCountsOffset(mapSectionOffset),
                     localCounts.data(),
                     static_cast<unsigned long long>(localCounts.extent(0)),
                     comm);

  const unsigned long long globalOffset = exclusiveScanUnsignedLongLong(localCount, comm);
  const auto gidsHost                   = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gidsDevice);
  writeArrayCollective(filename,
                       mapSectionPayloadOffset(mapSectionOffset, sectionHeader.numRanks),
                       gidsHost.data(),
                       localCount,
                       globalOffset,
                       comm);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readMapSection(const std::string& filename,
                                                                    const unsigned long long mapSectionOffset,
                                                                    const trcp_tcomm_t& comm) {
  const MapSectionHeader sectionHeader = readMapSectionHeader(filename, mapSectionOffset, comm);
  const unsigned long long globalCount = sectionHeader.numGlobalElements;

  TEUCHOS_TEST_FOR_EXCEPTION(sectionHeader.numRanks != static_cast<unsigned long long>(comm->getSize()),
                             std::runtime_error,
                             "Tpetra::BinaryIO: Map section was written for " << sectionHeader.numRanks
                                                                              << " ranks, but the read communicator has "
                                                                              << comm->getSize() << " ranks.");

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> localCounts("Tpetra::BinaryIO::mapLocalCounts",
                                                                   Details::binaryIOCheckedSize(sectionHeader.numRanks, "map section rank count"));
  readArrayFromRoot(filename,
                    mapSectionCountsOffset(mapSectionOffset),
                    localCounts.data(),
                    static_cast<unsigned long long>(localCounts.extent(0)),
                    comm);

  const unsigned long long localCount = localCounts(comm->getRank());
  unsigned long long globalOffset     = 0;
  for (int rank = 0; rank < comm->getRank(); ++rank) {
    globalOffset += localCounts(rank);
  }

  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> gids("Tpetra::BinaryIO::mapGids",
                                                       Details::binaryIOCheckedSize(localCount, "map section local element count"));
  if (localCount > 0) {
    readArrayCollective(filename,
                        mapSectionPayloadOffset(mapSectionOffset, sectionHeader.numRanks),
                        gids.data(),
                        localCount,
                        globalOffset,
                        comm);
  }

  using device_type = typename map_type::device_type;
  Kokkos::View<GlobalOrdinal*, device_type> gidsDevice("Tpetra::BinaryIO::mapGidsDevice", gids.extent(0));
  Kokkos::deep_copy(gidsDevice, gids);

  return Teuchos::rcp(new map_type(static_cast<Tpetra::global_size_t>(globalCount),
                                   gidsDevice,
                                   static_cast<GlobalOrdinal>(sectionHeader.indexBase),
                                   comm));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeMapFile(const std::string& filename,
                                                                       const map_type& map) {
  FileHeader header    = makeBaseHeader(MAP_OBJECT);
  header.numGlobalRows = static_cast<unsigned long long>(map.getGlobalNumElements());
  header.rowMapOffset  = static_cast<unsigned long long>(sizeof(FileHeader));
  writeHeaderToNewFile(filename, header, map.getComm());
  writeMapSection(filename, header.rowMapOffset, map, map.getComm());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readMapFile(const std::string& filename,
                                                                 const trcp_tcomm_t& comm) {
  return Tpetra::readBinaryMapFile<LocalOrdinal, GlobalOrdinal, Node>(filename, comm);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeDenseFile(const std::string& filename,
                                                                         const multivector_type& X) {
  const auto map       = X.getMap();
  FileHeader header    = makeBaseHeader(MULTIVECTOR_OBJECT);
  header.numGlobalRows = static_cast<unsigned long long>(X.getGlobalLength());
  header.numGlobalCols = static_cast<unsigned long long>(X.getNumVectors());
  header.numVectors    = static_cast<unsigned long long>(X.getNumVectors());
  header.rowMapOffset  = static_cast<unsigned long long>(sizeof(FileHeader));
  header.valuesOffset  = header.rowMapOffset + mapSectionSize(*map) + header.numGlobalRows * sizeof(GlobalOrdinal);
  writeHeaderToNewFile(filename, header, map->getComm());
  writeMapSection(filename, header.rowMapOffset, *map, map->getComm());

  const unsigned long long localCount   = static_cast<unsigned long long>(map->getLocalNumElements());
  const unsigned long long globalOffset = exclusiveScanUnsignedLongLong(localCount, map->getComm());
  const unsigned long long globalLength = static_cast<unsigned long long>(X.getGlobalLength());

  using host_execution_space = Kokkos::DefaultHostExecutionSpace;
  auto localXDevice          = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto localXHost            = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localXDevice);
  Kokkos::View<Scalar*, Kokkos::HostSpace> localValues("Tpetra::BinaryIO::denseLocalValues", Details::binaryIOCheckedSize(localCount, "dense local row count"));
  for (size_t j = 0; j < X.getNumVectors(); ++j) {
    const unsigned long long columnOffset = static_cast<unsigned long long>(j) * globalLength + globalOffset;
    Kokkos::parallel_for(
        "Tpetra::BinaryIO::copyDenseColumnToHost",
        Kokkos::RangePolicy<host_execution_space>(0, localValues.extent(0)),
        KOKKOS_LAMBDA(const size_t i) {
          localValues(i) = static_cast<Scalar>(localXHost(i, j));
        });
    writeArrayCollective(filename, header.valuesOffset, localValues.data(), localCount, columnOffset, map->getComm());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::multivector_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readDenseFile(const std::string& filename,
                                                                   const trcp_tcomm_t& comm) {
  const FileHeader header = readHeaderFromFile(filename, comm);
  validateHeader(header, MULTIVECTOR_OBJECT);
  auto map = readMapSection(filename, header.rowMapOffset, comm);
  return readDenseFile(filename, map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::multivector_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readDenseFile(const std::string& filename,
                                                                   const Teuchos::RCP<const map_type>& map) {
  const auto comm         = map->getComm();
  const FileHeader header = readHeaderFromFile(filename, comm);
  validateHeader(header, MULTIVECTOR_OBJECT);

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<unsigned long long>(map->getGlobalNumElements()) != header.numGlobalRows,
                             std::invalid_argument,
                             "Tpetra::BinaryIO: Requested dense read map has "
                                 << map->getGlobalNumElements() << " global elements, but file has "
                                 << header.numGlobalRows << " rows.");

  const unsigned long long globalLength = header.numGlobalRows;
  auto fileMap                          = readMapSection(filename, header.rowMapOffset, comm);
  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<unsigned long long>(fileMap->getGlobalNumElements()) != globalLength,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Dense file row-map length does not match header row count.");

  auto fileX                                = Teuchos::rcp(new multivector_type(fileMap, static_cast<size_t>(header.numVectors)));
  const unsigned long long fileLocalCount   = static_cast<unsigned long long>(fileMap->getLocalNumElements());
  const unsigned long long fileGlobalOffset = exclusiveScanUnsignedLongLong(fileLocalCount, comm);
  using impl_scalar_type                    = typename multivector_type::impl_scalar_type;
  using host_execution_space                = Kokkos::DefaultHostExecutionSpace;
  auto localFileXDevice                     = fileX->getLocalViewDevice(Tpetra::Access::OverwriteAll);
  auto localFileXHost                       = Kokkos::create_mirror_view(localFileXDevice);
  Kokkos::View<Scalar*, Kokkos::HostSpace> localValues("Tpetra::BinaryIO::denseLocalValues", Details::binaryIOCheckedSize(fileLocalCount, "dense local row count"));

  for (size_t j = 0; j < fileX->getNumVectors(); ++j) {
    const unsigned long long columnOffset = static_cast<unsigned long long>(j) * globalLength + fileGlobalOffset;
    if (fileLocalCount > 0) {
      readArrayCollective(filename, header.valuesOffset, localValues.data(), fileLocalCount, columnOffset, comm);
    }
    Kokkos::parallel_for(
        "Tpetra::BinaryIO::copyDenseColumnToDeviceHost",
        Kokkos::RangePolicy<host_execution_space>(0, localValues.extent(0)),
        KOKKOS_LAMBDA(const size_t i) {
          localFileXHost(i, j) = static_cast<impl_scalar_type>(localValues(i));
        });
  }
  Kokkos::deep_copy(localFileXDevice, localFileXHost);

  if (map->isSameAs(*fileMap)) {
    return fileX;
  }

  using import_type = Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  auto X            = Teuchos::rcp(new multivector_type(map, static_cast<size_t>(header.numVectors)));
  import_type importer(fileMap, map);
  X->doImport(*fileX, importer, Tpetra::INSERT);
  return X;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeVectorFile(const std::string& filename,
                                                                          const vector_type& X) {
  writeDenseFile(filename, static_cast<const multivector_type&>(X));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vector_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readVectorFile(const std::string& filename,
                                                                    const trcp_tcomm_t& comm) {
  auto X = readDenseFile(filename, comm);
  TEUCHOS_TEST_FOR_EXCEPTION(X->getNumVectors() != 1, std::runtime_error,
                             "Tpetra::BinaryIO: File does not contain a single vector.");
  return X->getVectorNonConst(0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vector_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readVectorFile(const std::string& filename,
                                                                    const Teuchos::RCP<const map_type>& map) {
  auto X = readDenseFile(filename, map);
  TEUCHOS_TEST_FOR_EXCEPTION(X->getNumVectors() != 1, std::runtime_error,
                             "Tpetra::BinaryIO: File does not contain a single vector.");
  return X->getVectorNonConst(0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeSparseFile(const std::string& filename,
                                                                          const sparse_matrix_type& A) {
  TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), std::invalid_argument,
                             "Tpetra::BinaryIO::writeSparseFile requires a fill-complete matrix.");

  const auto rowMap    = A.getRowMap();
  const auto domainMap = A.getDomainMap();
  const auto rangeMap  = A.getRangeMap();
  const auto colMap    = A.getColMap();
  TEUCHOS_TEST_FOR_EXCEPTION(colMap.is_null(), std::invalid_argument,
                             "Tpetra::BinaryIO::writeSparseFile requires a nonnull column map.");

  const unsigned long long globalNumRows    = static_cast<unsigned long long>(A.getGlobalNumRows());
  const unsigned long long globalNumEntries = static_cast<unsigned long long>(A.getGlobalNumEntries());
  const unsigned long long rowMapBytes      = mapSectionSize(*rowMap) + globalNumRows * sizeof(GlobalOrdinal);
  const unsigned long long domainMapBytes   = mapSectionSize(*domainMap) + static_cast<unsigned long long>(domainMap->getGlobalNumElements()) * sizeof(GlobalOrdinal);
  const unsigned long long rangeMapBytes    = mapSectionSize(*rangeMap) + static_cast<unsigned long long>(rangeMap->getGlobalNumElements()) * sizeof(GlobalOrdinal);

  FileHeader header          = makeBaseHeader(MATRIX_OBJECT);
  header.numGlobalRows       = globalNumRows;
  header.numGlobalCols       = static_cast<unsigned long long>(A.getGlobalNumCols());
  header.numGlobalEntries    = globalNumEntries;
  header.rowMapOffset        = static_cast<unsigned long long>(sizeof(FileHeader));
  header.domainMapOffset     = header.rowMapOffset + rowMapBytes;
  header.rangeMapOffset      = header.domainMapOffset + domainMapBytes;
  header.rowPtrOffset        = header.rangeMapOffset + rangeMapBytes;
  header.columnIndicesOffset = header.rowPtrOffset + (globalNumRows + 1ull) * sizeof(unsigned long long);
  header.valuesOffset        = header.columnIndicesOffset + globalNumEntries * sizeof(GlobalOrdinal);

  writeHeaderToNewFile(filename, header, rowMap->getComm());
  writeMapSection(filename, header.rowMapOffset, *rowMap, rowMap->getComm());
  writeMapSection(filename, header.domainMapOffset, *domainMap, domainMap->getComm());
  writeMapSection(filename, header.rangeMapOffset, *rangeMap, rangeMap->getComm());

  using device_type          = typename sparse_matrix_type::device_type;
  using execution_space      = typename device_type::execution_space;
  using host_execution_space = Kokkos::DefaultHostExecutionSpace;

  const size_t localNumRows = rowMap->getLocalNumElements();
  const auto localMatrix    = A.getLocalMatrixDevice();
  const auto rowPtrDevice   = localMatrix.graph.row_map;
  auto rowPtrDeviceHost     = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowPtrDevice);

  const unsigned long long localNnz = static_cast<unsigned long long>(rowPtrDeviceHost(localNumRows));
  const size_t localNnzSize         = Details::binaryIOCheckedSize(localNnz, "local sparse entry count");
  Kokkos::View<unsigned long long*, Kokkos::HostSpace> localRowPtr("Tpetra::BinaryIO::localRowPtr", localNumRows + 1);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::copyLocalRowPtr",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        localRowPtr(i) = static_cast<unsigned long long>(rowPtrDeviceHost(i));
      });

  Kokkos::View<GlobalOrdinal*, device_type> globalColumnsDevice("Tpetra::BinaryIO::globalColumns", localNnzSize);
  const auto localColumnsDevice = localMatrix.graph.entries;
  const auto localColMap        = colMap->getLocalMap();
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::convertColumnsToGlobal",
      Kokkos::RangePolicy<execution_space>(0, localNnzSize),
      KOKKOS_LAMBDA(const size_t i) {
        globalColumnsDevice(i) = localColMap.getGlobalElement(localColumnsDevice(i));
      });
  auto localColumns = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), globalColumnsDevice);

  auto localValuesImpl = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), localMatrix.values);
  Kokkos::View<Scalar*, Kokkos::HostSpace> localValues("Tpetra::BinaryIO::localValues", localNnzSize);
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::copyLocalValues",
      Kokkos::RangePolicy<host_execution_space>(0, localNnzSize),
      KOKKOS_LAMBDA(const size_t i) {
        localValues(i) = static_cast<Scalar>(localValuesImpl(i));
      });

  const unsigned long long globalRowOffset = exclusiveScanUnsignedLongLong(static_cast<unsigned long long>(localNumRows), rowMap->getComm());
  const unsigned long long globalNnzOffset = exclusiveScanUnsignedLongLong(localNnz, rowMap->getComm());
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::offsetLocalRowPtr",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        localRowPtr(i) += globalNnzOffset;
      });

  writeArrayCollective(filename, header.rowPtrOffset, localRowPtr.data(),
                       static_cast<unsigned long long>(localRowPtr.extent(0)), globalRowOffset, rowMap->getComm());
  writeArrayCollective(filename, header.columnIndicesOffset, localColumns.data(), localNnz, globalNnzOffset, rowMap->getComm());
  writeArrayCollective(filename, header.valuesOffset, localValues.data(), localNnz, globalNnzOffset, rowMap->getComm());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sparse_matrix_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readSparseFile(const std::string& filename,
                                                                    const trcp_tcomm_t& comm) {
  const FileHeader header = readHeaderFromFile(filename, comm);
  validateHeader(header, MATRIX_OBJECT);
  auto rowMap    = readMapSection(filename, header.rowMapOffset, comm);
  auto domainMap = readMapSection(filename, header.domainMapOffset, comm);
  auto rangeMap  = readMapSection(filename, header.rangeMapOffset, comm);
  return readSparseFile(filename, rowMap, Teuchos::null, domainMap, rangeMap, true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sparse_matrix_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readSparseFile(const std::string& filename,
                                                                    const Teuchos::RCP<const map_type>& rowMap,
                                                                    const Teuchos::RCP<const map_type>& colMap,
                                                                    const Teuchos::RCP<const map_type>& domainMap,
                                                                    const Teuchos::RCP<const map_type>& rangeMap,
                                                                    const bool callFillComplete) {
  const auto comm         = rowMap->getComm();
  const FileHeader header = readHeaderFromFile(filename, comm);
  validateHeader(header, MATRIX_OBJECT);

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<unsigned long long>(rowMap->getGlobalNumElements()) != header.numGlobalRows,
                             std::invalid_argument,
                             "Tpetra::BinaryIO: Requested sparse read row map has "
                                 << rowMap->getGlobalNumElements() << " global elements, but file has "
                                 << header.numGlobalRows << " rows.");

  const auto fileRowMap    = readMapSection(filename, header.rowMapOffset, comm);
  const auto fileDomainMap = readMapSection(filename, header.domainMapOffset, comm);
  const auto fileRangeMap  = readMapSection(filename, header.rangeMapOffset, comm);

  using host_execution_space = Kokkos::DefaultHostExecutionSpace;

  const size_t localNumRows                = fileRowMap->getLocalNumElements();
  const unsigned long long localNumRows64  = static_cast<unsigned long long>(localNumRows);
  const unsigned long long globalRowOffset = exclusiveScanUnsignedLongLong(localNumRows64, comm);

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> localRowPtr("Tpetra::BinaryIO::localRowPtr", localNumRows + 1);
  if (localNumRows > 0) {
    readArrayCollective(filename, header.rowPtrOffset, localRowPtr.data(),
                        localNumRows64 + 1, globalRowOffset, comm);
    TEUCHOS_TEST_FOR_EXCEPTION(localRowPtr(localNumRows) < localRowPtr(0),
                               std::runtime_error,
                               "Tpetra::BinaryIO: Sparse row pointers in file are not monotonic on this rank.");
  } else {
    localRowPtr(0) = 0;
  }
  const unsigned long long nnzStart          = localRowPtr(0);
  unsigned long long nonmonotonicRowPtrCount = 0;
  Kokkos::parallel_reduce(
      "Tpetra::BinaryIO::checkReadRowPtr",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows),
      KOKKOS_LAMBDA(const size_t i, unsigned long long& count) {
        if (localRowPtr(i + 1) < localRowPtr(i)) {
          ++count;
        }
      },
      nonmonotonicRowPtrCount);
  TEUCHOS_TEST_FOR_EXCEPTION(nonmonotonicRowPtrCount != 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Sparse row pointers in file are not monotonic on this rank.");
  Kokkos::parallel_for(
      "Tpetra::BinaryIO::normalizeReadRowPtr",
      Kokkos::RangePolicy<host_execution_space>(0, localNumRows + 1),
      KOKKOS_LAMBDA(const size_t i) {
        localRowPtr(i) -= nnzStart;
      });

  const unsigned long long localNnz64 = localRowPtr(localNumRows);
  const size_t localNnz               = Details::binaryIOCheckedSize(localNnz64, "local sparse entry count");
  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> globalColumns("Tpetra::BinaryIO::globalColumns", localNnz);
  Kokkos::View<Scalar*, Kokkos::HostSpace> values("Tpetra::BinaryIO::values", localNnz);
  if (localNnz64 > 0) {
    readArrayCollective(filename, header.columnIndicesOffset, globalColumns.data(), localNnz64, nnzStart, comm);
    readArrayCollective(filename, header.valuesOffset, values.data(), localNnz64, nnzStart, comm);
  }

  const auto dom = domainMap.is_null() ? fileDomainMap : domainMap;
  const auto ran = rangeMap.is_null() ? fileRangeMap : rangeMap;

  if (rowMap->isSameAs(*fileRowMap)) {
    auto matrix = Details::buildSparseMatrixFromLocalCrsViews<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rowMap, colMap, dom, localRowPtr, globalColumns, values);
    if (callFillComplete) {
      matrix->fillComplete(dom, ran);
    }
    return matrix;
  }

  auto fileMatrix = Details::buildSparseMatrixFromLocalCrsViews<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fileRowMap, Teuchos::null, fileDomainMap, localRowPtr, globalColumns, values);
  fileMatrix->fillComplete(fileDomainMap, fileRangeMap);

  using import_type = Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  import_type importer(fileRowMap, rowMap);

  if (callFillComplete && colMap.is_null()) {
    return Tpetra::importAndFillCompleteCrsMatrix<sparse_matrix_type>(fileMatrix,
                                                                      importer,
                                                                      dom,
                                                                      ran);
  }

  Teuchos::RCP<sparse_matrix_type> matrix;
  if (colMap.is_null()) {
    matrix = Teuchos::rcp(new sparse_matrix_type(rowMap,
                                                 static_cast<size_t>(fileMatrix->getGlobalMaxNumRowEntries())));
  } else {
    matrix = Teuchos::rcp(new sparse_matrix_type(rowMap,
                                                 colMap,
                                                 static_cast<size_t>(fileMatrix->getGlobalMaxNumRowEntries())));
  }
  matrix->doImport(*fileMatrix, importer, Tpetra::INSERT);
  if (callFillComplete) {
    matrix->fillComplete(dom, ran);
  }
  return matrix;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sparse_matrix_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readSparseFile(const std::string& filename,
                                                                    const Teuchos::RCP<const map_type>& rowMap,
                                                                    const Teuchos::RCP<const map_type>& domainMap,
                                                                    const Teuchos::RCP<const map_type>& rangeMap,
                                                                    const bool callFillComplete) {
  return readSparseFile(filename, rowMap, Teuchos::null, domainMap, rangeMap, callFillComplete);
}

}  // namespace Tpetra

//
// Explicit template instantiation macro for the Tpetra::BinaryIO class.
// BinaryIO is templated on Scalar, LocalOrdinal, GlobalOrdinal, and Node (the
// SC-LO-GO-Node scheme).  The scalar-free Tpetra::readBinaryMapFile reader is
// instantiated separately (see Tpetra_ReadBinaryMapFile_def.hpp); BinaryIO only
// sees its declaration, so no cross-scheme extern-template coupling is needed.
//
#define TPETRA_BINARYIO_INSTANT(SCALAR, LO, GO, NODE) \
  template class BinaryIO<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_BINARYIO_DEF_HPP
