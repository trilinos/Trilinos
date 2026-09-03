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

#include "Tpetra_Details_MpiTypeTraits.hpp"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace Tpetra {

namespace Details {

template <class GO, class Node>
Teuchos::RCP<const Tpetra::Map<int, GO, Node>>
makeContiguousMapForRead(const unsigned long long globalCount,
                         const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  using map_type         = Tpetra::Map<int, GO, Node>;
  using global_size_type = typename map_type::global_size_t;
  return Teuchos::rcp(new map_type(static_cast<global_size_type>(globalCount), static_cast<GO>(0), comm));
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
                             "Tpetra::BinaryIO: File does not have a recognized Tpetra binary I/O header.");
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
  TEUCHOS_TEST_FOR_EXCEPTION(header.scalarSize != sizeof(Scalar) || header.scalarFlags != typeFlags<Scalar>(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File scalar type does not match this BinaryIO instantiation.");
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
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::computeContiguousLocalCount(const unsigned long long globalCount,
                                                                                                    const int rank,
                                                                                                    const int size) {
  const unsigned long long quotient  = globalCount / static_cast<unsigned long long>(size);
  const unsigned long long remainder = globalCount % static_cast<unsigned long long>(size);
  return quotient + (static_cast<unsigned long long>(rank) < remainder ? 1ull : 0ull);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::computeContiguousStart(const unsigned long long globalCount,
                                                                                               const int rank,
                                                                                               const int size) {
  const unsigned long long quotient  = globalCount / static_cast<unsigned long long>(size);
  const unsigned long long remainder = globalCount % static_cast<unsigned long long>(size);
  const unsigned long long rankU     = static_cast<unsigned long long>(rank);
  return rankU * quotient + std::min(rankU, remainder);
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
  return header;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapSectionSize() {
  return static_cast<unsigned long long>(sizeof(MapSectionHeader));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
unsigned long long BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mapSectionPayloadOffset(const unsigned long long mapSectionOffset) {
  return mapSectionOffset + mapSectionSize();
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

    const char* current               = reinterpret_cast<const char*>(data);
    unsigned long long remaining      = byteCount;
    unsigned long long currentOffset  = byteOffset;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<int>::max());
    int writeErr                      = MPI_SUCCESS;
    while (remaining > 0 && writeErr == MPI_SUCCESS) {
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

    char* current                     = reinterpret_cast<char*>(data);
    unsigned long long remaining      = byteCount;
    unsigned long long currentOffset  = byteOffset;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<int>::max());
    int readErr                       = MPI_SUCCESS;
    while (remaining > 0 && readErr == MPI_SUCCESS) {
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

  auto gids                             = map.getLocalElementList();
  const unsigned long long localCount   = static_cast<unsigned long long>(gids.size());
  const unsigned long long globalOffset = exclusiveScanUnsignedLongLong(localCount, comm);
  writeArrayCollective(filename,
                       mapSectionPayloadOffset(mapSectionOffset),
                       gids.getRawPtr(),
                       localCount,
                       globalOffset,
                       comm);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::map_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readMapSection(const std::string& filename,
                                                                    const unsigned long long mapSectionOffset,
                                                                    const trcp_tcomm_t& comm) {
  const MapSectionHeader sectionHeader  = readMapSectionHeader(filename, mapSectionOffset, comm);
  const unsigned long long globalCount  = sectionHeader.numGlobalElements;
  const unsigned long long localCount   = computeContiguousLocalCount(globalCount, comm->getRank(), comm->getSize());
  const unsigned long long globalOffset = computeContiguousStart(globalCount, comm->getRank(), comm->getSize());

  Teuchos::Array<GlobalOrdinal> gids(static_cast<Teuchos::Array<GlobalOrdinal>::size_type>(localCount));
  if (localCount > 0) {
    readArrayCollective(filename,
                        mapSectionPayloadOffset(mapSectionOffset),
                        gids.getRawPtr(),
                        localCount,
                        globalOffset,
                        comm);
  }

  return Teuchos::rcp(new map_type(static_cast<Tpetra::global_size_t>(globalCount),
                                   gids(),
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
  const FileHeader header = readHeaderFromFile(filename, comm);
  validateHeader(header, MAP_OBJECT);
  return readMapSection(filename, header.rowMapOffset, comm);
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
  header.valuesOffset  = header.rowMapOffset + mapSectionSize() + header.numGlobalRows * sizeof(GlobalOrdinal);
  writeHeaderToNewFile(filename, header, map->getComm());
  writeMapSection(filename, header.rowMapOffset, *map, map->getComm());

  const unsigned long long localCount   = static_cast<unsigned long long>(map->getLocalNumElements());
  const unsigned long long globalOffset = exclusiveScanUnsignedLongLong(localCount, map->getComm());
  const unsigned long long globalLength = static_cast<unsigned long long>(X.getGlobalLength());

  for (size_t j = 0; j < X.getNumVectors(); ++j) {
    const auto data                       = X.getData(j);
    const unsigned long long columnOffset = static_cast<unsigned long long>(j) * globalLength + globalOffset;
    writeArrayCollective(filename, header.valuesOffset, data.getRawPtr(), localCount, columnOffset, map->getComm());
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

  for (size_t j = 0; j < fileX->getNumVectors(); ++j) {
    auto data                             = fileX->getDataNonConst(j);
    const unsigned long long columnOffset = static_cast<unsigned long long>(j) * globalLength + fileGlobalOffset;
    if (fileLocalCount > 0) {
      readArrayCollective(filename, header.valuesOffset, data.getRawPtr(), fileLocalCount, columnOffset, comm);
    }
  }

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
  const unsigned long long rowMapBytes      = mapSectionSize() + globalNumRows * sizeof(GlobalOrdinal);
  const unsigned long long domainMapBytes   = mapSectionSize() + static_cast<unsigned long long>(domainMap->getGlobalNumElements()) * sizeof(GlobalOrdinal);
  const unsigned long long rangeMapBytes    = mapSectionSize() + static_cast<unsigned long long>(rangeMap->getGlobalNumElements()) * sizeof(GlobalOrdinal);

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

  const size_t localNumRows = rowMap->getLocalNumElements();
  Teuchos::Array<unsigned long long> localRowPtr(localNumRows + 1);
  localRowPtr[0] = 0;
  std::vector<GlobalOrdinal> localColumns;
  std::vector<Scalar> localValues;
  localColumns.reserve(static_cast<size_t>(A.getLocalNumEntries()));
  localValues.reserve(static_cast<size_t>(A.getLocalNumEntries()));

  for (size_t lclRow = 0; lclRow < localNumRows; ++lclRow) {
    const GlobalOrdinal gblRow = rowMap->getGlobalElement(static_cast<LocalOrdinal>(lclRow));
    const size_t numEntries    = A.getNumEntriesInGlobalRow(gblRow);
    typename sparse_matrix_type::nonconst_global_inds_host_view_type indices("indices", numEntries);
    typename sparse_matrix_type::nonconst_values_host_view_type values("values", numEntries);
    size_t entriesRead = 0;
    A.getGlobalRowCopy(gblRow, indices, values, entriesRead);
    for (size_t k = 0; k < entriesRead; ++k) {
      localColumns.push_back(indices(k));
      localValues.push_back(values(k));
    }
    localRowPtr[lclRow + 1] = static_cast<unsigned long long>(localColumns.size());
  }

  const unsigned long long localNnz        = static_cast<unsigned long long>(localColumns.size());
  const unsigned long long globalRowOffset = exclusiveScanUnsignedLongLong(static_cast<unsigned long long>(localNumRows), rowMap->getComm());
  const unsigned long long globalNnzOffset = exclusiveScanUnsignedLongLong(localNnz, rowMap->getComm());
  for (size_t i = 0; i < localRowPtr.size(); ++i) {
    localRowPtr[i] += globalNnzOffset;
  }

  writeArrayCollective(filename, header.rowPtrOffset, localRowPtr.getRawPtr(),
                       static_cast<unsigned long long>(localRowPtr.size()), globalRowOffset, rowMap->getComm());
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
  return readSparseFile(filename, rowMap, domainMap, rangeMap, true);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sparse_matrix_type>
BinaryIO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::readSparseFile(const std::string& filename,
                                                                    const Teuchos::RCP<const map_type>& rowMap,
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

  const auto activeRowMap                  = fileRowMap;
  const unsigned long long localNumRows    = static_cast<unsigned long long>(activeRowMap->getLocalNumElements());
  const unsigned long long globalRowOffset = exclusiveScanUnsignedLongLong(localNumRows, comm);

  Teuchos::Array<unsigned long long> localRowPtr(localNumRows + 1);
  if (localNumRows > 0) {
    readArrayCollective(filename, header.rowPtrOffset, localRowPtr.getRawPtr(),
                        localNumRows + 1, globalRowOffset, comm);
  } else {
    localRowPtr[0] = 0;
  }
  const unsigned long long nnzStart = localRowPtr[0];
  for (size_t i = localRowPtr.size() - 1; i > 0; --i) {
    localRowPtr[i] -= nnzStart;
  }
  localRowPtr[0] = 0;

  const unsigned long long localNnz = localRowPtr[localNumRows];
  std::vector<GlobalOrdinal> globalColumns(static_cast<size_t>(localNnz));
  std::vector<Scalar> values(static_cast<size_t>(localNnz));
  if (localNnz > 0) {
    readArrayCollective(filename, header.columnIndicesOffset, globalColumns.data(), localNnz, nnzStart, comm);
    readArrayCollective(filename, header.valuesOffset, values.data(), localNnz, nnzStart, comm);
  }

  Teuchos::Array<size_t> rowLengths(static_cast<size_t>(localNumRows));
  for (size_t lclRow = 0; lclRow < static_cast<size_t>(localNumRows); ++lclRow) {
    rowLengths[lclRow] = static_cast<size_t>(localRowPtr[lclRow + 1] - localRowPtr[lclRow]);
  }

  auto fileMatrix = Teuchos::rcp(new sparse_matrix_type(fileRowMap, rowLengths()));
  for (size_t lclRow = 0; lclRow < static_cast<size_t>(localNumRows); ++lclRow) {
    const GlobalOrdinal gblRow = fileRowMap->getGlobalElement(static_cast<LocalOrdinal>(lclRow));
    const size_t rowBegin      = static_cast<size_t>(localRowPtr[lclRow]);
    const size_t rowEnd        = static_cast<size_t>(localRowPtr[lclRow + 1]);
    const size_t rowLen        = rowEnd - rowBegin;
    if (rowLen > 0) {
      Teuchos::ArrayView<const GlobalOrdinal> cols(globalColumns.data() + rowBegin, rowLen);
      Teuchos::ArrayView<const Scalar> vals(values.data() + rowBegin, rowLen);
      fileMatrix->insertGlobalValues(gblRow, cols, vals);
    }
  }

  if (rowMap->isSameAs(*fileRowMap)) {
    if (callFillComplete) {
      const auto dom = domainMap.is_null() ? fileDomainMap : domainMap;
      const auto ran = rangeMap.is_null() ? fileRangeMap : rangeMap;
      fileMatrix->fillComplete(dom, ran);
    }
    return fileMatrix;
  }

  fileMatrix->fillComplete(fileDomainMap, fileRangeMap);

  using import_type = Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  import_type importer(fileRowMap, rowMap);

  if (callFillComplete) {
    const auto dom = domainMap.is_null() ? fileDomainMap : domainMap;
    const auto ran = rangeMap.is_null() ? fileRangeMap : rangeMap;
    return Tpetra::importAndFillCompleteCrsMatrix<sparse_matrix_type>(fileMatrix,
                                                                       importer,
                                                                       dom,
                                                                       ran);
  }

  auto matrix = Teuchos::rcp(new sparse_matrix_type(rowMap,
                                                    static_cast<size_t>(fileMatrix->getGlobalMaxNumRowEntries())));
  matrix->doImport(*fileMatrix, importer, Tpetra::INSERT);
  return matrix;
}

}  // namespace Tpetra

#define TPETRA_BINARYIO_INSTANT(SCALAR, LO, GO, NODE) \
  template class BinaryIO<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_BINARYIO_DEF_HPP
