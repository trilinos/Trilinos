// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_READBINARYMAPFILE_DEF_HPP
#define TPETRA_READBINARYMAPFILE_DEF_HPP

#include "Tpetra_ReadBinaryMapFile_decl.hpp"

#include "Tpetra_BinaryIO_Helpers.hpp"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Kokkos_Core.hpp"

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <limits>
#include <type_traits>

namespace Tpetra {

namespace Details {

/// \brief On-disk header shared by every Tpetra binary I/O file.
///
/// This mirrors \c Tpetra::BinaryIO::FileHeader.  It is duplicated here (rather
/// than shared) so that the scalar-free map reader does not depend on the
/// Scalar-templated \c Tpetra::BinaryIO class.
struct BinaryIOFileHeader {
  char magic[8];
  unsigned long long version;
  unsigned long long byteOrderMarker;
  unsigned long long objectKind;
  unsigned long long scalarSize;
  unsigned long long scalarFlags;
  unsigned long long localOrdinalSize;
  unsigned long long localOrdinalFlags;
  unsigned long long globalOrdinalSize;
  unsigned long long globalOrdinalFlags;
  unsigned long long objectFlags;
  unsigned long long numGlobalRows;
  unsigned long long numGlobalCols;
  unsigned long long numGlobalEntries;
  unsigned long long numVectors;
  unsigned long long rowMapOffset;
  unsigned long long domainMapOffset;
  unsigned long long rangeMapOffset;
  unsigned long long rowPtrOffset;
  unsigned long long columnIndicesOffset;
  unsigned long long valuesOffset;
};

/// \brief On-disk header preceding each map section in a Tpetra binary file.
struct BinaryIOMapSectionHeader {
  unsigned long long numGlobalElements;
  long long indexBase;
  unsigned long long mapFlags;
  unsigned long long numRanks;
};

static const unsigned long long binaryIOFileVersion     = 2ull;
static const unsigned long long binaryIOByteOrderMarker = 0x0102030405060708ull;
static const unsigned long long binaryIOMapObject       = 1ull;

template <class T>
unsigned long long binaryIOTypeFlags() {
  unsigned long long flags = 0;
  if (Teuchos::ScalarTraits<T>::isOrdinal) {
    flags |= 1u << 0;
  }
  if (Teuchos::ScalarTraits<T>::isComplex) {
    flags |= 1u << 2;
  }
  if (std::is_floating_point<T>::value || Teuchos::ScalarTraits<T>::isComplex) {
    flags |= 1u << 1;
  }
  if (std::is_signed<T>::value || Teuchos::ScalarTraits<T>::isComplex || std::is_floating_point<T>::value) {
    flags |= 1u << 3;
  }
  return flags;
}

inline void binaryIOBroadcastBytesFromRoot(char data[],
                                           const unsigned long long byteCount,
                                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
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

inline void binaryIOCheckFileOpen(const bool success,
                                  const std::string& filename,
                                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                                  const char mode[]) {
  int opened = success ? 1 : 0;
  Teuchos::broadcast(*comm, 0, 1, &opened);
  TEUCHOS_TEST_FOR_EXCEPTION(opened == 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Failed to open file '" << filename << "' for " << mode << ".");
}

inline unsigned long long binaryIOCheckedByteCount(const unsigned long long count,
                                                   const size_t elementSize) {
  TEUCHOS_TEST_FOR_EXCEPTION(count > 0 &&
                                 elementSize > static_cast<size_t>(std::numeric_limits<unsigned long long>::max() / count),
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte count overflow while computing transfer size.");
  return count * static_cast<unsigned long long>(elementSize);
}

inline unsigned long long binaryIOCheckedByteOffset(const unsigned long long dataOffset,
                                                    const unsigned long long globalOffset,
                                                    const size_t elementSize) {
  const unsigned long long payloadOffset = binaryIOCheckedByteCount(globalOffset, elementSize);
  TEUCHOS_TEST_FOR_EXCEPTION(dataOffset > std::numeric_limits<unsigned long long>::max() - payloadOffset,
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte offset overflow while computing file offset.");
  return dataOffset + payloadOffset;
}

template <class T>
void binaryIOReadArrayFromRoot(const std::string& filename,
                               const unsigned long long dataOffset,
                               T* data,
                               const unsigned long long count,
                               const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount = binaryIOCheckedByteCount(count, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
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
      binaryIOBroadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
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
    binaryIOBroadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
  }
}

template <class T>
void binaryIOReadArrayCollective(const std::string& filename,
                                 const unsigned long long dataOffset,
                                 T* data,
                                 const unsigned long long count,
                                 const unsigned long long globalOffset,
                                 const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount  = binaryIOCheckedByteCount(count, sizeof(T));
  const unsigned long long byteOffset = binaryIOCheckedByteOffset(dataOffset, globalOffset, sizeof(T));
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
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
    binaryIOBroadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
  }
}

inline BinaryIOFileHeader binaryIOReadHeaderFromFile(const std::string& filename,
                                                     const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  BinaryIOFileHeader header;
  std::memset(&header, 0, sizeof(BinaryIOFileHeader));
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed for '" << filename << "'.");
    MPI_Status status;
    const int readErr = MPI_File_read_at_all(file, 0, &header, sizeof(BinaryIOFileHeader), MPI_BYTE, &status);
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
      in.read(reinterpret_cast<char*>(&header), sizeof(BinaryIOFileHeader));
      success = static_cast<bool>(in);
    }
  }
  binaryIOCheckFileOpen(success, filename, comm, "reading");
  Teuchos::broadcast(*comm, 0, static_cast<int>(sizeof(BinaryIOFileHeader)), reinterpret_cast<char*>(&header));
  return header;
}

template <class LocalOrdinal, class GlobalOrdinal>
void binaryIOValidateMapHeader(const BinaryIOFileHeader& header) {
  const char expectedMagic[8] = {'T', 'P', 'B', 'I', 'O', '0', '0', '1'};
  TEUCHOS_TEST_FOR_EXCEPTION(std::memcmp(header.magic, expectedMagic, sizeof(expectedMagic)) != 0,
                             std::runtime_error,
                             "Tpetra::BinaryIO: File does not have a recognized Tpetra binary I/O header. "
                                 << "Legacy Xpetra binary matrix files must be converted with "
                                 << "Xpetra::IO::ConvertLegacyBinaryToBinary before reading with Tpetra::BinaryIO.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.version != binaryIOFileVersion,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Unsupported file version " << header.version << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(header.byteOrderMarker != binaryIOByteOrderMarker,
                             std::runtime_error,
                             "Tpetra::BinaryIO: Unsupported byte order marker.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.objectKind != binaryIOMapObject,
                             std::runtime_error,
                             "Tpetra::BinaryIO: File object kind " << header.objectKind
                                                                   << " does not match expected kind " << binaryIOMapObject << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(header.localOrdinalSize != sizeof(LocalOrdinal) || header.localOrdinalFlags != binaryIOTypeFlags<LocalOrdinal>(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File local ordinal type does not match this map reader instantiation.");
  TEUCHOS_TEST_FOR_EXCEPTION(header.globalOrdinalSize != sizeof(GlobalOrdinal) || header.globalOrdinalFlags != binaryIOTypeFlags<GlobalOrdinal>(),
                             std::runtime_error,
                             "Tpetra::BinaryIO: File global ordinal type does not match this map reader instantiation.");
}

inline BinaryIOMapSectionHeader binaryIOReadMapSectionHeader(const std::string& filename,
                                                             const unsigned long long mapSectionOffset,
                                                             const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  BinaryIOMapSectionHeader header;
  std::memset(&header, 0, sizeof(BinaryIOMapSectionHeader));
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int err = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading map section header.");
    MPI_Status status;
    const int readErr = MPI_File_read_at_all(file, static_cast<MPI_Offset>(mapSectionOffset),
                                             &header, sizeof(BinaryIOMapSectionHeader), MPI_BYTE, &status);
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
    in.read(reinterpret_cast<char*>(&header), sizeof(BinaryIOMapSectionHeader));
    TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to read map section header from file '" << filename << "'.");
  }
  Teuchos::broadcast(*comm, 0, static_cast<int>(sizeof(BinaryIOMapSectionHeader)), reinterpret_cast<char*>(&header));
  return header;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
binaryIOReadMapSection(const std::string& filename,
                       const unsigned long long mapSectionOffset,
                       const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  using map_type                               = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  const BinaryIOMapSectionHeader sectionHeader = binaryIOReadMapSectionHeader(filename, mapSectionOffset, comm);
  const unsigned long long globalCount         = sectionHeader.numGlobalElements;

  TEUCHOS_TEST_FOR_EXCEPTION(sectionHeader.numRanks != static_cast<unsigned long long>(comm->getSize()),
                             std::runtime_error,
                             "Tpetra::BinaryIO: Map section was written for " << sectionHeader.numRanks
                                                                              << " ranks, but the read communicator has "
                                                                              << comm->getSize() << " ranks.");

  Kokkos::View<unsigned long long*, Kokkos::HostSpace> localCounts("Tpetra::BinaryIO::mapLocalCounts",
                                                                   binaryIOCheckedSize(sectionHeader.numRanks, "map section rank count"));
  binaryIOReadArrayFromRoot(filename,
                            mapSectionOffset + static_cast<unsigned long long>(sizeof(BinaryIOMapSectionHeader)),
                            localCounts.data(),
                            static_cast<unsigned long long>(localCounts.extent(0)),
                            comm);

  const unsigned long long localCount = localCounts(comm->getRank());
  unsigned long long globalOffset     = 0;
  for (int rank = 0; rank < comm->getRank(); ++rank) {
    globalOffset += localCounts(rank);
  }

  Kokkos::View<GlobalOrdinal*, Kokkos::HostSpace> gids("Tpetra::BinaryIO::mapGids",
                                                       binaryIOCheckedSize(localCount, "map section local element count"));
  if (localCount > 0) {
    binaryIOReadArrayCollective(filename,
                                mapSectionOffset + static_cast<unsigned long long>(sizeof(BinaryIOMapSectionHeader)) + sectionHeader.numRanks * sizeof(unsigned long long),
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

}  // namespace Details

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
readBinaryMapFile(const std::string& filename,
                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const Details::BinaryIOFileHeader header = Details::binaryIOReadHeaderFromFile(filename, comm);
  Details::binaryIOValidateMapHeader<LocalOrdinal, GlobalOrdinal>(header);
  return Details::binaryIOReadMapSection<LocalOrdinal, GlobalOrdinal, Node>(filename, header.rowMapOffset, comm);
}

}  // namespace Tpetra

//
// Explicit template instantiation macro for Tpetra::readBinaryMapFile.
// This function is templated only on LocalOrdinal, GlobalOrdinal, and Node
// (the LGN scheme), so it is instantiated separately from Tpetra::BinaryIO.
//
#define TPETRA_READBINARYMAPFILE_INSTANT(LO, GO, NODE)   \
  template Teuchos::RCP<const Tpetra::Map<LO, GO, NODE>> \
  readBinaryMapFile<LO, GO, NODE>(                       \
      const std::string&,                                \
      const Teuchos::RCP<const Teuchos::Comm<int>>&);

#endif  // TPETRA_READBINARYMAPFILE_DEF_HPP
