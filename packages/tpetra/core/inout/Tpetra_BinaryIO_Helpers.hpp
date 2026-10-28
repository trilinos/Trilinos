// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BINARYIO_HELPERS_HPP
#define TPETRA_BINARYIO_HELPERS_HPP

/// \file Tpetra_BinaryIO_Helpers.hpp
/// \brief Inline utilities shared by Tpetra's binary I/O readers and writers.
///
/// These helpers are intentionally header-only and free of explicit template
/// instantiation.  They are shared between the \c Tpetra::BinaryIO class
/// (see Tpetra_BinaryIO_decl.hpp) and the free function
/// \c Tpetra::readBinaryMapFile (see Tpetra_ReadBinaryMapFile_decl.hpp), which
/// live in separate translation units so that their explicit template
/// instantiations stay independent.

#include "TpetraCore_config.h"
#include "Tpetra_Details_extractMpiCommFromTeuchos.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace Tpetra {
namespace Details {

/// \brief Cast an unsigned 64-bit count to \c size_t, throwing if it does not fit.
inline size_t binaryIOCheckedSize(const unsigned long long count, const char label[]) {
  TEUCHOS_TEST_FOR_EXCEPTION(count > static_cast<unsigned long long>(std::numeric_limits<size_t>::max()),
                             std::overflow_error,
                             "Tpetra::BinaryIO: " << label << " does not fit in size_t.");
  return static_cast<size_t>(count);
}

/// \brief Cast a \c size_t count to \c LocalOrdinal, throwing if it does not fit.
template <class LocalOrdinal>
LocalOrdinal binaryIOCheckedLocalOrdinalCount(const size_t count, const char label[]) {
  TEUCHOS_TEST_FOR_EXCEPTION(count > static_cast<size_t>(std::numeric_limits<LocalOrdinal>::max()),
                             std::overflow_error,
                             "Tpetra::BinaryIO: " << label << " does not fit in LocalOrdinal.");
  return static_cast<LocalOrdinal>(count);
}

/// \brief Override byte chunks for unit tests; zero means use production limits.
inline unsigned long long& binaryIOMaxChunkOverrideForUnitTests() {
  static unsigned long long maxChunk = 0;
  return maxChunk;
}

/// \brief Set a small int-count chunk limit for BinaryIO unit tests.
inline void binaryIOSetMaxChunkForUnitTests(const unsigned long long maxChunk) {
  TEUCHOS_TEST_FOR_EXCEPTION(maxChunk > static_cast<unsigned long long>(std::numeric_limits<int>::max()),
                             std::logic_error,
                             "Tpetra::BinaryIO: Unit-test chunk limit must fit in int.");
  binaryIOMaxChunkOverrideForUnitTests() = maxChunk;
}

/// \brief Return the int-count chunk limit for byte transfers.
inline unsigned long long binaryIOIntCountMaxChunk() {
  const unsigned long long maxChunk = binaryIOMaxChunkOverrideForUnitTests();
  return maxChunk == 0 ? static_cast<unsigned long long>(std::numeric_limits<int>::max()) : maxChunk;
}

/// \brief Return the number of byte chunks needed for a transfer.
inline unsigned long long binaryIOChunkCount(const unsigned long long byteCount,
                                             const unsigned long long maxChunk) {
  TEUCHOS_TEST_FOR_EXCEPTION(maxChunk == 0,
                             std::logic_error,
                             "Tpetra::BinaryIO: Chunk size must be nonzero.");
  return byteCount == 0 ? 0 : 1ull + (byteCount - 1ull) / maxChunk;
}

/// \brief Return the size of one byte chunk in a chunked transfer.
inline unsigned long long binaryIOChunkSize(const unsigned long long byteCount,
                                            const unsigned long long maxChunk,
                                            const unsigned long long chunkIndex) {
  const unsigned long long numChunks = binaryIOChunkCount(byteCount, maxChunk);
  TEUCHOS_TEST_FOR_EXCEPTION(chunkIndex >= numChunks,
                             std::logic_error,
                             "Tpetra::BinaryIO: Chunk index is outside the chunked transfer.");
  const unsigned long long chunkOffset = chunkIndex * maxChunk;
  return std::min(byteCount - chunkOffset, maxChunk);
}

/// \brief Add a byte increment to a file offset, throwing on overflow.
inline unsigned long long binaryIOCheckedAddByteOffset(const unsigned long long offset,
                                                       const unsigned long long increment) {
  TEUCHOS_TEST_FOR_EXCEPTION(offset > std::numeric_limits<unsigned long long>::max() - increment,
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte offset overflow while advancing file offset.");
  return offset + increment;
}

/// \brief Reverse the byte order of a 64-bit value.
inline unsigned long long binaryIOByteSwap64(const unsigned long long value) {
  unsigned long long swapped = 0;
  for (int i = 0; i < 8; ++i) {
    swapped = (swapped << 8) | ((value >> (8 * i)) & 0xffull);
  }
  return swapped;
}

/// \brief Return the byte count of \c count elements of size \c elementSize, throwing on overflow.
inline unsigned long long binaryIOCheckedByteCount(const unsigned long long count,
                                                   const size_t elementSize) {
  TEUCHOS_TEST_FOR_EXCEPTION(count > 0 &&
                                 elementSize > static_cast<size_t>(std::numeric_limits<unsigned long long>::max() / count),
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte count overflow while computing transfer size.");
  return count * static_cast<unsigned long long>(elementSize);
}

/// \brief Return the file offset of element \c globalOffset of an array starting at \c dataOffset, throwing on overflow.
inline unsigned long long binaryIOCheckedByteOffset(const unsigned long long dataOffset,
                                                    const unsigned long long globalOffset,
                                                    const size_t elementSize) {
  const unsigned long long payloadOffset = binaryIOCheckedByteCount(globalOffset, elementSize);
  TEUCHOS_TEST_FOR_EXCEPTION(dataOffset > std::numeric_limits<unsigned long long>::max() - payloadOffset,
                             std::overflow_error,
                             "Tpetra::BinaryIO: Byte offset overflow while computing file offset.");
  return dataOffset + payloadOffset;
}

/// \brief Broadcast \c byteCount bytes from rank 0 in int-count chunks.
inline void binaryIOBroadcastBytesFromRoot(char data[],
                                           const unsigned long long byteCount,
                                           const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  unsigned long long remaining      = byteCount;
  char* current                     = data;
  const unsigned long long maxChunk = binaryIOIntCountMaxChunk();
  while (remaining > 0) {
    const unsigned long long chunk = std::min(remaining, maxChunk);
    Teuchos::broadcast(*comm, 0, static_cast<int>(chunk), current);
    current += static_cast<std::ptrdiff_t>(chunk);
    remaining -= chunk;
  }
}

/// \brief Throw unless the (non-MPI) communicator has exactly one process.
///
/// The non-MPI fallbacks of the collective array helpers access the file with
/// std::fstream on the calling process only, which is correct only when that
/// process is the entire communicator (e.g. \c Teuchos::SerialComm).
inline void binaryIORequireSingleProcessComm(const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() != 1,
                             std::logic_error,
                             "Tpetra::BinaryIO: Collective I/O on a non-MPI communicator requires exactly one process, but the "
                             "communicator has "
                                 << comm->getSize() << ".");
}

#ifdef HAVE_TPETRACORE_MPI
/// \brief Return true on every process if \c localFailed is true on any process.
inline bool binaryIOAnyProcessFailed(const bool localFailed,
                                     const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const int localFlag = localFailed ? 1 : 0;
  int globalFlag      = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localFlag, &globalFlag);
  return globalFlag != 0;
}

/// \brief Return true if the byte range [offset, offset + byteCount) lies inside the file.
///
/// MPI reads past end-of-file succeed with a short count, and some
/// implementations (e.g. ROMIO collective reads) report the requested count
/// regardless, so truncated files are detected by checking the file size.
inline bool binaryIOByteRangeInFile(MPI_File file,
                                    const unsigned long long offset,
                                    const unsigned long long byteCount) {
  if (byteCount == 0) {
    return true;
  }
  MPI_Offset fileSize = 0;
  if (MPI_File_get_size(file, &fileSize) != MPI_SUCCESS || fileSize < 0) {
    return false;
  }
  return offset <= static_cast<unsigned long long>(fileSize) &&
         byteCount <= static_cast<unsigned long long>(fileSize) - offset;
}
#endif

/// \brief Write \c count elements from rank 0 to \c filename at byte offset \c dataOffset.
template <class T>
void binaryIOWriteArrayFromRoot(const std::string& filename,
                                const unsigned long long dataOffset,
                                const T* data,
                                const unsigned long long count,
                                const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount = binaryIOCheckedByteCount(count, sizeof(T));
  if (byteCount > 0) {
    binaryIOCheckedAddByteOffset(dataOffset, byteCount - 1ull);
  }
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while writing root-owned array data.");
    int writeErr = MPI_SUCCESS;
    if (comm->getRank() == 0) {
      const char* current                = reinterpret_cast<const char*>(data);
      unsigned long long currentOffset   = dataOffset;
      const unsigned long long maxChunk  = binaryIOIntCountMaxChunk();
      const unsigned long long numChunks = binaryIOChunkCount(byteCount, maxChunk);
      for (unsigned long long chunkIndex = 0; chunkIndex < numChunks && writeErr == MPI_SUCCESS; ++chunkIndex) {
        const unsigned long long chunk = binaryIOChunkSize(byteCount, maxChunk, chunkIndex);
        MPI_Status status;
        writeErr = MPI_File_write_at(file,
                                     static_cast<MPI_Offset>(currentOffset),
                                     const_cast<char*>(current),
                                     static_cast<int>(chunk),
                                     MPI_BYTE,
                                     &status);
        current += static_cast<std::ptrdiff_t>(chunk);
        currentOffset += chunk;
      }
    }
    MPI_File_close(&file);
    Teuchos::broadcast(*comm, 0, 1, &writeErr);
    TEUCHOS_TEST_FOR_EXCEPTION(writeErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_write_at failed while writing root-owned array data.");
    comm->barrier();
    return;
  }
#endif
  if (comm->getRank() == 0) {
    std::fstream out(filename.c_str(), std::ios::binary | std::ios::in | std::ios::out);
    TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to open file '" << filename << "' while writing root-owned array data.");
    out.seekp(static_cast<std::streamoff>(dataOffset), std::ios::beg);

    const char* current               = reinterpret_cast<const char*>(data);
    unsigned long long remaining      = byteCount;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<std::streamsize>::max());
    while (remaining > 0) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      out.write(current, static_cast<std::streamsize>(chunk));
      TEUCHOS_TEST_FOR_EXCEPTION(!out.good(), std::runtime_error,
                                 "Tpetra::BinaryIO: Failed to write root-owned array data to file '" << filename << "'.");
      current += static_cast<std::ptrdiff_t>(chunk);
      remaining -= chunk;
    }
  }
  comm->barrier();
}

/// \brief Read \c count elements at byte offset \c dataOffset on rank 0 and broadcast them to all processes.
template <class T>
void binaryIOReadArrayFromRoot(const std::string& filename,
                               const unsigned long long dataOffset,
                               T* data,
                               const unsigned long long count,
                               const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount = binaryIOCheckedByteCount(count, sizeof(T));
  if (byteCount > 0) {
    binaryIOCheckedAddByteOffset(dataOffset, byteCount - 1ull);
  }
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading root-owned array data.");
    int readErr = MPI_SUCCESS;
    if (comm->getRank() == 0) {
      if (!binaryIOByteRangeInFile(file, dataOffset, byteCount)) {
        readErr = MPI_ERR_IO;
      }
      char* current                      = reinterpret_cast<char*>(data);
      unsigned long long currentOffset   = dataOffset;
      const unsigned long long maxChunk  = binaryIOIntCountMaxChunk();
      const unsigned long long numChunks = binaryIOChunkCount(byteCount, maxChunk);
      for (unsigned long long chunkIndex = 0; chunkIndex < numChunks && readErr == MPI_SUCCESS; ++chunkIndex) {
        const unsigned long long chunk = binaryIOChunkSize(byteCount, maxChunk, chunkIndex);
        MPI_Status status;
        readErr = MPI_File_read_at(file,
                                   static_cast<MPI_Offset>(currentOffset),
                                   current,
                                   static_cast<int>(chunk),
                                   MPI_BYTE,
                                   &status);
        current += static_cast<std::ptrdiff_t>(chunk);
        currentOffset += chunk;
      }
    }
    MPI_File_close(&file);
    Teuchos::broadcast(*comm, 0, 1, &readErr);
    TEUCHOS_TEST_FOR_EXCEPTION(readErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: Failed to read root-owned array data from file '"
                                   << filename << "' (file truncated or MPI_File_read_at failed).");
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

    char* current                     = reinterpret_cast<char*>(data);
    unsigned long long remaining      = byteCount;
    const unsigned long long maxChunk = static_cast<unsigned long long>(std::numeric_limits<std::streamsize>::max());
    while (remaining > 0) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      in.read(current, static_cast<std::streamsize>(chunk));
      TEUCHOS_TEST_FOR_EXCEPTION(!in.good(), std::runtime_error,
                                 "Tpetra::BinaryIO: Failed to read root-owned array data from file '" << filename << "'.");
      current += static_cast<std::ptrdiff_t>(chunk);
      remaining -= chunk;
    }
  }
  if (comm->getSize() > 1) {
    binaryIOBroadcastBytesFromRoot(reinterpret_cast<char*>(data), byteCount, comm);
  }
}

/// \brief Collectively write each process's \c count elements at element offset \c globalOffset of the array at \c dataOffset.
///
/// Every process must call this, including processes with \c count == 0.
template <class T>
void binaryIOWriteArrayCollective(const std::string& filename,
                                  const unsigned long long dataOffset,
                                  const T* data,
                                  const unsigned long long count,
                                  const unsigned long long globalOffset,
                                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount  = binaryIOCheckedByteCount(count, sizeof(T));
  const unsigned long long byteOffset = binaryIOCheckedByteOffset(dataOffset, globalOffset, sizeof(T));
  binaryIOCheckedAddByteOffset(byteOffset, byteCount);
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while writing array data.");

    const char* current                     = reinterpret_cast<const char*>(data);
    unsigned long long remaining            = byteCount;
    unsigned long long currentOffset        = byteOffset;
    const unsigned long long maxChunk       = binaryIOIntCountMaxChunk();
    const unsigned long long localNumChunks = binaryIOChunkCount(byteCount, maxChunk);
    unsigned long long numChunks            = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localNumChunks, &numChunks);
    // Every process makes every collective call, even after a local error or
    // once its own data is exhausted; errors are combined after the loop.
    int writeErr = MPI_SUCCESS;
    for (unsigned long long chunkIndex = 0; chunkIndex < numChunks; ++chunkIndex) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      MPI_Status status;
      const int err = MPI_File_write_at_all(file,
                                            static_cast<MPI_Offset>(currentOffset),
                                            const_cast<char*>(current),
                                            static_cast<int>(chunk),
                                            MPI_BYTE,
                                            &status);
      if (writeErr == MPI_SUCCESS) {
        writeErr = err;
      }
      current += static_cast<std::ptrdiff_t>(chunk);
      currentOffset += chunk;
      remaining -= chunk;
    }
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(binaryIOAnyProcessFailed(writeErr != MPI_SUCCESS, comm), std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_write_at_all failed on at least one process while writing array data to file '"
                                   << filename << "'.");
    return;
  }
#endif
  binaryIORequireSingleProcessComm(comm);
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

/// \brief Collectively read each process's \c count elements at element offset \c globalOffset of the array at \c dataOffset.
///
/// Every process must call this, including processes with \c count == 0.
template <class T>
void binaryIOReadArrayCollective(const std::string& filename,
                                 const unsigned long long dataOffset,
                                 T* data,
                                 const unsigned long long count,
                                 const unsigned long long globalOffset,
                                 const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  const unsigned long long byteCount  = binaryIOCheckedByteCount(count, sizeof(T));
  const unsigned long long byteOffset = binaryIOCheckedByteOffset(dataOffset, globalOffset, sizeof(T));
  binaryIOCheckedAddByteOffset(byteOffset, byteCount);
#ifdef HAVE_TPETRACORE_MPI
  if (teuchosCommIsAnMpiComm(*comm)) {
    MPI_Comm rawComm = extractMpiCommFromTeuchos(*comm);
    MPI_File file;
    const int openErr = MPI_File_open(rawComm, const_cast<char*>(filename.c_str()), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    TEUCHOS_TEST_FOR_EXCEPTION(openErr != MPI_SUCCESS, std::runtime_error,
                               "Tpetra::BinaryIO: MPI_File_open failed while reading array data.");

    const bool inFile                       = binaryIOByteRangeInFile(file, byteOffset, byteCount);
    char* current                           = reinterpret_cast<char*>(data);
    unsigned long long remaining            = byteCount;
    unsigned long long currentOffset        = byteOffset;
    const unsigned long long maxChunk       = binaryIOIntCountMaxChunk();
    const unsigned long long localNumChunks = binaryIOChunkCount(byteCount, maxChunk);
    unsigned long long numChunks            = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &localNumChunks, &numChunks);
    // Every process makes every collective call, even after a local error or
    // once its own data is exhausted; errors are combined after the loop.
    int readErr = MPI_SUCCESS;
    for (unsigned long long chunkIndex = 0; chunkIndex < numChunks; ++chunkIndex) {
      const unsigned long long chunk = std::min(remaining, maxChunk);
      MPI_Status status;
      const int err = MPI_File_read_at_all(file,
                                           static_cast<MPI_Offset>(currentOffset),
                                           current,
                                           static_cast<int>(chunk),
                                           MPI_BYTE,
                                           &status);
      if (readErr == MPI_SUCCESS) {
        readErr = err;
      }
      current += static_cast<std::ptrdiff_t>(chunk);
      currentOffset += chunk;
      remaining -= chunk;
    }
    MPI_File_close(&file);
    TEUCHOS_TEST_FOR_EXCEPTION(binaryIOAnyProcessFailed(!inFile || readErr != MPI_SUCCESS, comm), std::runtime_error,
                               "Tpetra::BinaryIO: Failed to read array data from file '"
                                   << filename << "' on at least one process (file truncated or MPI_File_read_at_all failed).");
    return;
  }
#endif
  binaryIORequireSingleProcessComm(comm);
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

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_BINARYIO_HELPERS_HPP
