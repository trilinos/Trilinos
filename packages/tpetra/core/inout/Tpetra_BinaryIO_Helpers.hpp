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
/// \brief Small inline utilities shared by Tpetra's binary I/O readers and
///   writers.
///
/// These helpers are intentionally header-only and free of explicit template
/// instantiation.  They are shared between the \c Tpetra::BinaryIO class
/// (see Tpetra_BinaryIO_decl.hpp) and the free function
/// \c Tpetra::readBinaryMapFile (see Tpetra_ReadBinaryMapFile_decl.hpp), which
/// live in separate translation units so that their explicit template
/// instantiations stay independent.

#include "Teuchos_TestForException.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>

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

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_BINARYIO_HELPERS_HPP
