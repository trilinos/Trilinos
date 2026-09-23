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

#include <cstddef>
#include <limits>

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

}  // namespace Details
}  // namespace Tpetra

#endif  // TPETRA_BINARYIO_HELPERS_HPP
