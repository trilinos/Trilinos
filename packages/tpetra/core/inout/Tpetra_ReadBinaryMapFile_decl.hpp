// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_READBINARYMAPFILE_DECL_HPP
#define TPETRA_READBINARYMAPFILE_DECL_HPP

/// \file Tpetra_ReadBinaryMapFile_decl.hpp
/// \brief Declaration of the scalar-free Tpetra::readBinaryMapFile reader.
///
/// \c Tpetra::readBinaryMapFile reads a \c Tpetra::Map from a file written in
/// Tpetra's scalable binary I/O format (see \c Tpetra::BinaryIO).  It is a free
/// function templated only on <tt>LocalOrdinal, GlobalOrdinal, Node</tt>
/// (there is no Scalar), so it is explicitly instantiated on the LGN template
/// scheme.  Keeping it in its own translation unit lets its explicit template
/// instantiation stay independent of \c Tpetra::BinaryIO, which is instantiated
/// on the Scalar-LO-GO-Node scheme.

#include "Tpetra_Map.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include <string>

namespace Tpetra {

/// \brief Read a Tpetra::Map from a file in Tpetra's binary I/O format.
///
/// \param filename [in] Name of the binary file to read.
/// \param comm [in] Communicator over which to distribute the map.  It must
///   have the same number of ranks as were used when the file was written.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
readBinaryMapFile(const std::string& filename,
                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

}  // namespace Tpetra

#endif  // TPETRA_READBINARYMAPFILE_DECL_HPP
