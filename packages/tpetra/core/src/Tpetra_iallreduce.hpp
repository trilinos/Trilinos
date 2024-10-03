// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_IALLREDUCE_HPP
#define TPETRA_IALLREDUCE_HPP

/// \file Tpetra_iallreduce.hpp
/// \brief Declaration of Tpetra::iallreduce
///
/// Tpetra::iallreduce wraps MPI_Iallreduce, with a nonblocking
/// fall-back implementation if MPI_Iallreduce doesn't exist (i.e., if
/// MPI_VERSION < 3).

#include "Tpetra_Details_iallreduce.hpp"

namespace Tpetra {

using ::Tpetra::Details::iallreduce;

} // namespace Tpetra

#endif // TPETRA_IALLREDUCE_HPP
