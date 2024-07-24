// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_ISINTERCOMM_HPP
#define TPETRA_DETAILS_ISINTERCOMM_HPP

/// \file Tpetra_Details_isInterComm.hpp
/// \brief Declaration of Tpetra::Details::isInterComm
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.
///
/// Tpetra::Details::isInterComm wraps MPI_Comm_test_inter.  The
/// Tpetra::Details::isInterComm function is the only thing in this
/// file upon which Tpetra developers should rely.  Tpetra developers
/// should not rely on anything else in this file.  <i>Users</i> may
/// not rely on <i>anything</i> in this file!

#include "TpetraCore_config.h"
#include "Teuchos_Comm.hpp"

namespace Tpetra {
namespace Details {

/// \brief Return true if and only if the input communicator wraps an
///   MPI intercommunicator.
///
/// The most common MPI communicators are intracommunicators
/// ("in<i>tra</i>," not "in<i>ter</i>").  This includes
/// MPI_COMM_WORLD, MPI_COMM_SELF, and the results of MPI_Comm_dup and
/// MPI_Comm_split.  Intercommunicators come from
/// MPI_Intercomm_create.
///
/// This distinction matters because only collectives over
/// intracommunicators may use MPI_IN_PLACE, to let the send and
/// receive buffers alias each other.  Collectives over
/// intercommunicators may <i>not</i> use MPI_IN_PLACE.
///
/// \param comm [in] The input communicator.
///
/// \return Whether the input communicator wraps an MPI
///   intercommunicator.
bool
isInterComm (const Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ISINTERCOMM_HPP
