// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_MPIISINITIALIZED
#define TPETRA_DETAILS_MPIISINITIALIZED

#include "TpetraCore_config.h"

namespace Tpetra {
  namespace Details {

    /// \brief Has MPI_Init been called (on this process)?
    ///
    /// If Tpetra was built with MPI support, then this wraps
    /// MPI_Initialized.  If Tpetra was not built with MPI support,
    /// then this always returns false, regardless of whether the user
    /// has built with MPI.
    ///
    /// MPI (at least 3.0) only permits MPI to be initialized once.
    /// After MPI_Init has been called on a process, MPI_Initialized
    /// always returns true on that process, regardless of whether
    /// MPI_Finalize has been called.
    ///
    /// If you want to know whether MPI_Finalize has been called on
    /// this process, use mpiIsFinalized() (see below).
    bool mpiIsInitialized ();

    /// \brief Has MPI_Finalize been called (on this process)?
    ///
    /// If Tpetra was built with MPI support, then this wraps
    /// MPI_Finalized.  If Tpetra was not built with MPI support, then
    /// this always returns false, regardless of whether the user has
    /// built with MPI.
    ///
    /// MPI (at least 3.0) only permits MPI_Init to be called at most
    /// once on a process.  After MPI_Finalize has been called
    /// successfully on a process, MPI_Finalized always returns true
    /// on that process.
    ///
    /// If you want to know whether MPI_Init has been called on this
    /// process, use mpiIsInitialized() (see above).
    bool mpiIsFinalized ();

  } // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MPIISINITIALIZED
