// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_mpiIsInitialized.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "mpi.h"
#  include <iostream>
#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {
  namespace Details {

    bool mpiIsInitialized ()
    {
#ifdef HAVE_TPETRACORE_MPI
      int isInitialized = 0;
      const int errCode = MPI_Initialized (&isInitialized);
      // If the call failed, then assume MPI wasn't implemented
      // correctly and return false.
      return errCode == MPI_SUCCESS && (isInitialized != 0);
#else
      return false; // Tpetra was not built with MPI support
#endif // HAVE_TPETRACORE_MPI
    }

    bool mpiIsFinalized ()
    {
#ifdef HAVE_TPETRACORE_MPI
      int isFinalized = 0;
      const int errCode = MPI_Finalized (&isFinalized);
      // If the call failed, then assume MPI wasn't implemented
      // correctly and return false.
      return errCode == MPI_SUCCESS && (isFinalized != 0);
#else
      return false; // Tpetra was not built with MPI support
#endif // HAVE_TPETRACORE_MPI
    }

  } // namespace Details
} // namespace Tpetra
