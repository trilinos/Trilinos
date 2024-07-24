// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_DefaultMpiComm.hpp>

// Only enable the contents of this file if building with MPI.
#ifdef HAVE_TEUCHOS_MPI

namespace Teuchos {

  std::string
  mpiErrorCodeToString (const int errCode)
  {
    if (errCode == MPI_SUCCESS) {
      return "MPI_SUCCESS";
    }
    else {
      char rawErrString[MPI_MAX_ERROR_STRING];
      int len = 0;
      int err = MPI_Error_string (errCode, rawErrString, &len);
      if (err != MPI_SUCCESS) {
        // Assume that the string wasn't written.  This means it might
        // not be null terminated, so make it a valid empty string by
        // writing the null termination character to it.
        if (MPI_MAX_ERROR_STRING > 0) {
          rawErrString[0] = '\0';
        }
      }
      return std::string (rawErrString);
    }
  }

  namespace details {
    void safeCommFree (MPI_Comm* comm) {
      // FIXME (mfh 08 Dec 2014) Use the MPI_Finalize hook trick to
      // call MPI_Comm_free at MPI_Finalize automatically, if it
      // hasn't already been called on the object.  Store the MPI_Comm
      // (by allocated pointer) as the value of the (key,value) pair
      // (used in the hook), and be sure to free the pair if the free
      // function is called before MPI_Finalize.
      int finalized = 0;
      const int err = MPI_Finalized (&finalized);
      // Just to be safe, don't do anything if calling MPI_Finalized
      // didn't succeed.  It's better to leak memory than to crash.
      if (err == MPI_SUCCESS && ! finalized) {
        // Don't throw an exception if MPI_Comm_free reports an error,
        // since we're likely to be in a destructor and destructors
        // shouldn't throw exceptions.
        (void) MPI_Comm_free (comm);
      }
    }

    int setCommErrhandler (MPI_Comm comm, MPI_Errhandler handler) {
#if MPI_VERSION >= 2
      return MPI_Comm_set_errhandler (comm, handler);
#else // MPI 1
      return MPI_Errhandler_set (comm, handler);
#endif // MPI_VERSION >= 2
    }
  } // namespace details

} // namespace Teuchos

#endif // HAVE_TEUCHOS_MPI
