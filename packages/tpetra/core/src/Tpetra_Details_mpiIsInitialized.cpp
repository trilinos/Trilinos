// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
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
