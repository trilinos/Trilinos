// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_printOnce.hpp"

#if defined(HAVE_TPETRACORE_MPI)
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_Comm.hpp"
#endif // defined(HAVE_TPETRACORE_MPI)

namespace { // (anonymous)
  bool mpiIsInitialized ()
  {
#if defined(HAVE_TPETRACORE_MPI)
    int isInitialized = 0;
    try {
      (void) MPI_Initialized (&isInitialized);
    }
    catch (...) {
      // Not sure if MPI_Initialized meets strong exception guarantee
      isInitialized = 0;
    }
    return isInitialized != 0;
#else
    return false;
#endif // defined(HAVE_TPETRACORE_MPI)
  }

  bool mpiIsFinalized ()
  {
#if defined(HAVE_TPETRACORE_MPI)
    int isFinalized = 0;
    try {
      (void) MPI_Finalized (&isFinalized);
    }
    catch (...) {
      // Not sure if MPI_Initialized meets strong exception guarantee
      isFinalized = 0;
    }
    return isFinalized != 0;
#else
    return false;
#endif // defined(HAVE_TPETRACORE_MPI)
  }

#if defined(HAVE_TPETRACORE_MPI)  
  bool isMpiComm (const Teuchos::Comm<int>& comm)
  {
    using mpi_comm_type = Teuchos::MpiComm<int>;
    return dynamic_cast<const mpi_comm_type* > (&comm) != nullptr;
  }
#else
  bool isMpiComm (const Teuchos::Comm<int>& /* comm */ )
  {
    return false;
  }
#endif // defined(HAVE_TPETRACORE_MPI)    
  
  int getRankHarmlessly (const Teuchos::Comm<int>& comm)
  {
    if (mpiIsInitialized () && ! mpiIsFinalized () && isMpiComm (comm)) {
      return comm.getRank ();
    }
    else {
      return 0;
    }
  }
} // namespace (anonymous)  

namespace Tpetra {
  namespace Details {
    void
    printOnce (std::ostream& out,
	       const std::string& s,
	       const Teuchos::Comm<int>* comm)
    {
      if (comm == nullptr || getRankHarmlessly (*comm) == 0) {
      	out << s;
      }
    }
  } // namespace Details
} // namespace Tpetra
