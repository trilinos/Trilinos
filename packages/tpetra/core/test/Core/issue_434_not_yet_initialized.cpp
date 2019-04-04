/*
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
*/

// Test the case where Tpetra::initialize is called, when neither
// Kokkos nor MPI have been initialized yet.  Tpetra::initialize must,
// in that case, initialize both Kokkos and MPI.  Furthermore,
// Tpetra::finalize must finalize both Kokkos and MPI.

#include "Tpetra_Core.hpp"
#include "Kokkos_Core.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "mpi.h" // need direct MPI calls
#endif // HAVE_TPETRACORE_MPI
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <stdexcept>

namespace { // (anonymous)

  // Is Kokkos initialized?
  bool kokkosInitialized ()
  {
    // Kokkos lacks a global is_initialized function.  However,
    // Kokkos::initialize always initializes the default execution
    // space, so we can check that.
    return Kokkos::is_initialized ();
  }
}

int
main (int argc, char* argv[])
{
  using std::endl;
  bool success = true;
  std::ostream& out = std::cerr; // for unbuffered output, unlike std::cout
  std::ostream& err = std::cerr;

  out << "About to call Tpetra::initialize" << endl;
  try {
    Tpetra::initialize (&argc, &argv);
  }
  catch (std::exception& e) {
    err << "*** Tpetra::getDefaultComm threw an std::exception: "
        << e.what () << endl;
    success = false;
    goto EndOfTest;
  }
  catch (...) {
    err << "*** Tpetra::initialize threw an exception not a subclass of "
      "std::exception" << endl;
    success = false;
    goto EndOfTest;
  }

  out << "Tpetra::initialize did not throw an exception" << endl;

  // Neither Kokkos nor MPI have been initialized yet.
  // Thus, Tpetra::initialize must initialize both of them.

#ifdef HAVE_TPETRACORE_MPI
  {
    int mpiIsInitialized = 0;
    const int mpiErr = MPI_Initialized (&mpiIsInitialized);

    if (mpiErr != MPI_SUCCESS) {
      err << "*** MPI_Initialized returned err != MPI_SUCCESS" << endl;
      success = false;
      // If MPI_Initialized fails, MPI is totally messed up.  The best
      // we can do at this point is stop the program.
      goto EndOfTest;
    }
    else {
      if (mpiIsInitialized == 0) {
        err << "*** MPI_Initialized says MPI was not initialized" << endl;
        success = false;
        // If MPI is not initialized, it's safe just to stop the program.
        goto EndOfTest;
      }
      else {
        out << "MPI is initialized" << endl;
      }
    }
  }
#endif // HAVE_TPETRACORE_MPI

  if (! kokkosInitialized ()) {
    err << "*** Kokkos is not initialized" << endl;
    success = false;
  }
  else {
    out << "Kokkos is initialized" << endl;
  }

  {
    // Must hide this from the label, else the compiler reports an error.
    Teuchos::RCP<const Teuchos::Comm<int> > comm;
    bool threw = false;
    try {
      comm = Tpetra::getDefaultComm ();
    }
    catch (std::exception& e) {
      err << "*** Tpetra::getDefaultComm throw an std::exception: "
          << e.what () << endl;
      threw = true;
      success = false;
    }
    catch (...) {
      err << "*** Tpetra::getDefaultComm throw an exception not a subclass of "
        "std::exception" << endl;
      threw = true;
      success = false;
    }

    if (! threw && comm.is_null ()) {
      err << "*** Tpetra::getDefaultComm returned null (but did not throw an "
        "exception) after Tpetra was initialized" << endl;
      success = false;
    }
    else {
      out << "Tpetra::getDefaultComm returned nonnull" << endl;
    }
  }

  {
    // Must hide this from the label, else the compiler reports an error.
    bool tpetraFinalizeThrew = false;
    try {
      Tpetra::finalize ();
    }
    catch (std::exception& e) {
      err << "*** Tpetra::finalize threw an std::exception: "
          << e.what () << endl;
      tpetraFinalizeThrew = true;
      success = false;
    }
    catch (...) {
      err << "*** Tpetra::finalize threw an exception not a subclass of "
        "std::exception" << endl;
      tpetraFinalizeThrew = true;
      success = false;
    }

    if (! tpetraFinalizeThrew) {
      out << "Tpetra::finalize did not throw an exception" << endl;
    }
  }

#ifdef HAVE_TPETRACORE_MPI
  {
    int mpiIsFinalized = 0;
    int mpiErr = MPI_Finalized (&mpiIsFinalized);

    if (mpiErr != MPI_SUCCESS) {
      err << "*** MPI_Finalized returned err != MPI_SUCCESS" << endl;
      success = false;
      // If MPI_Finalized fails, MPI is totally messed up.  The best
      // we can do at this point is stop the program.
      goto EndOfTest;
    }
    else {
      if (mpiIsFinalized == 0) {
        err << "*** MPI_Finalized says MPI was not finalized" << endl;
        success = false;
        // If MPI was initialized (it must be if we got this far), but
        // was not finalized, then we need to finalize MPI.  After
        // that, we can stop the program.  Don't try to check the
        // error code; what would we do at this point if it fails,
        // other than stop the program?
        (void) MPI_Finalize ();
        goto EndOfTest;
      }
    }
  }
#endif // HAVE_TPETRACORE_MPI

  // Kokkos doesn't have an "is_finalized" function.  Some execution
  // spaces (like OpenMP) don't ever shut down (though Kokkos might
  // change the number of OpenMP threads to 1).  Thus, we don't try to
  // test for Kokkos finalization here.

 EndOfTest:

  if (success) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}

