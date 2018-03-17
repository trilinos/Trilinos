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

#include "Tpetra_Core.hpp"
#include "Tpetra_Details_mpiIsInitialized.hpp"

#ifdef HAVE_TPETRACORE_MPI
#  include <Teuchos_DefaultMpiComm.hpp> // this includes mpi.h too
#else
#  include <Teuchos_DefaultSerialComm.hpp>
#endif // HAVE_TPETRACORE_MPI

#include <Kokkos_Core.hpp>

namespace Tpetra {
  namespace { // (anonymous)

    class HideOutputExceptOnProcess0 {
    public:
      HideOutputExceptOnProcess0 (std::ostream& stream,
				  const int myRank) :
        stream_ (stream),
        originalBuffer_ (stream.rdbuf ())
      {
        if (myRank != 0) {
          stream.rdbuf (blackHole_.rdbuf ());
        }
      }

      ~HideOutputExceptOnProcess0 () {
        stream_.rdbuf (originalBuffer_);
      }
    private:
      std::ostream& stream_;
      decltype (std::cout.rdbuf ()) originalBuffer_;
      Teuchos::oblackholestream blackHole_;
    };

    // Whether one of the Tpetra::initialize() functions has been called before.
    bool tpetraIsInitialized_ = false;

    // Whether Tpetra initialized Kokkos.  Tpetra::finalize only
    // finalizes Kokkos if it initialized Kokkos.  Otherwise,
    // something else initialized Kokkos and is responsible for
    // finalizing it.
    bool tpetraInitializedKokkos_ = false;

#ifdef HAVE_TPETRACORE_MPI
    // Whether Tpetra initialized MPI.  Tpetra::finalize only
    // finalizes MPI if it initialized MPI.  Otherwise, something else
    // initialized MPI and is responsible for finalizing it.
    bool tpetraInitializedMpi_ = false;
#endif // HAVE_TPETRACORE_MPI

    // Tpetra's default communicator, wrapped in a Teuchos wrapper.
    // After Tpetra::finalize() is called, this GOES AWAY (is set to null).
    Teuchos::RCP<const Teuchos::Comm<int> > wrappedDefaultComm_;

    // Initialize Kokkos, if it needs initialization.
    // This takes the same arguments as (the first two of) initialize().
    void initKokkos (int* argc, char*** argv)
    {
      if (! tpetraInitializedKokkos_) {
        // Kokkos doesn't have a global is_initialized().  However,
        // Kokkos::initialize() always initializes the default execution
        // space, so it suffices to check whether that was initialized.
        const bool kokkosIsInitialized =
          Kokkos::DefaultExecutionSpace::is_initialized ();
        if (! kokkosIsInitialized) {
          int myRank = 0;
          const bool mpiHasBeenInitialized = Details::mpiIsInitialized ();
          if (mpiHasBeenInitialized) { // always false if not built with MPI
            auto comm = getDefaultComm ();
            myRank = comm->getRank ();
          }

          HideOutputExceptOnProcess0 hideCerr (std::cerr, myRank);
          HideOutputExceptOnProcess0 hideCout (std::cout, myRank);

          // Unlike MPI_Init, Kokkos promises not to modify argc and argv.
          Kokkos::initialize (*argc, *argv);
          tpetraInitializedKokkos_ = true;
        }
      }

      const bool kokkosIsInitialized =
        Kokkos::DefaultExecutionSpace::is_initialized ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (! kokkosIsInitialized, std::logic_error, "At the end of initKokkos, "
         "Kokkos is not initialized.  Please report this bug to the Tpetra "
         "developers.");
    }

#ifdef HAVE_TPETRACORE_MPI
    // Initialize MPI, if needed, and check for errors.  This takes
    // the same arguments as MPI_Init and the first two arguments of
    // initialize().
    void initMpi (int* argc, char*** argv)
    {
      const bool mpiAlreadyInitialized = Details::mpiIsInitialized ();
      if (! mpiAlreadyInitialized) {
	// Tpetra doesn't currently need to call MPI_Init_thread,
	// since with Tpetra, only one thread ever calls MPI
	// functions.  If we ever want to explore
	// MPI_THREAD_MULTIPLE, here would be the place to call
	// MPI_Init_thread.
	const int err = MPI_Init (argc, argv);
        TEUCHOS_TEST_FOR_EXCEPTION
          (err != MPI_SUCCESS, std::runtime_error, "MPI_Init failed with "
           "error code " << err << " != MPI_SUCCESS.  If MPI was set up "
           "correctly, then this should not happen, since we have already "
	   "checked that MPI_Init (or MPI_Init_thread) has not yet been "
	   "called.  This may indicate that your MPI library is corrupted "
	   "or that it is incorrectly linked to your program.");
        tpetraInitializedMpi_ = true;
      }
    }
#else
    void initMpi (int* /* argc */, char*** /* argv */) {}
#endif // HAVE_TPETRACORE_MPI

  } // namespace (anonymous)

  bool isInitialized () {
    return tpetraIsInitialized_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm ()
  {
    // It's OK to call this function if Tpetra::initialize hasn't been
    // called.  Users aren't obligated to call Tpetra::initialize, as
    // long as they take responsibility for initializing Kokkos and
    // MPI.  However, MPI must be initialized.  (Kokkos need not be
    // initialized for this method to be called.)

#ifdef HAVE_TPETRACORE_MPI
    const bool mpiInitd = Details::mpiIsInitialized ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (! mpiInitd, std::runtime_error,
       "Tpetra::getDefaultComm: MPI has not been initialized.  Before "
       "calling this method, you must either initialize MPI (by calling "
       "MPI_Init), or you must call Tpetra::initialize (which initializes "
       "MPI, if it has not yet been initialized).");
#endif // HAVE_TPETRACORE_MPI

    // This serves as lazy initialization.  Thus, the various
    // Tpetra::initialize functions do not need to set up
    // wrappedDefaultComm_.
    if (wrappedDefaultComm_.is_null ()) {
      Teuchos::RCP<const Teuchos::Comm<int> > comm;
#ifdef HAVE_TPETRACORE_MPI
      comm = Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
      comm = Teuchos::rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_TPETRACORE_MPI
      wrappedDefaultComm_ = comm;
    }
    return wrappedDefaultComm_;
  }

  void initialize (int* argc, char*** argv)
  {
    if (! tpetraIsInitialized_) {
      initMpi (argc, argv); // initialize MPI, if needed
      initKokkos (argc, argv); // initialize Kokkos, if needed
    }
    tpetraIsInitialized_ = true;
  }

#ifdef HAVE_TPETRACORE_MPI
  void initialize (int* argc, char*** argv, MPI_Comm comm)
  {
    initialize (argc, argv);
    // Set the default communicator.  We set it here, after the above
    // initialize() call, just in case users have not yet initialized
    // MPI.  (This is legal if users pass in a predefined
    // communicator, like MPI_COMM_WORLD or MPI_COMM_SELF.)
    //
    // What if users have already called initialize() before, but with
    // a different default communicator?  There are two possible
    // things we could do here:
    //
    //   1. Test via MPI_Comm_compare whether comm differs from the
    //      raw MPI communicator in wrappedDefaultComm_ (if indeed it
    //      is an MpiComm).
    //   2. Accept that the user might want to change the default
    //      communicator, and let them do it.
    //
    // I prefer #2.  Perhaps it would be sensible to print a warning
    // here, but on which process?  Would we use the old or the new
    // communicator to find that process' rank?  We don't want to use
    // MPI_COMM_WORLD's Process 0, since neither communicator might
    // include that process.  Furthermore, in some environments, only
    // Process 0 in MPI_COMM_WORLD is allowed to do I/O.  Thus, we
    // just let the change go without a warning.
    wrappedDefaultComm_ = Teuchos::rcp (new Teuchos::MpiComm<int> (comm));
  }
#endif // HAVE_TPETRACORE_MPI

  void
  initialize (int* argc, char*** argv,
              const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
    initialize (argc, argv);
    // See notes above on why we set the default communicator after
    // calling two-argument initialize().
    wrappedDefaultComm_ = comm;
  }

  void finalize ()
  {
    if (! tpetraIsInitialized_) {
      return; // user didn't call initialize(), so do nothing at all
    }

    // Tpetra should only finalize Kokkos if it initialized Kokkos.
    // See Github Issue #434.
    if (tpetraInitializedKokkos_) {
      Kokkos::finalize ();
    }

    // Make sure that no outstanding references to the communicator
    // remain.  If users gave initialize() an MPI_Comm, _they_ are
    // responsible for freeing it before calling finalize().
    wrappedDefaultComm_ = Teuchos::null;

#ifdef HAVE_TPETRACORE_MPI
    // Tpetra should only finalize MPI if it initialized MPI.
    // See Github Issue #434.
    if (tpetraInitializedMpi_) {
      // finalize() is a kind of destructor, so it's a bad idea to
      // throw an exception on error.  MPI implementations do have
      // the option to throw on error, so let's catch that here.
      try {
	if (Details::mpiIsInitialized ()) {
	  // This must be called by the same thread that called
	  // MPI_Init or MPI_Init_thread (possibly, but not
	  // necessarily, in Tpetra::initialize()).
	  (void) MPI_Finalize ();
	}
      }
      catch (...) {}
    }
#endif // HAVE_TPETRACORE_MPI

    tpetraIsInitialized_ = false; // it's not anymore.
  }
} // namespace Tpetra
