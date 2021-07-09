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
// ************************************************************************
// @HEADER

#include "Tpetra_Core.hpp"
#include "Tpetra_Details_mpiIsInitialized.hpp"

#ifdef HAVE_TPETRACORE_MPI
#  include <Teuchos_DefaultMpiComm.hpp> // this includes mpi.h too
#endif // HAVE_TPETRACORE_MPI
#include <Teuchos_DefaultSerialComm.hpp>

#include <Kokkos_Core.hpp>
#include "Tpetra_Details_checkLaunchBlocking.hpp"

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

#if defined(HAVE_TPETRACORE_MPI)
    bool mpiIsInitializedAndNotFinalized ()
    {
      int isInitialized = 0;
      int isFinalized = 0;
      // Not sure if MPI_Initialized or MPI_Finalized meet the strong
      // exception guarantee.
      try {
        (void) MPI_Initialized (&isInitialized);
      }
      catch (...) {
        isInitialized = 0;
      }
      try {
        (void) MPI_Finalized (&isFinalized);
      }
      catch (...) {
        isFinalized = 0;
      }
      return isInitialized != 0 && isFinalized == 0;
    }

    int getRankHarmlessly (MPI_Comm comm)
    {
      int myRank = 0;
      if (mpiIsInitializedAndNotFinalized ()) {
        try {
          (void) MPI_Comm_rank (comm, &myRank);
        }
        catch (...) {
          // Not sure if MPI_Comm_rank meets strong exception guarantee
          myRank = 0;
        }
      }
      return myRank;
    }
#endif // defined(HAVE_TPETRACORE_MPI)


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

    // This takes the same arguments as (the first two of) initialize().
    void initKokkosIfNeeded (int* argc, char*** argv, const int myRank)
    {
      if (! tpetraInitializedKokkos_) {
        // Kokkos doesn't have a global is_initialized().  However,
        // Kokkos::initialize() always initializes the default execution
        // space, so it suffices to check whether that was initialized.
        const bool kokkosIsInitialized =
          Kokkos::is_initialized ();
        if (! kokkosIsInitialized) {
          HideOutputExceptOnProcess0 hideCerr (std::cerr, myRank);
          HideOutputExceptOnProcess0 hideCout (std::cout, myRank);

          // Unlike MPI_Init, Kokkos promises not to modify argc and argv.
          Kokkos::initialize (*argc, *argv);
          tpetraInitializedKokkos_ = true;
        }
      }
      Details::checkOldCudaLaunchBlocking();
      const bool kokkosIsInitialized =
        Kokkos::is_initialized ();
      TEUCHOS_TEST_FOR_EXCEPTION
        (! kokkosIsInitialized, std::logic_error, "At the end of "
         "initKokkosIfNeeded, Kokkos is not initialized.  "
         "Please report this bug to the Tpetra developers.");
    }

#ifdef HAVE_TPETRACORE_MPI
    // This takes the same arguments as MPI_Init and the first two
    // arguments of initialize().
    void initMpiIfNeeded (int* argc, char*** argv)
    {
      // Both MPI_Initialized and MPI_Finalized report true after
      // MPI_Finalize has been called.  It's not legal to call
      // MPI_Init after MPI_Finalize has been called (see MPI 3.0
      // Standard, Section 8.7).  It would be unusual for users to
      // want to use Tpetra after MPI_Finalize has been called, but
      // there's no reason why we should forbid it.  It just means
      // that Tpetra will need to run without MPI.

      const bool mpiReady = mpiIsInitializedAndNotFinalized ();
      if (! mpiReady) {
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
#endif // HAVE_TPETRACORE_MPI

  } // namespace (anonymous)

  bool isInitialized () {
    return tpetraIsInitialized_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm ()
  {
    // It's technically not correct to call this function if Tpetra
    // has not yet been initialized, but requiring that may break some
    // downstream tests.
    //
    // This function initializes wrappedDefaultComm_ lazily.
    // Tpetra::initialize should not set it up.
    if (wrappedDefaultComm_.is_null ()) {
      Teuchos::RCP<const Teuchos::Comm<int> > comm;
#ifdef HAVE_TPETRACORE_MPI
      // Teuchos::MpiComm's constructor used to invoke MPI collectives.
      // It still reserves the right to do so.  This means MPI must be
      // initialized and not finalized.
      const bool mpiReady = mpiIsInitializedAndNotFinalized ();
      if (mpiReady) {
        comm = Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
      }
      else {
        comm = Teuchos::rcp (new Teuchos::SerialComm<int> ());
      }
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
#if defined(HAVE_TPETRACORE_MPI)
      initMpiIfNeeded (argc, argv);
      // It's technically legal to initialize Tpetra after
      // MPI_Finalize has been called.  This means that we can't call
      // MPI_Comm_rank without first checking MPI_Finalized.
      const int myRank = getRankHarmlessly (MPI_COMM_WORLD);
#else
      const int myRank = 0;
#endif // defined(HAVE_TPETRACORE_MPI)
      initKokkosIfNeeded (argc, argv, myRank);
    }
    tpetraIsInitialized_ = true;
  }

#ifdef HAVE_TPETRACORE_MPI
  void initialize (int* argc, char*** argv, MPI_Comm comm)
  {
    if (! tpetraIsInitialized_) {
#if defined(HAVE_TPETRACORE_MPI)
      initMpiIfNeeded (argc, argv);
      // It's technically legal to initialize Tpetra after
      // MPI_Finalize has been called.  This means that we can't call
      // MPI_Comm_rank without first checking MPI_Finalized.
      const int myRank = getRankHarmlessly (comm);
#else
      const int myRank = 0;
#endif // defined(HAVE_TPETRACORE_MPI)
      initKokkosIfNeeded (argc, argv, myRank);
    }
    tpetraIsInitialized_ = true;

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
    if (! tpetraIsInitialized_) {
#if defined(HAVE_TPETRACORE_MPI)
      initMpiIfNeeded (argc, argv);
#endif // defined(HAVE_TPETRACORE_MPI)
      // It's technically legal to initialize Tpetra after
      // MPI_Finalize has been called.  This means that we can't call
      // MPI_Comm_rank without first checking MPI_Finalized.
      const int myRank = comm->getRank ();
      initKokkosIfNeeded (argc, argv, myRank);
    }
    tpetraIsInitialized_ = true;
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

  ScopeGuard::ScopeGuard (int* argc, char*** argv)
  {
    ::Tpetra::initialize (argc, argv);
  }

#ifdef HAVE_TPETRA_MPI
  ScopeGuard::ScopeGuard (int* argc, char*** argv, MPI_Comm comm)
  {
    ::Tpetra::initialize (argc, argv, comm);
  }
#endif // HAVE_TPETRA_MPI

  ScopeGuard::~ScopeGuard ()
  {
    ::Tpetra::finalize ();
  }

} // namespace Tpetra
