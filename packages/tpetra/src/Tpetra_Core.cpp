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

#include <Tpetra_Core.hpp>
#ifdef HAVE_TPETRA_MPI
#  include <Teuchos_DefaultMpiComm.hpp>// this includes mpi.h too
#else
#  include <Teuchos_DefaultSerialComm.hpp>
#endif // HAVE_TPETRA_MPI
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
#  include <Kokkos_Core.hpp>
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

namespace Tpetra {
  namespace { // (anonymous)
    // Whether one of the Tpetra::initialize() functions has been called before.
    bool tpetraIsInitialized_ = false;

    // Tpetra's default communicator, wrapped in a Teuchos wrapper.
    // After Tpetra::finalize() is called, this GOES AWAY (is set to null).
    Teuchos::RCP<const Teuchos::Comm<int> > wrappedDefaultComm_;

#ifdef HAVE_TPETRA_MPI
    // Initialize MPI, if needed, and check for errors.
    void initMpi (int* argc, char*** argv)
    {
      int isInitialized = 0;
      int err = MPI_Initialized (&isInitialized);

      TEUCHOS_TEST_FOR_EXCEPTION(
        err != MPI_SUCCESS, std::runtime_error, "MPI_Initialized failed with "
        "error code " << err << " != MPI_SUCCESS.  This probably indicates "
        "that your MPI library is corrupted or that it is incorrectly linked "
        "to your program, since this function should otherwise always "
        "succeed.");

      if (isInitialized == 0) { // MPI not yet initialized
        // Tpetra doesn't currently need to call MPI_Init_thread, since
        // with Tpetra, only one thread ever calls MPI functions.  If we
        // ever want to explore MPI_THREAD_MULTIPLE, here would be the
        // place to call MPI_Init_thread.
        err = MPI_Init (argc, argv);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
        err != MPI_SUCCESS, std::runtime_error, "MPI_Init failed with error "
        "code " << err << " != MPI_SUCCESS.  If MPI was set up correctly, this "
        "should not happen, since we have already checked that MPI_Init (or "
        "MPI_Init_thread) has not yet been called.  This may indicate that "
        "your MPI library is corrupted or that it is incorrectly linked to "
        "your program.");
    }
#endif // HAVE_TPETRA_MPI

  } // namespace (anonymous)

  bool isInitialized () {
    return tpetraIsInitialized_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      tpetraIsInitialized_, std::runtime_error, "Tpetra::getDefaultComm: "
      "You must call Tpetra::initialize before you may get a default "
      "communicator.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      wrappedDefaultComm_.is_null (), std::logic_error,
      "Tpetra::getDefaultComm: wrappedDefaultComm_ is null!  "
      "This should never happen, because Tpetra claims that it has been "
      "initialized.  Please report this bug to the Tpetra developers.");
    return wrappedDefaultComm_;
  }


  void initialize (int* argc, char*** argv)
  {
#ifdef HAVE_TPETRA_MPI
    initMpi (argc, argv);
#endif // HAVE_TPETRA_MPI

    if (! tpetraIsInitialized_) {
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
      // Unlike MPI_Init, Kokkos promises not to modify argc and argv.
      Kokkos::initialize (*argc, *argv);
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

#ifdef HAVE_TPETRA_MPI
      wrappedDefaultComm_ =
        Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
      wrappedDefaultComm_ = Teuchos::rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_TPETRA_MPI
      tpetraIsInitialized_ = true;
    }
  }


#ifdef HAVE_TPETRA_MPI
  void initialize (int* argc, char*** argv, MPI_Comm comm)
  {
#ifdef HAVE_TPETRA_MPI
    initMpi (argc, argv);
#endif // HAVE_TPETRA_MPI

    // Set the default communicator.  What if users have already
    // called initialize() before, but with a different default
    // communicator?  There are two possible things we could do here:
    //
    //   1. Test via MPI_Comm_compare whether comm differs from the
    //      raw MPI communicator in wrappedDefaultComm_ (if indeed it
    //      is an MpiComm).
    //   2. Accept that the user might want to change the default
    //      communicator, and let them do it.
    //
    // I prefer #2.  Perhaps it would be sensible to print a warning
    // here, but on which process?  Would we use the old or the new
    // communicator to find that process' rank?  We don't want to
    // use MPI_COMM_WORLD's Process 0, since neither communicator
    // might include that process.  Furthermore, in some
    // environments, only Process 0 in MPI_COMM_WORLD is allowed to
    // do I/O.  Thus, we just let the change go without a warning.
    wrappedDefaultComm_ = Teuchos::rcp (new Teuchos::MpiComm<int> (comm));

    if (! tpetraIsInitialized_) {
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
      // Unlike MPI_Init, Kokkos promises not to modify argc and argv.
      Kokkos::initialize (*argc, *argv);
#endif // TPETRA_HAVE_KOKKOS_REFACTOR
      tpetraIsInitialized_ = true;
    }
  }
#endif // HAVE_TPETRA_MPI


  void
  initialize (int* argc, char*** argv,
              const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
  {
#ifdef HAVE_TPETRA_MPI
    initMpi (argc, argv);
#endif // HAVE_TPETRA_MPI

    // Set the default communicator.  What if users have already
    // called initialize() before, but with a different default
    // communicator?  There are two possible things we could do here:
    //
    //   1. Test via MPI_Comm_compare whether comm differs from the
    //      raw MPI communicator in wrappedDefaultComm_ (if indeed it
    //      is an MpiComm).
    //   2. Accept that the user might want to change the default
    //      communicator, and let them do it.
    //
    // I prefer #2.  Perhaps it would be sensible to print a warning
    // here, but on which process?  Would we use the old or the new
    // communicator to find that process' rank?  We don't want to
    // use MPI_COMM_WORLD's Process 0, since neither communicator
    // might include that process.  Furthermore, in some
    // environments, only Process 0 in MPI_COMM_WORLD is allowed to
    // do I/O.  Thus, we just let the change go without a warning.
    wrappedDefaultComm_ = comm;

    if (! tpetraIsInitialized_) {
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
      // Unlike MPI_Init, Kokkos promises not to modify argc and argv.
      Kokkos::initialize (*argc, *argv);
#endif // TPETRA_HAVE_KOKKOS_REFACTOR
      tpetraIsInitialized_ = true;
    }
  }


  void finalize ()
  {
    using std::cerr;
    using std::endl;

    if (! tpetraIsInitialized_) {
      // If the user didn't call initialize(), then we shouldn't call
      // MPI_Finalize() or otherwise shut anything else off.  For
      // example, the user might want to use Kokkos without Tpetra, so
      // they might have called Kokkos::initialize() without calling
      // Tpetra::initialize().  In that case, we shouldn't call
      // Kokkos::initialize() here.
      return;
    }

#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
    Kokkos::finalize ();
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

    // Make sure that no outstanding references to the communicator
    // remain.  If users gave initialize() an MPI_Comm, _they_ are
    // responsible for freeing it before calling finalize().
    wrappedDefaultComm_ = Teuchos::null;

#ifdef HAVE_TPETRA_MPI
    int isInitialized = 0;
    int err = MPI_Initialized (&isInitialized);

    // finalize() is a kind of destructor, so it's a bad idea to throw
    // an exception.  Instead, just print out an error message and
    // hope for the best.
    if (err != MPI_SUCCESS) {
      cerr << "MPI_Initialized failed with error code " << err << " != "
        "MPI_SUCCESS.  This probably indicates that your MPI library is "
        "corrupted or that it is incorrectly linked to your program, since "
        "this function should otherwise always succeed." << endl;
    }
    else if (isInitialized != 0) {
      // This must be called by the same thread that called MPI_Init
      // (possibly, but not necessarily, in Tpetra::initialize()).
      err = MPI_Finalize ();

      // finalize() is a kind of destructor, so it's a bad idea to
      // throw an exception.  Instead, just print out an error message
      // and hope for the best.
      if (err != MPI_SUCCESS) {
        cerr << "MPI_Finalize failed with error code " << err << " != "
          "MPI_SUCCESS.  I'm not sure how to recover from this safely, "
          "so I will keep going and hope for the best." << endl;
      }
    }
#endif // HAVE_TPETRA_MPI

    tpetraIsInitialized_ = false; // it's not anymore.
  }
} // namespace Tpetra
