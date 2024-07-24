// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_HPP
#define TPETRA_CORE_HPP

/// \file Tpetra_Core.hpp
/// \brief Functions for initializing and finalizing Tpetra.
///
/// This file declares functions for initializing (setting up) and
/// finalizing (tearing down) Tpetra.  All overloads of initialize()
/// automatically initialize MPI (if Trilinos was built with MPI
/// enabled) and Kokkos, if necessary.  The finalize() function
/// finalizes MPI (by calling MPI_Finalize) and Kokkos, if necessary.

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_Comm.hpp>
#ifdef HAVE_TPETRACORE_MPI
#  include "mpi.h"
#endif // HAVE_TPETRACORE_MPI

namespace Tpetra {

  /// \brief Get Tpetra's default communicator.
  ///
  /// \pre One of the Tpetra::initialize() functions has been called.
  ///
  /// \return If one of the versions of initialize() was called that
  ///   takes a default communicator, this function returns that
  ///   communicator.  Otherwise, this function returns MPI_COMM_WORLD
  ///   (wrapped in a Teuchos wrapper) if Trilinos was built with MPI
  ///   enabled, or a Teuchos::SerialComm instance otherwise.
  Teuchos::RCP<const Teuchos::Comm<int> > getDefaultComm ();

  /// \brief Whether Tpetra is in an initialized state.
  ///
  /// Initialize Tpetra by calling one of the versions of
  /// initialize().  After initialize() returns, Tpetra is
  /// initialized.  Once finalize() returns, Tpetra is no longer
  /// initialized.
  bool isInitialized ();

  /// \brief Initialize Tpetra.
  ///
  /// This initializes the following if they have not already been
  /// initialized:
  /// <ul>
  ///   <li> MPI (if Trilinos was built with MPI enabled) </li>
  ///   <li> Kokkos </li>
  /// </ul>
  ///
  /// If Trilinos was built with MPI enabled, this function sets the
  /// default communicator to MPI_COMM_WORLD (wrapped in a Teuchos
  /// wrapper).  Otherwise, it sets the default communicator to a
  /// Teuchos::SerialComm instance.
  ///
  /// \param argc [in/out] Same as first argument of MPI_Init()
  /// \param argv [in/out] Same as second argument of MPI_Init()
  ///
  /// The \c argc and \c argv arguments are both passed by pointer, in
  /// order to match MPI_Init's interface.  MPI_Init() reserves the
  /// right to modify command-line arguments, e.g., by reading and
  /// removing those that pertain to MPI.  Thus, in main(), one would
  /// write
  /// \code
  /// Tpetra::initialize (&argc, &argc);
  /// \endcode
  void initialize (int* argc, char*** argv);

#ifdef HAVE_TPETRA_MPI
  /// \brief Initialize Tpetra.
  ///
  /// This version of initialize() only exists if Trilinos was built
  /// with MPI enabled.
  ///
  /// This function initializes the following if they have not already
  /// been initialized:
  /// <ul>
  ///   <li> MPI (if Trilinos was built with MPI enabled) </li>
  ///   <li> Kokkos </li>
  /// </ul>
  ///
  /// \param argc [in/out] Same as first argument of MPI_Init()
  /// \param argv [in/out] Same as second argument of MPI_Init()
  /// \param comm [in] Tpetra's default communicator
  ///
  /// The \c argc and \c argv arguments are both passed by pointer, in
  /// order to match MPI_Init's interface.  MPI_Init() reserves the
  /// right to modify command-line arguments, e.g., by reading and
  /// removing those that pertain to MPI.  Thus, in main(), one would
  /// write
  /// \code
  /// MPI_Comm comm = ...; // whatever you want
  /// Tpetra::initialize (&argc, &argc, comm);
  /// \endcode
  void initialize (int* argc, char*** argv, MPI_Comm comm);
#endif // HAVE_TPETRA_MPI

  /// \brief Initialize Tpetra.
  ///
  /// \warning It is NOT legal to create a Teuchos::MpiComm until
  ///   after MPI has been initialized.  Thus, you may <i>not</i> call
  ///   this function with a Teuchos::MpiComm instance unless MPI has
  ///   been initialized.
  ///
  /// This initializes the following if they have not already been
  /// initialized:
  /// <ul>
  ///   <li> MPI (if Trilinos was built with MPI enabled) </li>
  ///   <li> Kokkos </li>
  /// </ul>
  ///
  /// \param argc [in/out] Same as first argument of MPI_Init()
  /// \param argv [in/out] Same as second argument of MPI_Init()
  /// \param comm [in] Tpetra's default communicator, wrapped in a
  ///   Teuchos wrapper.  This may be either a Teuchos::MpiComm or a
  ///   Teuchos::SerialComm instance.
  ///
  /// The \c argc and \c argv arguments are both passed by pointer, in
  /// order to match MPI_Init's interface.  MPI_Init() reserves the
  /// right to modify command-line arguments, e.g., by reading and
  /// removing those that pertain to MPI.  Thus, in main(), one would
  /// write
  /// \code
  /// Teuchos::RCP<const Teuchos::Comm<int> > comm = ...; // whatever you want
  /// Tpetra::initialize (&argc, &argc, comm);
  /// \endcode
  void
  initialize (int* argc, char*** argv,
              const Teuchos::RCP<const Teuchos::Comm<int> >& comm);

  /// \brief Finalize Tpetra.
  ///
  /// If Tpetra::initialize initialized Kokkos, finalize Kokkos.  If
  /// Tpetra::initialize initialized MPI, finalize MPI.  Don't call
  /// this unless you first call Tpetra::initialize.
  ///
  /// If you (the user) initialized Kokkos resp. MPI before
  /// Tpetra::initialize was called, then this function does NOT
  /// finalize Kokkos resp. MPI.  In that case, you (the user) are
  /// responsible for finalizing Kokkos resp. MPI.
  void finalize ();

  /// \brief Scope guard whose destructor automatically calls
  ///   Tpetra::finalize for you.
  ///
  /// This class' constructor does the same thing as
  /// Tpetra::initialize (see above).  Its destructor automatically
  /// calls Tpetra::finalize.  This ensures correct Tpetra
  /// finalization even if intervening code throws an exception.
  ///
  /// Compare to Kokkos::ScopeGuard and Teuchos::GlobalMPISession.
  ///
  /// Always give the ScopeGuard instance a name.  Otherwise, you'll
  /// create a temporary object whose destructor will be called right
  /// away.  That's not what you want.
  ///
  /// Here is an example of how to use this class:
  /// \code
  /// #include "Tpetra_Core.hpp"
  /// #include "Tpetra_Map.hpp"
  ///
  /// int main (int argc, char* argv[]) {
  ///   Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  ///
  ///   // Never create Tpetra or Kokkos objects (other than
  ///   // the ScopeGuard object itself) at main scope.
  ///   // Otherwise, their destructors will be called after
  ///   // MPI_Finalize and Kokkos::finalize are called, which
  ///   // is forbidden.
  ///   {
  ///     auto comm = Tpetra::getDefaultComm ();
  ///     using GO = Tpetra::Map<>::global_ordinal_type;
  ///     const GO gblNumInds = 1000;
  ///     const GO indexBase = 0;
  ///     Tpetra::Map<> map (gblNumInds, indexBase, comm,
  ///                        Tpetra::GloballyDistributed);
  ///     // ... code that uses map ...
  ///   }
  ///   return EXIT_SUCCESS;
  /// }
  /// \endcode
  class ScopeGuard {
  public:
    /// \brief Initialize Tpetra.
    ///
    /// If MPI_Init has not yet been called, then call it.  If
    /// Kokkos::initialize has not yet been called, then call it.
    ///
    /// \param argc [in/out] Address of the first argument to main().
    /// \param argv [in/out] Address of the second argument to main().
    ScopeGuard (int* argc, char*** argv);

#ifdef HAVE_TPETRA_MPI
    /// \brief Initialize Tpetra, and set Tpetra's default MPI communicator.
    ///
    /// Assume that MPI_Init has been called.  If Kokkos::initialize
    /// has not yet been called, then call it.  Make the input MPI
    /// communicator Tpetra's default MPI communicator.
    ///
    /// \param argc [in/out] Address of the first argument to main().
    /// \param argv [in/out] Address of the second argument to main().
    /// \param comm [in] Default MPI communicator for Tpetra.  The
    ///   caller is responsible for calling MPI_Comm_free on this
    ///   after use, if necessary.
    ScopeGuard (int* argc, char*** argv, MPI_Comm comm);
#endif // HAVE_TPETRA_MPI

    /// \brief Default constructor (FORBIDDEN)
    ///
    /// You must give ScopeGuard's constructor argc and argv.
    /// Use the constructor above this one.
    ScopeGuard () = delete;

    /// \brief Finalize Tpetra.
    ///
    /// If the constructor called Kokkos::initialize, then call
    /// Kokkos::finalize.  If the constructor called MPI_Init, then
    /// call MPI_Finalize.
    ~ScopeGuard ();
  };

} // namespace Tpetra

#endif // TPETRA_CORE_HPP
