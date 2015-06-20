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
  ///
  ///   - MPI (if Trilinos was built with MPI enabled)
  ///   - The default communicator (returned by getDefaultComm())
  ///   - Kokkos
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
  ///
  ///   - MPI (if Trilinos was built with MPI enabled)
  ///   - The default communicator (returned by getDefaultComm())
  ///   - Kokkos
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
  /// This initializes the following if they have not already been
  /// initialized:
  ///
  ///   - MPI (if Trilinos was built with MPI enabled)
  ///   - The default communicator (returned by getDefaultComm())
  ///   - Kokkos
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

} // namespace Tpetra
