// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_MPI_COMM_HPP
#define TEUCHOS_MPI_COMM_HPP

/// \file Teuchos_DefaultMpiComm.hpp
/// \brief Implementation of Teuchos wrappers for MPI.
///
/// \warning It only makes sense to include this file if MPI is enabled.

#include <Teuchos_ConfigDefs.hpp>

// If MPI is not enabled, disable the contents of this file.
#ifdef HAVE_TEUCHOS_MPI

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommUtilities.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_MpiReductionOpSetter.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Assert.hpp"
#include <mpi.h>
#include <iterator>
#include "Teuchos_DefaultMpiComm_decl.hpp"

// This must be defined globally for the whole program!
//#define TEUCHOS_MPI_COMM_DUMP

#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif

namespace Teuchos {

//! Human-readable string version of the given MPI error code.
TEUCHOSCOMM_LIB_DLL_EXPORT std::string
mpiErrorCodeToString (const int err);

namespace details {
  /// \brief Give \c comm to MPI_Comm_free, if MPI is not yet finalized.
  ///
  /// This function "frees" the given communicator by giving it to
  /// MPI_Comm_free.  It only does so if MPI_Finalize has not yet been
  /// called.  If MPI_Finalize has been called, this function does
  /// nothing, since it is not legal to call most MPI functions after
  /// MPI_Finalize has been called.  This function also ignores any
  /// errors returned by MPI_Finalize, which makes it suitable for use
  /// in a destructor.
  ///
  /// \note This function may allow a memory leak in your program, if
  ///   you have allowed the MPI_Comm to persist after MPI_Finalize
  ///   has been called.
  TEUCHOSCOMM_LIB_DLL_EXPORT void safeCommFree (MPI_Comm* comm);

  /// Set the given communicator's error handler to \c handler.
  ///
  /// If the MPI version is >= 2, this calls MPI_Comm_set_handler().
  /// If the MPI version is 1, this calls MPI_Errhandler_set().
  TEUCHOSCOMM_LIB_DLL_EXPORT int setCommErrhandler (MPI_Comm comm, MPI_Errhandler handler);

} // namespace details

#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal, typename T>
void dumpBuffer(
  const std::string &funcName, const std::string &buffName
  ,const Ordinal bytes, const T buff[]
  )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out
    << "\n" << funcName << "::" << buffName << ":\n";
  tab.incrTab();
  for( Ordinal i = 0; i < bytes; ++i ) {
    *out << buffName << "[" << i << "] = '" << buff[i] << "'\n";
  }
  *out << "\n";
}
#endif // TEUCHOS_MPI_COMM_DUMP

}

#endif // HAVE_TEUCHOS_MPI
#endif // TEUCHOS_MPI_COMM_HPP

