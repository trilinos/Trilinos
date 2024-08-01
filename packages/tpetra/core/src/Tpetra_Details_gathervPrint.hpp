// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tpetra_Details_gathervPrint.hpp
/// \brief Declaration of a function that prints strings from each process.
///
/// \warning This is an implementation detail of Tpetra.  Users may
///   not rely on this header file or any declarations or definitions
///   in it.  They may disappear or change at any time.

#ifndef TPETRA_DETAILS_GATHERPRINT_HPP
#define TPETRA_DETAILS_GATHERPRINT_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_Comm.hpp>
#include <string>
#include <ostream>

namespace Tpetra {
namespace Details {

/// \brief On Process 0 in the given communicator, print strings from
///   each process in that communicator, in rank order.
///
/// For each process in the given communicator \c comm, send its
/// string \c s to Process 0 in that communicator.  Process 0 prints
/// the strings in rank order.
///
/// This is a collective over the given communicator \c comm.  Process
/// 0 promises not to store all the strings in its memory.  This
/// function's total memory usage on any process is proportional to
/// the calling process' string length, plus the max string length
/// over any process.  This does NOT depend on the number of processes
/// in the communicator.  Thus, we call this a "memory-scalable"
/// operation.  While the function's name suggests MPI_Gatherv, the
/// implementation may NOT use MPI_Gather or MPI_Gatherv, because
/// neither of those are not memory scalable.
///
/// Process 0 prints nothing other than what is in the string.  It
/// does not add an endline after each string, nor does it identify
/// each string with its owning process' rank.  If you want either of
/// those in the string, you have to put it there yourself.
///
/// \param out [out] The output stream to which to write.  ONLY
///   Process 0 in the given communicator will write to this.  Thus,
///   this stream need only be valid on Process 0.
///
/// \param s [in] The string to write.  Each process in the given
///   communicator has its own string.  Strings may be different on
///   different processes.  Zero-length strings are OK.
///
/// \param comm [in] The communicator over which this operation is a
///   collective.
void
gathervPrint (std::ostream& out,
              const std::string& s,
              const Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GATHERPRINT_HPP
