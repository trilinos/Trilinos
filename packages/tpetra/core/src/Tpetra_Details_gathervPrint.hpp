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
// ************************************************************************
// @HEADER
*/

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
