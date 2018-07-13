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

#ifndef TPETRA_DETAILS_ISINTERCOMM_HPP
#define TPETRA_DETAILS_ISINTERCOMM_HPP

/// \file Tpetra_Details_isInterComm.hpp
/// \brief Declaration of Tpetra::Details::isInterComm
///
/// \warning This file and its contents are implementation details of
///   Tpetra.  Users must not rely on them.
///
/// Tpetra::Details::isInterComm wraps MPI_Comm_test_inter.  The
/// Tpetra::Details::isInterComm function is the only thing in this
/// file upon which Tpetra developers should rely.  Tpetra developers
/// should not rely on anything else in this file.  <i>Users</i> may
/// not rely on <i>anything</i> in this file!

#include "TpetraCore_config.h"
#include "Teuchos_Comm.hpp"

namespace Tpetra {
namespace Details {

/// \brief Return true if and only if the input communicator wraps an
///   MPI intercommunicator.
///
/// The most common MPI communicators are intracommunicators
/// ("in<i>tra</i>," not "in<i>ter</i>").  This includes
/// MPI_COMM_WORLD, MPI_COMM_SELF, and the results of MPI_Comm_dup and
/// MPI_Comm_split.  Intercommunicators come from
/// MPI_Intercomm_create.
///
/// This distinction matters because only collectives over
/// intracommunicators may use MPI_IN_PLACE, to let the send and
/// receive buffers alias each other.  Collectives over
/// intercommunicators may <i>not</i> use MPI_IN_PLACE.
///
/// \param comm [in] The input communicator.
///
/// \return Whether the input communicator wraps an MPI
///   intercommunicator.
bool
isInterComm (const Teuchos::Comm<int>& comm);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ISINTERCOMM_HPP
