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

#ifndef TPETRA_DETAILS_ASSUMEMPIISCUDAAWARE_HPP
#define TPETRA_DETAILS_ASSUMEMPIISCUDAAWARE_HPP

// Forward declarations cause build errors, depending on what other
// header files get included first.  Thus, we must include the header
// file here.
#include <ostream>

namespace Tpetra {
namespace Details {

/// \brief Whether to assume that MPI is CUDA aware.
///
/// See Trilinos GitHub issues #1571 and #1088 to learn what it means
/// for an MPI implementation to be "CUDA aware," and why this matters
/// for performance.
///
/// \param out [out] If not NULL, print human-readable output to
///   <tt>*out</tt>, describing what this function is doing.  Users
///   are responsible for managing output with multiple MPI processes
///   (e.g., only make this non-NULL on Process 0).
///
/// \return Whether to assume that the MPI implementation that
///   Trilinos uses is CUDA aware.
bool
assumeMpiIsCudaAware (std::ostream* out);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ASSUMEMPIISCUDAAWARE_HPP
