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

#ifndef TPETRA_DETAILS_GETENTRYONHOST_HPP
#define TPETRA_DETAILS_GETENTRYONHOST_HPP

/// \file Tpetra_Details_getEntryOnHost.hpp
/// \brief Declaration and definition of
///   Tpetra::Details::getEntryOnHost.
/// \warning The contents of this file are implementation details of
///   Tpetra.  We make no promises of backwards compatibility.

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

template<class ViewType,
         class IndexType>
typename ViewType::non_const_value_type
getEntryOnHost (const ViewType& x,
                const IndexType ind)
{
  static_assert (ViewType::Rank == 1, "x must be a rank-1 Kokkos::View.");
  // Get a 0-D subview of the entry of the array, and copy to host scalar.
  typename ViewType::non_const_value_type val;
  Kokkos::deep_copy(val, Kokkos::subview(x, ind));
  return val;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_GETENTRYONHOST_HPP
