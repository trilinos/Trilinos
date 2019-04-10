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

#ifndef TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DECL_HPP
#define TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DECL_HPP

#include "Tpetra_BlockMultiVector_decl.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_packBlockMultiVector_decl.hpp
/// \brief Declaration of function that packs the entries of a
///   Tpetra::BlockMultiVector for communication.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \brief Implementation of BlockMultiVector::copyAndPermuteNew.
///
/// \tparam LO The type of local indices.
/// \tparam GO The type of global indices.
/// \tparam NT The Node type.
///
/// \param tgt [in] "Target" object of an Export or Import operation
///   on a BlockMultiVector.
/// \param src [in] "Source" object of an Export or Import operation
///   on a BlockMultiVector.
/// \param numSameIDs [in] See DistObject::copyAndPermuteNew.
/// \param permuteToLIDs [in] See DistObject::copyAndPermuteNew.
/// \param permuteFromLIDs [in] See DistObject::copyAndPermuteNew.
///
/// \return (Error code, pointer to error message).  If no error, then
///   error code is zero and error message pointer is null.
template<class SC, class LO, class GO, class NT>
std::pair<int, std::unique_ptr<std::string>>
copyAndPermuteBlockMultiVector
  (BlockMultiVector<SC, LO, GO, NT>& tgt,
   const BlockMultiVector<SC, LO, GO, NT>& src,
   const size_t numSameIDs,
   const Kokkos::DualView<
     const LO*,
     typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
   >& permuteToLIDs,
   const Kokkos::DualView<
     const LO*,
     typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
   >& permuteFromLIDs);

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DECL_HPP
