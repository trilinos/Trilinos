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

#ifndef TPETRA_TRANSFERREQUEST_DEF_HPP
#define TPETRA_TRANSFERREQUEST_DEF_HPP

/// \file Tpetra_TransferRequest_def.hpp
/// \brief Definition of the Tpetra::TransferRequest class
///
/// If you want to use Tpetra::TransferRequest, include
/// "Tpetra_TransferRequest.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::TransferRequest,
/// include "Tpetra_TransferRequest_decl.hpp".


#include "Tpetra_DistObject.hpp"

namespace Tpetra {

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  TransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  TransferRequest (Teuchos::RCP<dist_object_type>& dest,
                   Teuchos::RCP<const transfer_type>& transfer,
                   rev_op_type revOp,
                   CombineMode CM,
                   bool needCommunication) :
    dest_ (dest),
    transfer_ (transfer),
    revOp_ (revOp),
    CM_ (CM),
    needCommunication_ (needCommunication),
    hasBeenProcessed_ (false)
  {
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  TransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node>::
  process()
  {
    if (!hasBeenProcessed_) {
      dest_->finishTransfer(*transfer_, revOp_, CM_, needCommunication_);
      dest_ = Teuchos::null;
      transfer_ = Teuchos::null;
      hasBeenProcessed_ = true;
    }
  }

// Explicit instantiation macro for general TransferRequest.
#define TPETRA_TRANSFERREQUEST_INSTANT(SCALAR, LO, GO, NODE) \
  template class TransferRequest< SCALAR , LO , GO , NODE >;

// Explicit instantiation macro for TransferRequest<char, ...>.
// The "SLGN" stuff above doesn't work for Packet=char.
#define TPETRA_TRANSFERREQUEST_INSTANT_CHAR(LO, GO, NODE) \
  template class TransferRequest< char , LO , GO , NODE >;


} // namespace Tpetra

#endif // TPETRA_DISTOBJECT_DEF_HPP
