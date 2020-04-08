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

#ifndef TPETRA_TRANSFERREQUEST_DECL_HPP
#define TPETRA_TRANSFERREQUEST_DECL_HPP

/// \file Tpetra_TransferRequest_decl.hpp
/// \brief Declaration of the Tpetra::TransferRequest class
///
/// If you want to use Tpetra::TransferRequest, include
/// "Tpetra_TransferRequest.hpp" (a file which CMake generates and installs
/// for you).  If you only want the declaration of Tpetra::TransferRequest,
/// include this file (Tpetra_TransferRequest_decl.hpp).

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_DistObject_fwd.hpp"
#include "Kokkos_ArithTraits.hpp"


namespace Tpetra {

  template <class Packet,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class TransferRequest {
  public:
    //! @name Typedefs
    //@{

    /// \brief The type of each datum being sent or received in an Import or Export.
    ///
    /// Note that this type does not always correspond to the
    /// <tt>Scalar</tt> template parameter of subclasses.
    using packet_type = typename ::Kokkos::Details::ArithTraits<Packet>::val_type;
    //! The type of local indices.
    using local_ordinal_type = LocalOrdinal;
    //! The type of global indices.
    using global_ordinal_type = GlobalOrdinal;
    //! The Node type.  If you don't know what this is, don't use it.
    using node_type = Node;

    using dist_object_type = DistObject<Packet, LocalOrdinal, GlobalOrdinal, Node>;

    using transfer_type = ::Tpetra::Details::Transfer<local_ordinal_type, global_ordinal_type, node_type>;

    using rev_op_type = typename dist_object_type::ReverseOption;

    TransferRequest(Teuchos::RCP<dist_object_type>& dest,
                    Teuchos::RCP<const transfer_type>& transfer,
                    rev_op_type revOp,
                    CombineMode CM,
                    bool needCommunication);

    void process();

    virtual ~TransferRequest () = default;

  private:


    Teuchos::RCP<dist_object_type> dest_;
    Teuchos::RCP<const transfer_type> transfer_;
    rev_op_type revOp_;
    CombineMode CM_;
    bool needCommunication_;
    bool hasBeenProcessed_;


  }; // class TransferRequest
} // namespace Tpetra

#endif // TPETRA_TRANSFERREQUEST_DECL_HPP
