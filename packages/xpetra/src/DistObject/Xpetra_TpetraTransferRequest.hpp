// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_TPETRATRANSFERREQUEST_DECL_HPP
#define XPETRA_TPETRATRANSFERREQUEST_DECL_HPP

#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_DistObject.hpp"
#include "Tpetra_TransferRequest.hpp"


namespace Xpetra {


  template <class Packet,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class TpetraTransferRequest
    : public virtual TransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node >
  {

    // The following typedef are used by the XPETRA_DYNAMIC_CAST() macro.
    typedef TpetraTransferRequest<Packet,LocalOrdinal,GlobalOrdinal,Node> TpetraTransferRequestClass;

  public:

    //! @name Constructors and destructor
    //@{

    //! Basic constuctor.
    TpetraTransferRequest(Teuchos::RCP<Tpetra::TransferRequest<Packet,LocalOrdinal, GlobalOrdinal, Node > > &tr)
      : tr_ (tr)
    {}

    virtual ~TpetraTransferRequest() = default;

    void process() {tr_->process(); };


  private:
    //! The Tpetra::TransferRequest which this class wraps.
    RCP< Tpetra::TransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node> > tr_;

  }; // TpetraTransferRequest class


  // Things we actually need

  // Things we actually need
  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<TransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node > > toXpetra(RCP<Tpetra::TransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node > > tr) {
    if (!tr.is_null())
      return rcp(new TpetraTransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node >(tr));

    return Teuchos::null;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const TransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node > > toXpetra(RCP<const Tpetra::TransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node > > tr) {
    if (!tr.is_null())
      return rcp(new TpetraTransferRequest<Packet, LocalOrdinal, GlobalOrdinal, Node >(tr));

    return Teuchos::null;
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  const Tpetra::TransferRequest< Packet,LocalOrdinal, GlobalOrdinal, Node> & toTpetra(const TransferRequest< Packet,LocalOrdinal, GlobalOrdinal, Node> &xTr) {
    typedef TpetraTransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node > TpetraTransferRequestClass;
    XPETRA_DYNAMIC_CAST(const TpetraTransferRequestClass, xTr, tTr, "toTpetra");
    return *tTr.getTpetra_TransferRequest();
  }

  template <class Packet, class LocalOrdinal, class GlobalOrdinal, class Node>
  Tpetra::TransferRequest< Packet,LocalOrdinal, GlobalOrdinal, Node> & toTpetra(TransferRequest< Packet,LocalOrdinal, GlobalOrdinal, Node> &xTr) {
    typedef TpetraTransferRequest< Packet, LocalOrdinal, GlobalOrdinal, Node > TpetraTransferRequestClass;
    XPETRA_DYNAMIC_CAST(      TpetraTransferRequestClass, xTr, tTr, "toTpetra");
    return *tTr.getTpetra_TransferRequest();
  }

} // Xpetra namespace

#define XPETRA_TPETRATRANSFERREQUEST_SHORT
#endif // XPETRA_TPETRATRANSFERREQUEST_HPP
