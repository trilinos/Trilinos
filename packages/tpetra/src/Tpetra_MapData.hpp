// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_MAPDATA_HPP
#define TPETRA_MAPDATA_HPP

#include "Tpetra_MapDataDecl.hpp"

namespace Tpetra {

  template<class LocalOrdinal,class GlobalOrdinal>
  MapData<LocalOrdinal,GlobalOrdinal>::MapData(
            Teuchos_Ordinal indexBase, 
            GlobalOrdinal numGlobalEntries,
            LocalOrdinal numMyEntries,
            GlobalOrdinal minAllGID,
            GlobalOrdinal maxAllGID,
            GlobalOrdinal minMyGID,
            GlobalOrdinal maxMyGID,
            const Teuchos::ArrayRCP<GlobalOrdinal> &lgMap,
            const std::map<GlobalOrdinal,LocalOrdinal>  &glMap,
            bool contiguous,
            Teuchos::RCP<const Teuchos::Comm<int> > comm) 
      : comm_(comm)
      , numGlobalEntries_(numGlobalEntries)
      , indexBase_(indexBase)
      , numMyEntries_(numMyEntries)
      , minMyGID_(minMyGID)
      , maxMyGID_(maxMyGID)
      , minAllGID_(minAllGID)
      , maxAllGID_(maxAllGID)
      , contiguous_(contiguous)
      , distributed_(checkIsDist())
      , lgMap_(lgMap)
      , glMap_(glMap)
  {}

  template<class LocalOrdinal,class GlobalOrdinal>
  MapData<LocalOrdinal,GlobalOrdinal>::MapData(
            Teuchos_Ordinal indexBase, 
            GlobalOrdinal numGlobalEntries,
            LocalOrdinal numMyEntries,
            GlobalOrdinal minAllGID,
            GlobalOrdinal maxAllGID,
            GlobalOrdinal minMyGID,
            GlobalOrdinal maxMyGID,
            const Teuchos::ArrayRCP<GlobalOrdinal> &lgMap,
            const std::map<GlobalOrdinal,LocalOrdinal>  &glMap,
            bool contiguous,
            Teuchos::RCP<const Teuchos::Comm<int> > comm,
            bool isLocal)
      : comm_(comm)
      , numGlobalEntries_(numGlobalEntries)
      , indexBase_(indexBase)
      , numMyEntries_(numMyEntries)
      , minMyGID_(minMyGID)
      , maxMyGID_(maxMyGID)
      , minAllGID_(minAllGID)
      , maxAllGID_(maxAllGID)
      , contiguous_(contiguous)
      , distributed_(!isLocal)
      , lgMap_(lgMap)
      , glMap_(glMap)
    {}

  template<class LocalOrdinal,class GlobalOrdinal>
  MapData<LocalOrdinal,GlobalOrdinal>::~MapData() {}

  template<class LocalOrdinal,class GlobalOrdinal>
  bool MapData<LocalOrdinal,GlobalOrdinal>::checkIsDist() {
    bool global = false;
    if(comm_->getSize() > 1) {
      char localRep = 0;
      if (numGlobalEntries_ == numMyEntries_) {
        localRep = 1;
      }
      char allLocalRep;
      Teuchos::reduceAll<int>(*comm_,Teuchos::REDUCE_MIN,localRep,&allLocalRep);
      if (allLocalRep != 1) {
        global = true;
      }
    }
    return global;
  }

} // namespace Tpetra

#endif // TPETRA_MAPDATA_HPP

