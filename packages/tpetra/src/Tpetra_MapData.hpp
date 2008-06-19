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

  template<typename OrdinalType>
  MapData<OrdinalType>::MapData(
            OrdinalType indexBase, 
            OrdinalType numGlobalEntries,
            OrdinalType numMyEntries,
            OrdinalType minAllGID,
            OrdinalType maxAllGID,
            OrdinalType minMyGID,
            OrdinalType maxMyGID,
            const map<OrdinalType, OrdinalType>& lgMap,
            const map<OrdinalType, OrdinalType>& glMap,
            bool contiguous,
            Teuchos::RCP< Platform<OrdinalType, OrdinalType> > platform,
            Teuchos::RCP< Comm<OrdinalType, OrdinalType> > comm)
      : Teuchos::Object("Tpetra::MapData")
      , Platform_(platform)
      , Comm_(comm)
      , numGlobalEntries(numGlobalEntries)
      , numMyEntries(numMyEntries)
      , indexBase_(indexBase)
      , minLID_(Teuchos::OrdinalTraits<OrdinalType>::zero())
      , maxLID_(minLID_ + numMyEntries_ - Teuchos::OrdinalTraits<OrdinalType>::one())
      , minMyGID_(minMyGID)
      , maxMyGID_(maxMyGID)
      , minAllGID_(minAllGID)
      , maxAllGID_(maxAllGID)
      , contiguous_(contiguous)
      , global_(checkGlobalness())
      , haveDirectory_(false)
      , lgMap_(lgMap)
      , glMap_(glMap)
      , myGlobalEntries_()
      , Directory_()
    {}

  template<typename OrdinalType>
  MapData<OrdinalType>::~MapData() {}

} // namespace Tpetra

#endif // TPETRA_MAPDATA_HPP

