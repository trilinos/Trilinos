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

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAPFACTORY_DEF_HPP
#define XPETRA_STRIDEDMAPFACTORY_DEF_HPP

#include "Xpetra_StridedMapFactory_decl.hpp"

#include "Xpetra_Exceptions.hpp"

// This factory creates Xpetra::Map. User have to specify the exact class of
// object that he want to create (ie: a Xpetra::TpetraMap or a Xpetra::EpetraMap).



namespace Xpetra {



template<class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
StridedMapFactory()
{
}


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
template<class LocalOrdinal, class GlobalOrdinal, class Node>
TPETRA_DEPRECATED
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                 lib,
      global_size_t                                 numGlobalElements,
      GlobalOrdinal                                 indexBase,
      std::vector<size_t>&                          stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
      LocalOrdinal                                  stridedBlockId,
      GlobalOrdinal                                 offset,
      LocalGlobal                                   lg,
      const Teuchos::RCP<Node>&                     /* node */)
{
    return Build(lib, numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg);
}
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                 lib,
      global_size_t                                 numGlobalElements,
      GlobalOrdinal                                 indexBase,
      std::vector<size_t>&                          stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
      LocalOrdinal                                  stridedBlockId,
      GlobalOrdinal                                 offset,
      LocalGlobal                                   lg)
{
    return rcp(new Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>(lib, numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg));
}


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
template<class LocalOrdinal, class GlobalOrdinal, class Node>
TPETRA_DEPRECATED
RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                 lib,
      global_size_t                                 numGlobalElements,
      size_t                                        numLocalElements,
      GlobalOrdinal                                 indexBase,
      std::vector<size_t>&                          stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
      LocalOrdinal                                  stridedBlockId,
      GlobalOrdinal                                 offset,
      const Teuchos::RCP<Node>&                     /* node */)
{
    return Build(lib, numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset);
}
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                 lib,
      global_size_t                                 numGlobalElements,
      size_t                                        numLocalElements,
      GlobalOrdinal                                 indexBase,
      std::vector<size_t>&                          stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
      LocalOrdinal                                  stridedBlockId,
      GlobalOrdinal                                 offset)
{
    return rcp(new StridedMap(lib, numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset));
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(const RCP<const Map>& map,
      std::vector<size_t>&  stridingInfo,
      LocalOrdinal          stridedBlockId,
      GlobalOrdinal         offset)
{
    return rcp(new StridedMap(map, stridingInfo, map->getIndexBase(), stridedBlockId, offset));
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(const RCP<const StridedMap>& map, LocalOrdinal stridedBlockId)
{
    TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId < 0, 
                               Exceptions::RuntimeError, 
                               "Xpetra::StridedMapFactory::Build: constructor expects stridedBlockId > -1.");

    TEUCHOS_TEST_FOR_EXCEPTION(map->getStridedBlockId() != -1,
                               Exceptions::RuntimeError,
                               "Xpetra::StridedMapFactory::Build: constructor expects a full map (stridedBlockId == -1).");

    std::vector<size_t> stridingInfo = map->getStridingData();

    Teuchos::ArrayView<const GlobalOrdinal> dofGids = map->getNodeElementList();
    // std::sort(dofGids.begin(),dofGids.end()); // TODO: do we need this?

    // determine nStridedOffset
    size_t nStridedOffset = 0;
    for(int j = 0; j < map->getStridedBlockId(); j++) 
    { 
        nStridedOffset += stridingInfo[ j ]; 
    }

    size_t numMyBlockDofs = (stridingInfo[ stridedBlockId ] * map->getNodeNumElements()) / map->getFixedBlockSize();

    std::vector<GlobalOrdinal> subBlockDofGids(numMyBlockDofs);

    // TODO fill vector with dofs
    LocalOrdinal ind = 0;
    for(typename Teuchos::ArrayView<const GlobalOrdinal>::iterator it = dofGids.begin(); it != dofGids.end(); ++it)
    {
        if(map->GID2StridingBlockId(*it) == Teuchos::as<size_t>(stridedBlockId))
        {
            subBlockDofGids[ ind++ ] = *it;
        }
    }

    const Teuchos::ArrayView<const GlobalOrdinal> subBlockDofGids_view(&subBlockDofGids[ 0 ], subBlockDofGids.size());

    return rcp(new StridedMap(map->lib(),
                              Teuchos::OrdinalTraits<global_size_t>::invalid(),
                              subBlockDofGids_view,
                              map->getIndexBase(),
                              stridingInfo,
                              map->getComm(),
                              stridedBlockId));
}


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(const StridedMap& map)
{
    XPETRA_MONITOR("MapFactory::Build");

    LocalOrdinal                            N           = map.getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> oldElements = map.getNodeElementList();
    Teuchos::Array<GlobalOrdinal>           newElements(map.getNodeNumElements());

    for(LocalOrdinal i = 0; i < N; i++) 
    { 
        newElements[ i ] = oldElements[ i ]; 
    }

    std::vector<size_t> strData = map.getStridingData();
    return rcp(new StridedMap(map.lib(), map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId()));

    // XPETRA_FACTORY_END;
}


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
template<class LocalOrdinal, class GlobalOrdinal, class Node>
TPETRA_DEPRECATED
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                  lib,
      global_size_t                                  numGlobalElements,
      const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
      GlobalOrdinal                                  indexBase,
      std::vector<size_t>&                           stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
      LocalOrdinal                                   stridedBlockId,      // FIXME (mfh 03 Sep 2014) This breaks if LocalOrdinal is unsigned
      GlobalOrdinal                                  offset,
      const Teuchos::RCP<Node>&                      /* node */)
{
    return Build(lib, numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId, offset);
}
#endif      // TPETRA_ENABLE_DEPRECATED_CODE


template<class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node>>
StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::
Build(UnderlyingLib                                  lib,
      global_size_t                                  numGlobalElements,
      const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
      GlobalOrdinal                                  indexBase,
      std::vector<size_t>&                           stridingInfo,
      const Teuchos::RCP<const Teuchos::Comm<int>>&  comm,
      LocalOrdinal                                   stridedBlockId,      // FIXME (mfh 03 Sep 2014) This breaks if LocalOrdinal is unsigned
      GlobalOrdinal                                  /* offset */)
{
    return rcp(new StridedMap(lib, numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId));
}



}      // namespace Xpetra



#endif      // __XPETRA_STRIDEDMAPFACTORY_DEF_HPP__

// TODO: removed unused methods
