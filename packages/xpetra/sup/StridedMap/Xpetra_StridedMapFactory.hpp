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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAPFACTORY_HPP
#define XPETRA_STRIDEDMAPFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_StridedMap.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_StridedTpetraMap.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_StridedEpetraMap.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

// This factory creates Xpetra::Map. User have to specify the exact class of object that he want to create (ie: a Xpetra::TpetraMap or a Xpetra::EpetraMap).

namespace Xpetra {

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class StridedMapFactory {

  private:
    //! Private constructor. This is a static class.
    StridedMapFactory() {}

  public:

    //! Map constructor with Xpetra-defined contiguous uniform distribution.
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, LocalGlobal lg=Xpetra::GloballyDistributed, const Teuchos::RCP<Node> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return Teuchos::rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution.
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, const Teuchos::RCP<Node> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) {
        Teuchos::RCP<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> > tmap = Teuchos::rcp_dynamic_cast<const TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
        TEUCHOS_TEST_FOR_EXCEPTION(tmap == Teuchos::null, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: bad cast. map is not a TpetraMap.");
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (tmap->getTpetra_Map(), stridingInfo, stridedBlockId, offset) );
      }
#endif
      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    // special constructor for generating a given subblock of a strided map
    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(const RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >& map, LocalOrdinal stridedBlockId) {
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId < 0, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: constructor expects stridedBlockId > -1.");
#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra) {
        TEUCHOS_TEST_FOR_EXCEPTION(map->getStridedBlockId() != -1, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: constructor expects a full map (stridedBlockId == -1).");
        std::vector<size_t> stridingInfo = map->getStridingData();

        /////////////////////////////////////////////
        Teuchos::ArrayView< const GlobalOrdinal > dofGids = map->getNodeElementList();
        //std::sort(dofGids.begin(),dofGids.end()); // TODO: do i need this?

        // determine nStridedOffset
        size_t nStridedOffset = 0;
        for(int j=0; j<map->getStridedBlockId(); j++) {
          nStridedOffset += stridingInfo[j];
        }

        size_t numMyBlockDofs = (stridingInfo[stridedBlockId] * map->getNodeNumElements()) / map->getFixedBlockSize();
        std::vector<GlobalOrdinal> subBlockDofGids(numMyBlockDofs);

        // TODO fill vector with dofs
        typename Teuchos::ArrayView< const GlobalOrdinal >::iterator it;
        for(it = dofGids.begin(); it!=dofGids.end(); ++it) {
          if(map->GID2StridingBlockId( *it ) == stridedBlockId) {
            subBlockDofGids.push_back( *it );
          }
        }

        const Teuchos::ArrayView<const LocalOrdinal> subBlockDofGids_view(&subBlockDofGids[0],subBlockDofGids.size());

        // call constructor for TpetraMap
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (subBlockDofGids.size(), subBlockDofGids_view, map->getIndexBase(), stridingInfo, map->getComm(), stridedBlockId, map->getNode()) );
        ////////////////////////////////////////

      }
#endif
      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
    }

    //! Create copy of existing map (this just creates a copy of your map, it's not a clone in the sense of Tpetra)
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(const StridedMap<LocalOrdinal,GlobalOrdinal, Node>& map) {
      XPETRA_MONITOR("MapFactory::Build");

      LocalOrdinal N = map.getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map.getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map.getNodeNumElements());
      for (LocalOrdinal i = 0; i < N; i++)
          newElements[i] = oldElements[i];

#ifdef HAVE_XPETRA_TPETRA
      if (map.lib() == UseTpetra) {
        std::vector<size_t> strData = map.getStridingData();
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId(), map.getOffset(), map.getNode()) );
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map.lib() == UseEpetra) {
        std::vector<size_t> strData = map.getStridingData();
        return rcp( new StridedEpetraMap (map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId(), map.getOffset(), map.getNode()) );
      }

#endif

      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution. (for experts only. There is no special check whether the generated strided maps are valid)
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, const Teuchos::RCP<Node> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

  };

  template <>
  class StridedMapFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

  private:
    //! Private constructor. This is a static class.
    StridedMapFactory() {}

  public:

    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, int indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<KokkosClassic::DefaultNode::DefaultNodeType> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, lg, node) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, int indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, const Teuchos::RCP<KokkosClassic::DefaultNode::DefaultNodeType> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& map, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) {
	Teuchos::RCP<const TpetraMap<LocalOrdinal,GlobalOrdinal,Node> > tmap = Teuchos::rcp_dynamic_cast<const TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
	TEUCHOS_TEST_FOR_EXCEPTION(tmap == Teuchos::null, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: bad cast. map is not a TpetraMap.");
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (tmap->getTpetra_Map(), stridingInfo, stridedBlockId, offset) );
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra) {
	Teuchos::RCP<const EpetraMap> emap = Teuchos::rcp_dynamic_cast<const EpetraMap>(map);
	TEUCHOS_TEST_FOR_EXCEPTION(emap == Teuchos::null, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: bad cast. map is not a EpetraMap.");
        return rcp( new StridedEpetraMap(Teuchos::rcpFromRef(emap->getEpetra_Map()), stridingInfo,stridedBlockId, offset) ); // ugly: EpetraMap and TpetraMap differ here
      }
#endif

      XPETRA_FACTORY_END;
    }

    // special constructor for generating a given subblock of a strided map
    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(const RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >& map, LocalOrdinal stridedBlockId) {
      TEUCHOS_TEST_FOR_EXCEPTION(stridedBlockId < 0, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: constructor expects stridedBlockId > -1.");
      //typename Teuchos::ArrayView< const GlobalOrdinal >::iterator it;
#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra) {
        TEUCHOS_TEST_FOR_EXCEPTION(map->getStridedBlockId() != -1, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: constructor expects a full map (stridedBlockId == -1).");
        std::vector<size_t> stridingInfo = map->getStridingData();

        /////////////////////////////////////////////
        Teuchos::ArrayView< const GlobalOrdinal > dofGids = map->getNodeElementList();
        //std::sort(dofGids.begin(),dofGids.end()); // TODO: do i need this?

        // determine nStridedOffset
        size_t nStridedOffset = 0;
        for(int j=0; j<map->getStridedBlockId(); j++) {
          nStridedOffset += stridingInfo[j];
        }

        size_t numMyBlockDofs = stridingInfo[stridedBlockId] / map->getFixedBlockSize() * map->getNodeNumElements();
        std::vector<GlobalOrdinal> subBlockDofGids(numMyBlockDofs);

        GlobalOrdinal offset = map->getOffset();

        // fill vector with dofs

        for(Teuchos::ArrayView< const GlobalOrdinal >::iterator it = dofGids.begin(); it!=dofGids.end(); ++it) {
          if(map->GID2StridingBlockId( *it ) == Teuchos::as<size_t>(stridedBlockId)) {
            subBlockDofGids.push_back( *it );
          }
        }

        const Teuchos::ArrayView<const LocalOrdinal> subBlockDofGids_view(&subBlockDofGids[0],subBlockDofGids.size());

        // call constructor for TpetraMap
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (/*subBlockDofGids.size()*/Teuchos::OrdinalTraits<global_size_t>::invalid(), subBlockDofGids_view, map->getIndexBase(), stridingInfo, map->getComm(), stridedBlockId, offset, map->getNode()) );
        ////////////////////////////////////////

      }
#endif
#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra) {
        TEUCHOS_TEST_FOR_EXCEPTION(map->getStridedBlockId() != -1, Exceptions::RuntimeError,"Xpetra::StridedMapFactory::Build: constructor expects a full map (stridedBlockId == -1).");
        std::vector<size_t> stridingInfo = map->getStridingData();

        /////////////////////////////////////////////
        Teuchos::ArrayView< const GlobalOrdinal > dofGids = map->getNodeElementList();
        //std::sort(dofGids.begin(),dofGids.end()); // TODO: do i need this?

        // determine nStridedOffset
        size_t nStridedOffset = 0;
        for(int j=0; j<map->getStridedBlockId(); j++) {
          nStridedOffset += stridingInfo[j];
        }

        size_t numMyBlockDofs = stridingInfo[stridedBlockId] / map->getFixedBlockSize() * map->getNodeNumElements();
        std::vector<GlobalOrdinal> subBlockDofGids(numMyBlockDofs);

        GlobalOrdinal offset = map->getOffset();

        // TODO fill vector with dofs
        //Teuchos::ArrayView< const GlobalOrdinal >::iterator it;
        for(Teuchos::ArrayView< const GlobalOrdinal >::iterator it = dofGids.begin(); it!=dofGids.end(); ++it) {
          if(map->GID2StridingBlockId( *it ) == Teuchos::as<size_t>(stridedBlockId)) {
            subBlockDofGids.push_back( *it );
          }
        }

        const Teuchos::ArrayView<const LocalOrdinal> subBlockDofGids_view(&subBlockDofGids[0],subBlockDofGids.size());

        // call constructor for TpetraMap
        return rcp( new StridedEpetraMap(Teuchos::OrdinalTraits<global_size_t>::invalid(), subBlockDofGids_view, map->getIndexBase(), stridingInfo, map->getComm(), stridedBlockId, offset, map->getNode()) );
        ////////////////////////////////////////

      }
#endif
      XPETRA_FACTORY_END;
    }

    //! Create copy of existing map (this just creates a copy of your map, it's not a clone in the sense of Tpetra)
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(const StridedMap<LocalOrdinal,GlobalOrdinal, Node>& map) {
      XPETRA_MONITOR("MapFactory::Build");

      LocalOrdinal N = map.getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map.getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map.getNodeNumElements());
      for (LocalOrdinal i = 0; i < N; i++)
          newElements[i] = oldElements[i];

#ifdef HAVE_XPETRA_TPETRA
      if (map.lib() == UseTpetra) {
        std::vector<size_t> strData = map.getStridingData();
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId(), map.getOffset(), map.getNode()) );
      }
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map.lib() == UseEpetra) {
        std::vector<size_t> strData = map.getStridingData();
        return rcp( new StridedEpetraMap (map.getGlobalNumElements(), newElements, map.getIndexBase(), strData, map.getComm(), map.getStridedBlockId(), map.getOffset(), map.getNode()) );
      }

#endif

      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution. (for experts only. There is no special check whether the generated strided maps are valid)
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, GlobalOrdinal offset = 0, const Teuchos::RCP<Node> &node = KokkosClassic::Details::getNode<Node> ()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, elementList, indexBase, stridingInfo, comm, stridedBlockId, offset, node) );
#endif
      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_STRIDEDMAPFACTORY_SHORT
#endif
//TODO: removed unused methods
