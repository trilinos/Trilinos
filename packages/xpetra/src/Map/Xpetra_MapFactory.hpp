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
#ifndef XPETRA_MAPFACTORY_HPP
#define XPETRA_MAPFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMap.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMap.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

// This factory creates Xpetra::Map. User have to specify the exact class of object that he want to create (ie: a Xpetra::TpetraMap or a Xpetra::EpetraMap).

namespace Xpetra {

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MapFactory {

  private:
    //! Private constructor. This is a static class.
    MapFactory() {}

  public:

    //! Map constructor with Xpetra-defined contiguous uniform distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=Xpetra::GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return Teuchos::rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor transforming degrees of freedom
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map, LocalOrdinal numDofPerNode) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      LocalOrdinal N = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map->getNodeNumElements()*numDofPerNode);
      for (LocalOrdinal i = 0; i < N; i++)
        for (LocalOrdinal j = 0; j < numDofPerNode; j++)
          newElements[i*numDofPerNode + j] = oldElements[i]*numDofPerNode + j;
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
    }

    //! Create a locally replicated Map with the default node.
    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a locally replicated Map with a specified node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a uniform, contiguous Map with a user-specified node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode(UnderlyingLib lib, global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a uniform, contiguous Map with the default node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with the default node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createContigMap<LocalOrdinal,GlobalOrdinal>(numElements, localNumElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with a user-specified node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

  };

  template <>
  class MapFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  private:
    //! Private constructor. This is a static class.
    MapFactory() {}

  public:

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, indexBase, comm, lg, node) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::Details::getNode<Node> ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new EpetraMap(numGlobalElements, elementList, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_END;
    }

    //! Map constructor transforming degrees of freedom
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map, LocalOrdinal numDofPerNode) {
      XPETRA_MONITOR("MapFactory::Build");

      LocalOrdinal N = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map->getNodeNumElements()*numDofPerNode);
      for (LocalOrdinal i = 0; i < N; i++)
        for (LocalOrdinal j = 0; j < numDofPerNode; j++)
          newElements[i*numDofPerNode + j] = oldElements[i]*numDofPerNode + j;

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraMap(map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createLocalMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, node)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap((Xpetra::global_size_t)numElements, // num elements, global and local
                                            0,                                   // index base is zero
                                            comm, LocallyReplicated, node));
          return map.getConst();
        }
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode(UnderlyingLib lib, global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, comm, node)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap(numElements,        // num elements, global and local
                                            0,                  //index base is zero
                                            comm, GloballyDistributed, node));
          return map.getConst();
        }
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createUniformContigMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMap<int,int>(numElements, localNumElements, comm)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createContigMapWithNode(lib, numElements, localNumElements, comm, Kokkos::DefaultNode::getDefaultNode() );
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<int,int,Kokkos::DefaultNode::DefaultNodeType>(numElements, localNumElements, comm, node)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMap > map;
          map = Teuchos::rcp( new EpetraMap(numElements,localNumElements,
                                            0,  // index base is zero
                                            comm, node) );
          return map.getConst();
        }
#endif

      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_MAPFACTORY_SHORT
#endif
//TODO: removed unused methods
