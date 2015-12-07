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
#ifndef XPETRA_MAPFACTORY_HPP
#define XPETRA_MAPFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"

#ifdef HAVE_XPETRA_TPETRA
#  include "Xpetra_TpetraMap.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#  include "Xpetra_EpetraMap.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  /// \class MapFactory
  /// \brief Create an Xpetra::Map instance.
  ///
  /// Users must specify the exact class of the object that they want
  /// to create (either an Xpetra::TpetraMap or an Xpetra::EpetraMap).
  template <class LocalOrdinal = Map<>::local_ordinal_type,
            class GlobalOrdinal = typename Map<LocalOrdinal>::global_ordinal_type,
            class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
  class MapFactory {

  private:
    //! Private constructor. This is a static class.
    MapFactory() {}

  public:
    static Teuchos::RCP<Node> defaultArgNode() {
      // Workaround function for a deferred visual studio bug
      //
      // http://connect.microsoft.com/VisualStudio/feedback/details/719847/erroneous-error-c2783-could-not-deduce-template-argument
      //
      // Use this function for default arguments rather than calling
      // what is the return value below.  Also helps in reducing
      // duplication in various constructors.
      return KokkosClassic::Details::getNode<Node>();
    }

    //! Map constructor with Xpetra-defined contiguous uniform distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           GlobalOrdinal indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg = Xpetra::GloballyDistributed,
           const Teuchos::RCP<Node> &node = defaultArgNode ())
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return Teuchos::rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           GlobalOrdinal indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node> &node = defaultArgNode ())
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
           GlobalOrdinal indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node> &node = defaultArgNode ())
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor transforming degrees of freedom
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& map,
           LocalOrdinal numDofPerNode)
    {
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
    createLocalMap (UnderlyingLib lib,
                    size_t numElements,
                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a locally replicated Map with a specified node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode (UnderlyingLib lib,
                            size_t numElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                            const Teuchos::RCP< Node > &node)
    {
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
    createUniformContigMapWithNode (UnderlyingLib lib,
                                    global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                    const Teuchos::RCP< Node > &node)
    {
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
    createUniformContigMap (UnderlyingLib lib,
                            global_size_t numElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with the default node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap (UnderlyingLib lib,
                     global_size_t numElements,
                     size_t localNumElements,
                     const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal,Node>(Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with a user-specified node.
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode (UnderlyingLib lib,
                             global_size_t numElements,
                             size_t localNumElements,
                             const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                             const Teuchos::RCP< Node > &node)
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }
  };

  // TODO fix the guards (if Epetra and Tpetra enabled, ETI on etc...)
//#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (!defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))
  template <>
  class MapFactory<int, int, EpetraNode> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MapFactory() {}

  public:

    static Teuchos::RCP<Node> defaultArgNode () {
      // Workaround function for a deferred visual studio bug
      //
      // http://connect.microsoft.com/VisualStudio/feedback/details/719847/erroneous-error-c2783-could-not-deduce-template-argument
      //
      // Use this function for default arguments rather than calling
      // what is the return value below.  Also helps in reducing
      // duplication in various constructors.
      return KokkosClassic::Details::getNode<Node>();
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg=GloballyDistributed,
           const Teuchos::RCP<Node>& node = defaultArgNode ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      // TODO adapt to EpetraNode
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, indexBase, comm, lg, node) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node>& node = defaultArgNode())
    {
      XPETRA_MONITOR("MapFactory::Build");
#ifdef HAVE_XPETRA_TPETRA
      // adapt to node EpetraNode
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node>& node = defaultArgNode ()) {
      XPETRA_MONITOR("MapFactory::Build");
#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, elementList, indexBase, comm, node) );
#endif
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
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (map->lib() == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra) {
        Teuchos::RCP< EpetraMapT<int, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int, Node>((Xpetra::global_size_t)numElements, // num elements, global and local
                                                    0,                                   // index base is zero
                                                    comm, LocallyReplicated));
        return map.getConst();
       }
        //return MapFactory<int, GlobalOrdinal, Node>::createLocalMapWithNode(lib, numElements, comm, Teuchos::rcp(typename Node));
#endif
#endif

      XPETRA_FACTORY_END;
    }

    // TODO remove this
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Node> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node> (numElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<int, Node> > map;
          map = Teuchos::rcp( new EpetraMapT<int, Node>((Xpetra::global_size_t)numElements, // num elements, global and local
                                            0,                                   // index base is zero
                                            comm, LocallyReplicated, node));
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

    // TODO remove this
    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode (UnderlyingLib lib, global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                    const Teuchos::RCP<Node>& node)
    {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<int,GlobalOrdinal,Node> (numElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<int,Node> > map;
          map = Teuchos::rcp( new EpetraMapT<int,Node>(numElements,        // num elements, global and local
                                            0,                  //index base is zero
                                            comm, GloballyDistributed, node));
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      //if (lib == UseEpetra)
      //  return MapFactory<int, GlobalOrdinal, Node>::createUniformContigMapWithNode (lib, numElements, comm);
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<int,Node> > map;
          map = Teuchos::rcp( new EpetraMapT<int,Node>(numElements,        // num elements, global and local
                                            0,                  //index base is zero
                                            comm, GloballyDistributed));
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMap<int,GlobalOrdinal>(numElements, localNumElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm, defaultArgNode() );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Node> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_INT)))
  // specialization for Tpetra Map on EpetraNode and GO=int

      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<int,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<int, Node> > map;
          map = Teuchos::rcp( new EpetraMapT<int, Node>(numElements,localNumElements,
                                            0,  // index base is zero
                                            comm, node) );
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

  };
#endif

  // we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))
  template <>
  class MapFactory<int, long long, EpetraNode> {

    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    MapFactory() {}

  public:

    static Teuchos::RCP<Node> defaultArgNode() {
      // Workaround function for a deferred visual studio bug
      //
      // http://connect.microsoft.com/VisualStudio/feedback/details/719847/erroneous-error-c2783-could-not-deduce-template-argument
      //
      // Use this function for default arguments rather than calling
      // what is the return value below.  Also helps in reducing
      // duplication in various constructors.
      return KokkosClassic::Details::getNode<Node>();
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg=GloballyDistributed,
           const Teuchos::RCP<Node>& node = defaultArgNode ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
  // specialization for Tpetra Map on EpetraNode and GO=int

      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, indexBase, comm, lg, node) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> >
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<Node>& node = defaultArgNode())
    {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node>& node = defaultArgNode ()) {
      XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, elementList, indexBase, comm, node) );
#endif
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
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (map->lib() == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm(), map->getNode()) );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createLocalMapWithNode (lib, numElements, comm,
                                                             defaultArgNode());
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Node> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node> (numElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<long long, Node> > map;
          map = Teuchos::rcp( new EpetraMapT<long long, Node>((Xpetra::global_size_t)numElements, // num elements, global and local
                                            0,                                   // index base is zero
                                            comm, LocallyReplicated, node));
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode (UnderlyingLib lib, global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                    const Teuchos::RCP<Node>& node)
    {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<int,GlobalOrdinal,Node> (numElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<long long, Node> > map;
          map = Teuchos::rcp( new EpetraMapT<long long, Node>(numElements,        // num elements, global and local
                                            0,                  //index base is zero
                                            comm, GloballyDistributed, node));
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createUniformContigMapWithNode (lib, numElements, comm,
                                                                           defaultArgNode());
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMap<int,GlobalOrdinal>(numElements, localNumElements, comm)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm, defaultArgNode() );
#endif
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Node> &node) {
       XPETRA_MONITOR("MapFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(HAVE_TPETRA_INST_SERIAL)) && (defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<int,GlobalOrdinal,Node>(numElements, localNumElements, comm, node)));
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
      if (lib == UseEpetra)
        {

          Teuchos::RCP< EpetraMapT<long long, Node> > map;
          map = Teuchos::rcp( new EpetraMapT<long long, Node>(numElements,localNumElements,
                                            0,  // index base is zero
                                            comm, node) );
          return map.getConst();
        }
#endif
#endif

      XPETRA_FACTORY_END;
    }

  };
#endif

}

#define XPETRA_MAPFACTORY_SHORT
#endif
//TODO: removed unused methods
