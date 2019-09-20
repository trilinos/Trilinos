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
#include "Xpetra_MapFactory.hpp"

#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_TpetraMap.hpp"

namespace Xpetra {

#if defined(HAVE_XPETRA_EPETRA)


#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)


    MapFactory<int, int, EpetraNode>::
    MapFactory()
    {
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t                                 numGlobalElements,
           int                                           indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal                                   lg,
           const Teuchos::RCP<EpetraNode>& /* node */ )
    {
      return Build(lib, numGlobalElements, indexBase, comm, lg);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, int, EpetraNode>>
    MapFactory<int, int, EpetraNode>::
    Build (UnderlyingLib                                 lib,
           global_size_t                                 numGlobalElements,
           int                                           indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal                                   lg)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg) );
#endif

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, indexBase, comm, lg) );

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<EpetraNode>& /* node */ )
    {
      return Build(lib, numGlobalElements, numLocalElements, indexBase, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm) );
#endif

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, numLocalElements, indexBase, comm) );

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const int> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          const Teuchos::RCP<EpetraNode>& /* node */ )
    {
      return Build(lib, numGlobalElements, elementList, indexBase, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const int> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    {
      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

      XPETRA_MONITOR("MapFactory::Build");
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm) );
#endif  // HAVE_XPETRA_TPETRA

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<int, Node>(numGlobalElements, elementList, indexBase, comm) );

      XPETRA_FACTORY_END;
    }


    //! Map constructor transforming degrees of freedom
    //! for numDofPerNode this acts like a deep copy
    Teuchos::RCP<Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    Build(const Teuchos::RCP<const Map<int, int, EpetraNode> >& map,
          int                                                   numDofPerNode)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null())
      {
        TEUCHOS_TEST_FOR_EXCEPTION(numDofPerNode!=1, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::Build: When provided a BlockedMap numDofPerNode must set to be one. It is set to " << numDofPerNode << ".");
        return rcp(new BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(*bmap));
      }

      LocalOrdinal N = Teuchos::as<LocalOrdinal>(map->getNodeNumElements());
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map->getNodeNumElements()*numDofPerNode);
      for (LocalOrdinal i = 0; i < N; i++)
      {
        for (LocalOrdinal j = 0; j < numDofPerNode; j++)
        {
          newElements[i*numDofPerNode + j] = oldElements[i]*numDofPerNode + j;
        }
      }

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
      {
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm()) );
      }
#endif // HAVE_XPETRA_TPETRA

      if (map->lib() == UseEpetra)
      {
        return rcp( new EpetraMapT<int, Node>(map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm()) );
      }

      XPETRA_FACTORY_END;
    }


    Teuchos::RCP<const Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    createLocalMap(UnderlyingLib                                    lib,
                   size_t                                           numElements,
                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createLocalMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

      if (lib == UseEpetra) {
        Teuchos::RCP< EpetraMapT<int, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int, Node>( (Xpetra::global_size_t)numElements, // num elements, global and local
                                                       0,                                   // index base is zero
                                                       comm, LocallyReplicated)
                          );
        return map.getConst();
      }

      XPETRA_FACTORY_END;
    }


    // TODO remove this
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib                                    lib,
                           size_t                                           numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                           const Teuchos::RCP<Node> & /* node */)
    {
      return createLocalMapWithNode(lib, numElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE



    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node> (numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createLocalMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

      if (lib == UseEpetra) {
        Teuchos::RCP< EpetraMapT<int, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int, Node>((Xpetra::global_size_t)numElements, // num elements, global and local
                                          0,                                   // index base is zero
                                          comm, LocallyReplicated));
        return map.getConst();
      }

      XPETRA_FACTORY_END;
    }



    // TODO remove this
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createUniformContigMapWithNode (UnderlyingLib                                    lib,
                                    global_size_t                                    numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                    const Teuchos::RCP<Node>& /* node */)
    {
      return createUniformContigMapWithNode(lib, numElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE



    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createUniformContigMapWithNode (UnderlyingLib lib,
                                    global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<int,GlobalOrdinal,Node> (numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

      if (lib == UseEpetra) {
        Teuchos::RCP< EpetraMapT<int,Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int,Node>(numElements,        // num elements, global and local
            0,                  //index base is zero
            comm, GloballyDistributed));
        return map.getConst();
      }

      XPETRA_FACTORY_END;
    }


    Teuchos::RCP< const Map<int, int, EpetraNode> >
    MapFactory<int, int, EpetraNode>::
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
         return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#else
         TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
           "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif  // HAVE_XPETRA_TPETRA

      if (lib == UseEpetra)
      {
        Teuchos::RCP< EpetraMapT<int,Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int,Node>(numElements,        // num elements, global and local
            0,                  //index base is zero
            comm, GloballyDistributed));
        return map.getConst();
      }
      XPETRA_FACTORY_END;
    }


    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createContigMap(UnderlyingLib lib,
                    global_size_t numElements,
                    size_t localNumElements,
                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

       #ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
           #if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
                (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
               return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<int,GlobalOrdinal,Node>(numElements, localNumElements, comm)));
           #else
               TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                    "Xpetra::MapFactory::createContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
            #endif
        #endif

      if (lib == UseEpetra)
      {
        return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm);
      }

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                            const Teuchos::RCP<Node> & /* node */)
    {
      return createContigMapWithNode(lib, numElements, localNumElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    Teuchos::RCP< const Map<int, int, EpetraNode>  >
    MapFactory<int, int, EpetraNode>::
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = int;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createContigMapWithNode<int,GlobalOrdinal,Node>(numElements, localNumElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=int");
#endif
#endif

      if (lib == UseEpetra)
      {
        Teuchos::RCP< EpetraMapT<int, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<int, Node>(numElements,localNumElements,
            0,  // index base is zero
            comm) );
        return map.getConst();
      }
      XPETRA_FACTORY_END;
    }



#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)








// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)



    //! Private constructor. This is a static class.
    MapFactory<int, long long, EpetraNode>::
    MapFactory()
    {
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal lg,
           const Teuchos::RCP<Node>& /* node */)
    {
      return Build(lib, numGlobalElements, indexBase, comm, lg);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, long long, EpetraNode>>
    MapFactory<int, long long, EpetraNode>::
    Build (UnderlyingLib                                 lib,
           global_size_t                                 numGlobalElements,
           int                                           indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           LocalGlobal                                   lg)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg) );
#endif

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, indexBase, comm, lg) );

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
           const Teuchos::RCP<EpetraNode>& /* node */)
    {
      return Build(lib, numGlobalElements, numLocalElements, indexBase, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    Build (UnderlyingLib lib,
           global_size_t numGlobalElements,
           size_t numLocalElements,
           int indexBase,
           const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm) );
#endif

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, numLocalElements, indexBase, comm) );

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    RCP<Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const long long> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
          const Teuchos::RCP<EpetraNode>& /* node */)
    {
      return Build(lib, numGlobalElements, elementList, indexBase, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    RCP<Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    Build(UnderlyingLib lib,
          global_size_t numGlobalElements,
          const Teuchos::ArrayView<const long long> &elementList,
          int indexBase,
          const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm) );
#endif

      if (lib == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(numGlobalElements, elementList, indexBase, comm) );

      XPETRA_FACTORY_END;
    }

    //! Map constructor transforming degrees of freedom
    //! for numDofPerNode this acts like a deep copy
     Teuchos::RCP<Map<int, long long, EpetraNode> >
     MapFactory<int, long long, EpetraNode>::
     Build(const Teuchos::RCP<const Map<int, long long, EpetraNode> >& map,
           int numDofPerNode)
    {
      XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

      RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> > bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node> >(map);
      if(!bmap.is_null()) {
        TEUCHOS_TEST_FOR_EXCEPTION(numDofPerNode!=1, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::Build: When provided a BlockedMap numDofPerNode must set to be one. It is set to " << numDofPerNode << ".");
        return rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(*bmap));
      }

      LocalOrdinal N = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> oldElements = map->getNodeElementList();
      Teuchos::Array<GlobalOrdinal> newElements(map->getNodeNumElements()*numDofPerNode);
      for (LocalOrdinal i = 0; i < N; i++)
        for (LocalOrdinal j = 0; j < numDofPerNode; j++)
          newElements[i*numDofPerNode + j] = oldElements[i]*numDofPerNode + j;

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm()) );
#endif

      if (map->lib() == UseEpetra)
        return rcp( new EpetraMapT<long long, Node>(map->getGlobalNumElements()*numDofPerNode, newElements, map->getIndexBase(), map->getComm()) );

      XPETRA_FACTORY_END;
    }


    Teuchos::RCP<const Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    createLocalMap(UnderlyingLib lib,
                   size_t numElements,
                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
        "Xpetra::MapFactory::createLocalMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createLocalMapWithNode (lib, numElements, comm);

      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, long long, EpetraNode>  >
    MapFactory<int, long long, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                           const Teuchos::RCP<EpetraNode> & /* node */)
    {
      return createLocalMapWithNode(lib, numElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    Teuchos::RCP< const Map<int, long long, EpetraNode>  >
    MapFactory<int, long long, EpetraNode>::
    createLocalMapWithNode(UnderlyingLib lib,
                           size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<int, GlobalOrdinal, Node> (numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createLocalMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

      if (lib == UseEpetra)
      {
        Teuchos::RCP< EpetraMapT<long long, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<long long, Node>((Xpetra::global_size_t)numElements, // num elements, global and local
            0,                                   // index base is zero
            comm, LocallyReplicated));
        return map.getConst();
      }
      XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, long long, EpetraNode>  >
    MapFactory<int, long long, EpetraNode>::
    createUniformContigMapWithNode (UnderlyingLib lib,
                                    global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                                    const Teuchos::RCP<EpetraNode>& /* node */)
    {
      return createUniformContigMapWithNode(lib, numElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    Teuchos::RCP< const Map<int, long long, EpetraNode>  >
    MapFactory<int, long long, EpetraNode>::
    createUniformContigMapWithNode (UnderlyingLib lib,
                                    global_size_t numElements,
                                    const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp (new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<int,GlobalOrdinal,Node> (numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createUniformContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

      if (lib == UseEpetra)
      {
        Teuchos::RCP< EpetraMapT<long long, Node> > map;
        map = Teuchos::rcp( new EpetraMapT<long long, Node>(numElements,        // num elements, global and local
            0,                  //index base is zero
            comm, GloballyDistributed));
        return map.getConst();
      }
      XPETRA_FACTORY_END;
    }


    Teuchos::RCP< const Map<int, long long, EpetraNode> >
    MapFactory<int, long long, EpetraNode>::
    createUniformContigMap(UnderlyingLib lib,
                           global_size_t numElements,
                           const Teuchos::RCP< const Teuchos::Comm< int > > &comm)
    {
       XPETRA_MONITOR("MapFactory::Build");

      using LocalOrdinal  = int;
      using GlobalOrdinal = long long;
      using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
       if (lib == UseTpetra)
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal, Node>(numElements, comm)));
#else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
          "Xpetra::MapFactory::createUniformContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

      if (lib == UseEpetra)
        return MapFactory<int, GlobalOrdinal, Node>::createUniformContigMapWithNode (lib, numElements, comm);

      XPETRA_FACTORY_END;
    }


    Teuchos::RCP<const Map<int, long long, EpetraNode>>
    MapFactory<int, long long, EpetraNode>::createContigMap(UnderlyingLib                                 lib,
                                                            global_size_t                                 numElements,
                                                            size_t                                        localNumElements,
                                                            const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
    {
        XPETRA_MONITOR("MapFactory::Build");

        using LocalOrdinal  = int;
        using GlobalOrdinal = long long;
        using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
        if(lib == UseTpetra)
#if((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) \
    || (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
            return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(
              Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
            TEUCHOS_TEST_FOR_EXCEPTION(true,
                                       Xpetra::Exceptions::RuntimeError,
                                       "Xpetra::MapFactory::createContigMap: Cannot create Xpetra::TpetraMap, since Tpetra is not instantiated on "
                                       "EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif

        if(lib == UseEpetra)
            return MapFactory<int, GlobalOrdinal, Node>::createContigMapWithNode(lib, numElements, localNumElements, comm);

        XPETRA_FACTORY_END;
    }


#ifdef TPETRA_ENABLE_DEPRECATED_CODE
    TPETRA_DEPRECATED
    Teuchos::RCP< const Map<int, long long, EpetraNode>  >
    MapFactory<int, long long, EpetraNode>::
    createContigMapWithNode(UnderlyingLib lib,
                            global_size_t numElements,
                            size_t localNumElements,
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm,
                            const Teuchos::RCP<EpetraNode> & /* node */ )
    {
      return createContigMapWithNode(lib, numElements, localNumElements, comm);
    }
#endif // TPETRA_ENABLE_DEPRECATED_CODE


    Teuchos::RCP<const Map<int, long long, EpetraNode>>
    MapFactory<int, long long, EpetraNode>::
    createContigMapWithNode(UnderlyingLib                                 lib,
                            global_size_t                                 numElements,
                            size_t                                        localNumElements,
                            const Teuchos::RCP<const Teuchos::Comm<int>>& comm)
    {
        XPETRA_MONITOR("MapFactory::Build");

        using LocalOrdinal  = int;
        using GlobalOrdinal = long long;
        using Node          = EpetraNode;

#ifdef HAVE_XPETRA_TPETRA
        if(lib == UseTpetra)
#if((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) \
    || (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
            return rcp(new TpetraMap<LocalOrdinal, GlobalOrdinal, Node>(
              Tpetra::createContigMapWithNode<int, GlobalOrdinal, Node>(numElements, localNumElements, comm)));
#else
            TEUCHOS_TEST_FOR_EXCEPTION(true,
                                       Xpetra::Exceptions::RuntimeError,
                                       "Xpetra::MapFactory::createContigMapWithNode: Cannot create Xpetra::TpetraMap, since Tpetra is not "
                                       "instantiated on EpetraNode (Serial or OpenMP, depending on configuration) and/or GO=long long");
#endif
#endif      // HAVE_XPETRA_TPETRA

        if(lib == UseEpetra)
        {
            Teuchos::RCP<EpetraMapT<long long, Node>> map;
            map = Teuchos::rcp(new EpetraMapT<long long, Node>(numElements,
                                                               localNumElements,
                                                               0,      // index base is zero
                                                               comm));
            return map.getConst();
        }

        XPETRA_FACTORY_END;
    }



#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)


#endif  // #if defined(HAVE_XPETRA_EPETRA)



}   // namespace Xpetra

// EOF
