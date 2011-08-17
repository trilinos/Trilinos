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
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=Xpetra::GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return Teuchos::rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, comm, lg, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)       
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }
        
    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    static Teuchos::RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a locally replicated Map with the default node.
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createLocalMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a locally replicated Map with a specified node.
    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a uniform, contiguous Map with a user-specified node.
    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createUniformContigMapWithNode(UnderlyingLib lib, global_size_t numElements,
                                   const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp(new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements, comm, node)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a uniform, contiguous Map with the default node.
    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createUniformContigMap(UnderlyingLib lib, global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createUniformContigMap<LocalOrdinal,GlobalOrdinal>(numElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with the default node.
    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp(new Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal>(Tpetra::createContigMap<LocalOrdinal,GlobalOrdinal>(numElements, localNumElements, comm)));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Create a (potentially) non-uniform, contiguous Map with a user-specified node.
    Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMapWithNode(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, 
                            const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
 
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
    
    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {

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

    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {

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
        
    static RCP<Map<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {

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

    static Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal, Node> >
    createLocalMap(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createLocalMap<int,int>(numElements, comm)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createLocalMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createLocalMapWithNode(UnderlyingLib lib, size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node) {
 
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
 
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp( new TpetraMap<LocalOrdinal,GlobalOrdinal, Node> (Tpetra::createUniformContigMap<int,int>(numElements, comm)));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return MapFactory<int, int, Kokkos::DefaultNode::DefaultNodeType>::createUniformContigMapWithNode(lib, numElements, comm, Kokkos::DefaultNode::getDefaultNode());
#endif

      XPETRA_FACTORY_END;
    }

    static Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal, Node>  >
    createContigMap(UnderlyingLib lib, global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
 
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
