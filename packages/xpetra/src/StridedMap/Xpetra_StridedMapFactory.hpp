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

  template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class StridedMapFactory {
    
  private:
    //! Private constructor. This is a static class. 
    StridedMapFactory() {}
    
  public:
    
    //! Map constructor with Xpetra-defined contiguous uniform distribution.
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, LocalGlobal lg=Xpetra::GloballyDistributed, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return Teuchos::rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, lg, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map constructor with a user-defined contiguous distribution.
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, GlobalOrdinal indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId=-1, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)       
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

#if 0  // TODO
    //! Map constructor with user-defined non-contiguous (arbitrary) distribution.
    static Teuchos::RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }
#endif


  };

  template <>
  class StridedMapFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    
  private:
    //! Private constructor. This is a static class. 
    StridedMapFactory() {}
    
  public:
    
    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, int indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId, LocalGlobal lg=GloballyDistributed, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, lg, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, indexBase, stridingInfo, comm, stridedBlockId, lg, node) );
#endif

      XPETRA_FACTORY_END;
    }

    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, size_t numLocalElements, int indexBase, std::vector<size_t>& stridingInfo, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, LocalOrdinal stridedBlockId, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {

#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra)       
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, numLocalElements, indexBase, stridingInfo, comm, stridedBlockId, node) );
#endif

      XPETRA_FACTORY_END;
    }

#if 0 // TODO
    static RCP<StridedMap<LocalOrdinal,GlobalOrdinal, Node> > Build(UnderlyingLib lib, global_size_t numGlobalElements, const Teuchos::ArrayView<const int> &elementList, int indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Kokkos::DefaultNode::DefaultNodeType> &node = Kokkos::DefaultNode::getDefaultNode()) {
#ifdef HAVE_XPETRA_TPETRA
      if (lib == UseTpetra) 
        return rcp( new StridedTpetraMap<LocalOrdinal,GlobalOrdinal, Node> (numGlobalElements, elementList, indexBase, comm, node) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (lib == UseEpetra)
        return rcp( new StridedEpetraMap(numGlobalElements, elementList, indexBase, comm, node) );
#endif
      XPETRA_FACTORY_END;
    }
#endif

 
  };

}

#define XPETRA_STRIDEDMAPFACTORY_SHORT
#endif
//TODO: removed unused methods
