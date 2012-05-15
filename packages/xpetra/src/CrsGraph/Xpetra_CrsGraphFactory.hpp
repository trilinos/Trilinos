#ifndef XPETRA_CRSGRAPHFACTORY_HPP
#define XPETRA_CRSGRAPHFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsGraph.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsGraph.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsGraph.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class CrsGraphFactory {

  private:
    //! Private constructor. This is a static class. 
    CrsGraphFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, ProfileType pftype=DynamicProfile) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, pftype) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
      return null;
    }

  };

  template <>
  class CrsGraphFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
  
  private:
    //! Private constructor. This is a static class. 
    CrsGraphFactory() {}
    
  public:
    
    static RCP<CrsGraph<LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, ProfileType pftype=DynamicProfile) {
      XPETRA_MONITOR("CrsGraphFactory::Build");

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, pftype) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraCrsGraph(map, NumVectors, pftype) );
#endif

      XPETRA_FACTORY_END;
      return null;
    }

  };

}

#define XPETRA_CRSGRAPHFACTORY_SHORT
#endif
