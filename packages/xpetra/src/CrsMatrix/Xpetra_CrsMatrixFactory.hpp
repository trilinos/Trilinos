#ifndef XPETRA_CRSMATRIX_FACTORY_HPP
#define XPETRA_CRSMATRIX_FACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_CrsMatrix.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraCrsMatrix.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraCrsMatrix.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class CrsMatrixFactory {

  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile) {

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      XPETRA_FACTORY_END;
    }
    
  };

  template <>
  class CrsMatrixFactory<double, int, int> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile) {

#ifdef HAVE_XPETRA_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, maxNumEntriesPerRow, pftype) );
#endif

      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_CRSMATRIX_FACTORY_SHORT
#endif
