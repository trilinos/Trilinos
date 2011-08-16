#ifndef CTHULHU_CRSMATRIX_FACTORY_HPP
#define CTHULHU_CRSMATRIX_FACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_CrsMatrix.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraCrsMatrix.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraCrsMatrix.hpp"
#endif

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class CrsMatrixFactory {

  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) {

#ifdef HAVE_CTHULHU_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype) );
#endif

      CTHULHU_FACTORY_ERROR_IF_EPETRA(rowMap->lib());
      CTHULHU_FACTORY_END;
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
    
    static RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) {

#ifdef HAVE_CTHULHU_TPETRA
      if (rowMap->lib() == UseTpetra)
        return rcp( new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(rowMap, maxNumEntriesPerRow, pftype) );
#endif

#ifdef HAVE_CTHULHU_EPETRA
      if (rowMap->lib() == UseEpetra)
        return rcp( new EpetraCrsMatrix(rowMap, maxNumEntriesPerRow, pftype) );
#endif

      CTHULHU_FACTORY_END;
    }

  };

}

#define CTHULHU_CRSMATRIX_FACTORY_SHORT
#endif
