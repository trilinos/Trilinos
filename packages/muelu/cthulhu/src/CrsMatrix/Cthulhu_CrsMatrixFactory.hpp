#ifndef CTHULHU_CRSMATRIX_FACTORY_DECL_HPP
#define CTHULHU_CRSMATRIX_FACTORY_DECL_HPP

#include "Cthulhu_Classes.hpp"

#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_TpetraCrsMatrix.hpp"
//#include "Cthulhu_EpetraCrsMatrix.hpp"

//#include "Cthulhu_Map.hpp"
#include "Cthulhu_TpetraMap.hpp"
//#include "Cthulhu_EpetraMap.hpp"

#include "Cthulhu_Debug.hpp"

// This factory creates Cthulhu::CrsMatrix. User don't have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraCrsMatrix or a Cthulhu::EpetraCrsMatrix).
// Each Build() method takes at least one Cthulhu object in argument (a Map or another CrsMatrix) and Build() methods return a CrsMatrix created by using the same underlying library (Epetra or Tpetra).

namespace Cthulhu {
  
  template <class ScalarType, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<ScalarType,LocalOrdinal,Node>::SparseOps>

  class CrsMatrixFactory {
    
#include "Cthulhu_UseShortNames.hpp"
    
  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<CrsMatrix> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Tpetra::ProfileType pftype = Tpetra::DynamicProfile) {

      const RCP<const TpetraMap> &tRowMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(rowMap);
      if (tRowMap != null)
        return rcp( new TpetraCrsMatrix(rowMap, maxNumEntriesPerRow, pftype) );

      //TODO
      //       const RCP<const EpetraMap> &eRowMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(rowMap);
      //       if (eRowMap != null)
      //         return rcp( new EpetraCrsMatrix(rowMap, maxNumEntriesPerRow, pftype) );
      
      return null;
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying the number of non-zeros for each row.
    Build(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Tpetra::ProfileType pftype = Tpetra::DynamicProfile) { }

    //! Constructor specifying a column map and the number of non-zeros for all rows.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, size_t maxNumEntriesPerRow, Tpetra::ProfileType pftype = Tpetra::DynamicProfile) { }
    
    //! Constructor specifying a column map and the number of non-zeros for each row.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Tpetra::ProfileType pftype = Tpetra::DynamicProfile) { } 

    //! Constructor specifying a pre-constructed graph.
    // TODO: need a CrsGraph
    // Build() {}

    //! Constructor specifying a Tpetra CrsMatrix
    Build(const Teuchos::RCP<const Tpetra::CrsMatrix> &mtx) { }

    //! Constructor specifying a Epetra CrsMatrix
    Build(const Teuchos::RCP<const Epetra_CrsMatrix> &mtx) { }
#endif
    
  };

}

#endif
