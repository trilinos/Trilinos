#ifndef CTHULHU_CRSMATRIX_FACTORY_HPP
#define CTHULHU_CRSMATRIX_FACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_CrsMatrix.hpp"
#include "Cthulhu_Map.hpp"
#include "Cthulhu_Exceptions.hpp"

//
#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraCrsMatrix.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraCrsMatrix.hpp"
#endif

#include "Cthulhu_Debug.hpp"

// This factory creates Cthulhu::CrsMatrix. User don't have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraCrsMatrix or a Cthulhu::EpetraCrsMatrix).
// Each Build() method takes at least one Cthulhu object in argument (a Map or another CrsMatrix) and Build() methods return a CrsMatrix created by using the same underlying library (Epetra or Tpetra).

namespace Cthulhu {
  
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>

  class CrsMatrixFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Cthulhu::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
    typedef TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif

  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<CrsMatrix> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tRowMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(rowMap);
      if (tRowMap != null)
        return rcp( new TpetraCrsMatrix(rowMap, maxNumEntriesPerRow) ); //TODO: convert pftype
#endif
#ifdef HAVE_CTHULHU_EPETRA
      // TODO: explicit error message
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map. The exact type of the Map 'rowMap' is unknown.");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying the number of non-zeros for each row.
    Build(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { }

    //! Constructor specifying a column map and the number of non-zeros for all rows.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { }
    
    //! Constructor specifying a column map and the number of non-zeros for each row.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { } 

    //! Constructor specifying a pre-constructed graph.
    // TODO: need a CrsGraph
    // Build() {}

    //! Constructor specifying a Epetra CrsMatrix
    // ?? Build(const Teuchos::RCP<const Cthulhu::CrsMatrix> &mtx) { }

    //! Constructor specifying a Epetra CrsMatrix
    Build(const Teuchos::RCP<const Epetra_CrsMatrix> &mtx) { }
#endif
    
  };

  template <>
  class CrsMatrixFactory<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> Map;
    typedef Cthulhu::CrsMatrix<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> CrsMatrix;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMap;
    typedef TpetraCrsMatrix<double, int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> TpetraCrsMatrix;
#endif

  private:
    //! Private constructor. This is a static class. 
    CrsMatrixFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<CrsMatrix> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tRowMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(rowMap);
      if (tRowMap != null)
        return rcp( new TpetraCrsMatrix(rowMap, maxNumEntriesPerRow) ); //TODO: convert pftype
#endif
#ifdef HAVE_CTHULHU_EPETRA
      const RCP<const EpetraMap> &eRowMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(rowMap);
      if (eRowMap != null)
        return rcp( new EpetraCrsMatrix(rowMap, maxNumEntriesPerRow) ); //TODO: convert pftype
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map. The exact type of the Map 'rowMap' is unknown.");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Constructor specifying the number of non-zeros for each row.
    Build(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { }

    //! Constructor specifying a column map and the number of non-zeros for all rows.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, size_t maxNumEntriesPerRow, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { }
    
    //! Constructor specifying a column map and the number of non-zeros for each row.
    /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
     */
    Build(const RCP<const Map> &rowMap, const RCP<const Map> &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, Cthulhu::ProfileType pftype = Cthulhu::DynamicProfile) { } 

    Build(const Teuchos::RCP<const Epetra_CrsMatrix> &mtx) { }
#endif
    
  };

}

#define CTHULHU_CRSMATRIX_FACTORY_SHORT
#endif
