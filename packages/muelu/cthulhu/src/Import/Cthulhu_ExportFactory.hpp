#ifndef CTHULHU_EXPORTFACTORY_HPP
#define CTHULHU_EXPORTFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_Export.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraExport.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraExport.hpp"
#endif

namespace Cthulhu {
  
  template <class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>

  class ExportFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
    typedef Export<LocalOrdinal, GlobalOrdinal, Node> ExportClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraExport<LocalOrdinal, GlobalOrdinal, Node> TpetraExportClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    ExportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<ExportClass> Build(const Teuchos::RCP<const MapClass> &source, const Teuchos::RCP<const MapClass> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMapClass> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(source);
      const RCP<const TpetraMapClass> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(target);
      if (tSource != null && tTarget != null)
        return rcp( new TpetraExportClass(tSource, tTarget) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap.");
    }
    
  };

  template <>
  class ExportFactory<int, int, Kokkos::DefaultNode::DefaultNodeType> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> MapClass;
    typedef Export<int, int, Kokkos::DefaultNode::DefaultNodeType> ExportClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMapClass;
    typedef TpetraExport<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraExportClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    ExportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<ExportClass> Build(const Teuchos::RCP<const MapClass> &source, const Teuchos::RCP<const MapClass> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      {
        const RCP<const TpetraMapClass> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(source);
        const RCP<const TpetraMapClass> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(target);
        if (tSource != null && tTarget != null)
          return rcp( new TpetraExportClass(tSource, tTarget) );
      }
#endif
#ifdef HAVE_CTHULHU_EPETRA
      {
        const RCP<const EpetraMap> &eSource = Teuchos::rcp_dynamic_cast<const EpetraMap>(source);
        const RCP<const EpetraMap> &eTarget = Teuchos::rcp_dynamic_cast<const EpetraMap>(target);
        if (eSource != null && eTarget != null)
          return rcp( new EpetraExport(eSource, eTarget) );
      }
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap.");
    }

  };


}

#define CTHULHU_EXPORTFACTORY_SHORT
#endif
