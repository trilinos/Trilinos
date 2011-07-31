#ifndef CTHULHU_IMPORTFACTORY_HPP
#define CTHULHU_IMPORTFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_Import.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraImport.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraImport.hpp"
#endif

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {
  
  template <class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>

  class ImportFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
    typedef Import<LocalOrdinal, GlobalOrdinal, Node> ImportClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraImport<LocalOrdinal, GlobalOrdinal, Node> TpetraImportClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<ImportClass> Build(const Teuchos::RCP<const MapClass> &source, const Teuchos::RCP<const MapClass> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMapClass> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(source);
      const RCP<const TpetraMapClass> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(target);
      if (tSource != null && tTarget != null)
        return rcp( new TpetraImportClass(tSource, tTarget) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap.");
    }
    
  };

  template <>
  class ImportFactory<int, int, Kokkos::DefaultNode::DefaultNodeType> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> MapClass;
    typedef Import<int, int, Kokkos::DefaultNode::DefaultNodeType> ImportClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMapClass;
    typedef TpetraImport<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraImportClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<ImportClass> Build(const Teuchos::RCP<const MapClass> &source, const Teuchos::RCP<const MapClass> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      {
        const RCP<const TpetraMapClass> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(source);
        const RCP<const TpetraMapClass> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(target);
        if (tSource != null && tTarget != null)
          return rcp( new TpetraImportClass(tSource, tTarget) );
      }
#endif
#ifdef HAVE_CTHULHU_EPETRA
      {
        const RCP<const EpetraMap> &eSource = Teuchos::rcp_dynamic_cast<const EpetraMap>(source);
        const RCP<const EpetraMap> &eTarget = Teuchos::rcp_dynamic_cast<const EpetraMap>(target);
        if (eSource != null && eTarget != null)
          return rcp( new EpetraImport(eSource, eTarget) );
      }
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap.");
    }

  };


}

#define CTHULHU_IMPORTFACTORY_SHORT
#endif
