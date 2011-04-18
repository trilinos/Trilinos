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
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Import<LocalOrdinal, GlobalOrdinal, Node> Import;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
    typedef TpetraImport<LocalOrdinal, GlobalOrdinal, Node> TpetraImport;
#endif

  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Import> Build(const Teuchos::RCP<const Map> &source, const Teuchos::RCP<const Map> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMap>(source);
      const RCP<const TpetraMap> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMap>(target);
      if (tSource != null && tTarget != null)
        return rcp( new TpetraImport(tSource, tTarget) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap.");
    }
    
  };

  template <>
  class ImportFactory<int, int, Kokkos::DefaultNode::DefaultNodeType> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> Map;
    typedef Import<int, int, Kokkos::DefaultNode::DefaultNodeType> Import;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMap;
    typedef TpetraImport<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraImport;
#endif

  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Import> Build(const Teuchos::RCP<const Map> &source, const Teuchos::RCP<const Map> &target) {
#ifdef HAVE_CTHULHU_TPETRA
      {
        const RCP<const TpetraMap> &tSource = Teuchos::rcp_dynamic_cast<const TpetraMap>(source);
        const RCP<const TpetraMap> &tTarget = Teuchos::rcp_dynamic_cast<const TpetraMap>(target);
        if (tSource != null && tTarget != null)
          return rcp( new TpetraImport(tSource, tTarget) );
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
