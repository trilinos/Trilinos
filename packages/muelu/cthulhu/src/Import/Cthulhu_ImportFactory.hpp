#ifndef CTHULHU_IMPORTFACTORY_HPP
#define CTHULHU_IMPORTFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Import.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraImport.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraImport.hpp"
#endif

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {
  
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ImportFactory {
    
  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}

  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() == target->lib(), Cthulhu::Exceptions::RuntimeError, "");

#ifdef HAVE_CTHULHU_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

      CTHULHU_FACTORY_ERROR_IF_EPETRA(source->lib());
      CTHULHU_FACTORY_END;
    }
    
  };

  template <>
  class ImportFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    
  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}
    
  public:
    
    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() == target->lib(), Cthulhu::Exceptions::RuntimeError, "");
    
#ifdef HAVE_CTHULHU_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif
    
#ifdef HAVE_CTHULHU_EPETRA
      if (source->lib() == UseEpetra)
        return rcp( new EpetraImport(source, target));
#endif

      CTHULHU_FACTORY_END;
    }
    
  };

}

#define CTHULHU_IMPORTFACTORY_SHORT
#endif
