#ifndef XPETRA_IMPORTFACTORY_HPP
#define XPETRA_IMPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Import.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraImport.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraImport.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ImportFactory {
    
  private:
    //! Private constructor. This is a static class. 
    ImportFactory() {}

  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
      XPETRA_FACTORY_END;
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
      TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");
    
#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif
    
#ifdef HAVE_XPETRA_EPETRA
      if (source->lib() == UseEpetra)
        return rcp( new EpetraImport(source, target));
#endif

      XPETRA_FACTORY_END;
    }
    
  };

}

#define XPETRA_IMPORTFACTORY_SHORT
#endif
