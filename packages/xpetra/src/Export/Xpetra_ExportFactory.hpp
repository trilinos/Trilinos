#ifndef XPETRA_EXPORTFACTORY_HPP
#define XPETRA_EXPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Export.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraExport.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraExport.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ExportFactory {
    
  private:
    //! Private constructor. This is a static class. 
    ExportFactory() {}

  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
      XPETRA_FACTORY_END;
    }
    
  };

  template <>
  class ExportFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    
  private:
    //! Private constructor. This is a static class. 
    ExportFactory() {}
    
  public:
    
    static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");
    
#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif
    
#ifdef HAVE_XPETRA_EPETRA
      if (source->lib() == UseEpetra)
        return rcp( new EpetraExport(source, target));
#endif

      XPETRA_FACTORY_END;
    }
    
  };

}

#define XPETRA_EXPORTFACTORY_SHORT
#endif
