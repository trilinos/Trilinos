#ifndef CTHULHU_EXPORTFACTORY_HPP
#define CTHULHU_EXPORTFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Export.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraExport.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraExport.hpp"
#endif

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ExportFactory {
    
  private:
    //! Private constructor. This is a static class. 
    ExportFactory() {}

  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      TEST_FOR_EXCEPTION(source->lib() == target->lib(), Cthulhu::Exceptions::RuntimeError, "");

#ifdef HAVE_CTHULHU_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

      CTHULHU_FACTORY_ERROR_IF_EPETRA(source->lib());
      CTHULHU_FACTORY_END;
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
      TEST_FOR_EXCEPTION(source->lib() == target->lib(), Cthulhu::Exceptions::RuntimeError, "");
    
#ifdef HAVE_CTHULHU_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif
    
#ifdef HAVE_CTHULHU_EPETRA
      if (source->lib() == UseEpetra)
        return rcp( new EpetraExport(source, target));
#endif

      CTHULHU_FACTORY_END;
    }
    
  };

}

#define CTHULHU_EXPORTFACTORY_SHORT
#endif
