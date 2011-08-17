#ifndef CTHULHU_MULTIVECTOR_FACTORY_HPP
#define CTHULHU_MULTIVECTOR_FACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_MultiVector.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMultiVector.hpp"
#endif

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMultiVector.hpp"
#endif

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {
  
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MultiVectorFactory {

  private:
    //! Private constructor. This is a static class. 
    MultiVectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut=true) {

#ifdef HAVE_CTHULHU_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      CTHULHU_FACTORY_ERROR_IF_EPETRA(map->lib());
      CTHULHU_FACTORY_END;
    }

  };

  template <>
  class MultiVectorFactory<double, int, int> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
  
  private:
    //! Private constructor. This is a static class. 
    MultiVectorFactory() {}
    
  public:
    
    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut=true) {

#ifdef HAVE_CTHULHU_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

#ifdef HAVE_CTHULHU_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVector(map, NumVectors, zeroOut) );
#endif

      CTHULHU_FACTORY_END;
    }

  };

}

#define CTHULHU_MULTIVECTORFACTORY_SHORT
#endif
