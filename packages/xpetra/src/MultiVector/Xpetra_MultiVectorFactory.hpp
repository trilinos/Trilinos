#ifndef XPETRA_MULTIVECTOR_FACTORY_HPP
#define XPETRA_MULTIVECTOR_FACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_MultiVector.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraMultiVector.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MultiVectorFactory {

  private:
    //! Private constructor. This is a static class. 
    MultiVectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut=true) {

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
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

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, NumVectors, zeroOut) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraMultiVector(map, NumVectors, zeroOut) );
#endif

      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_MULTIVECTORFACTORY_SHORT
#endif
