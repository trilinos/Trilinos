#ifndef XPETRA_VECTORFACTORY_HPP
#define XPETRA_VECTORFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Vector.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraVector.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraVector.hpp"
#include "Xpetra_EpetraIntVector.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class VectorFactory {
    
  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut=true) {

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
      XPETRA_FACTORY_END;
    }
    
  };

  template <>
  class VectorFactory<double, int, int> {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    static RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut=true) {

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraVector(map, zeroOut) );
#endif

      XPETRA_FACTORY_END;
    }
    
  };

  template <>
  class VectorFactory<int, int, int> {
    
    typedef int Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    static RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut=true) {

#ifdef HAVE_XPETRA_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntVector(map, zeroOut) );
#endif

      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_VECTORFACTORY_SHORT
#endif
// TODO: one factory for both Vector and MultiVector ?
