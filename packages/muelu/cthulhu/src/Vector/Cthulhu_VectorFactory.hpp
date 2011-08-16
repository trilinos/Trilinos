#ifndef CTHULHU_VECTORFACTORY_HPP
#define CTHULHU_VECTORFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Vector.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraVector.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraVector.hpp"
#include "Cthulhu_EpetraIntVector.hpp"
#endif

#include "Cthulhu_Exceptions.hpp"

namespace Cthulhu {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class VectorFactory {
    
  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, bool zeroOut=true) {

#ifdef HAVE_CTHULHU_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

      CTHULHU_FACTORY_ERROR_IF_EPETRA(map->lib());
      CTHULHU_FACTORY_END;
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

#ifdef HAVE_CTHULHU_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

#ifdef HAVE_CTHULHU_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraVector(map, zeroOut) );
#endif

      CTHULHU_FACTORY_END;
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

#ifdef HAVE_CTHULHU_TPETRA
      if (map->lib() == UseTpetra)
        return rcp( new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> (map, zeroOut) );
#endif

#ifdef HAVE_CTHULHU_EPETRA
      if (map->lib() == UseEpetra)
        return rcp( new EpetraIntVector(map, zeroOut) );
#endif

      CTHULHU_FACTORY_END;
    }

  };

}

#define CTHULHU_VECTORFACTORY_SHORT
#endif
// TODO: one factory for both Vector and MultiVector ?
