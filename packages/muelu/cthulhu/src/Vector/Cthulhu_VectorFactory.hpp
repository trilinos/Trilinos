#ifndef CTHULHU_VECTORFACTORY_HPP
#define CTHULHU_VECTORFACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_Vector.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraVector.hpp"
#endif
#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraVector.hpp"
#include "Cthulhu_EpetraIntVector.hpp"
#endif

// This factory creates Cthulhu::Vector. User don't have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraVector or a Cthulhu::EpetraVector).
// Each Build() method takes at least one Cthulhu object in argument (a Map or another Vector) and Build() methods return a Vector created by using the same underlying library (Epetra or Tpetra).

namespace Cthulhu {
  
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>

  class VectorFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
    typedef Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> VectorClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVectorClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<VectorClass> Build(const Teuchos::RCP<const MapClass> &map, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(map);
      if (tMap != null)
        return rcp( new TpetraVectorClass(map, zeroOut) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap. The exact type of the Map 'map' is unknown");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
    
  };

  template <>
  class VectorFactory<double, int, int, Kokkos::DefaultNode::DefaultNodeType> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> MapClass;
    typedef Vector<double, int, int, Kokkos::DefaultNode::DefaultNodeType> VectorClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMapClass;
    typedef TpetraVector<double, int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraVectorClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<VectorClass> Build(const Teuchos::RCP<const MapClass> &map, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(map);
      if (tMap != null)
        return rcp( new TpetraVectorClass(map, zeroOut) );
#endif
#ifdef HAVE_CTHULHU_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(map);
      if (eMap != null)
        return rcp( new EpetraVector(map, zeroOut) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap. The exact type of the Map 'map' is unknown");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
    
  };

  template <>
  class VectorFactory<int, int, int, Kokkos::DefaultNode::DefaultNodeType> {
    
    typedef Map<int, int, Kokkos::DefaultNode::DefaultNodeType> MapClass;
    typedef Vector<int, int, int, Kokkos::DefaultNode::DefaultNodeType> VectorClass;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraMapClass;
    typedef TpetraVector<int, int, int, Kokkos::DefaultNode::DefaultNodeType> TpetraVectorClass;
#endif

  private:
    //! Private constructor. This is a static class. 
    VectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<VectorClass> Build(const Teuchos::RCP<const MapClass> &map, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(map);
      if (tMap != null)
        return rcp( new TpetraVectorClass(map, zeroOut) );
#endif
#ifdef HAVE_CTHULHU_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(map);
      if (eMap != null)
        return rcp( new EpetraIntVector(map, zeroOut) );
#endif
      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map to an EpetraMap or a TpetraMap. The exact type of the Map 'map' is unknown");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
    
  };

}

// TODO: Only one factory for Vector and MultiVector ?? -> Yes

#define CTHULHU_VECTORFACTORY_SHORT
#endif
