#ifndef CTHULHU_MULTIVECTOR_FACTORY_HPP
#define CTHULHU_MULTIVECTOR_FACTORY_HPP

#include "Cthulhu_ConfigDefs.hpp"

#include "Cthulhu_Map.hpp"
#include "Cthulhu_MultiVector.hpp"

#ifdef HAVE_CTHULHU_TPETRA
#include "Cthulhu_TpetraMap.hpp"
#include "Cthulhu_TpetraMultiVector.hpp"
#endif // HAVE_CTHULHU_TPETRA

#ifdef HAVE_CTHULHU_EPETRA
#include "Cthulhu_EpetraMultiVector.hpp"
#include "Cthulhu_EpetraMap.hpp"
#endif // HAVE_CTHULHU_EPETRA

#include "Cthulhu_Debug.hpp"

// This factory creates Cthulhu::MultiVector. User don't have to specify the exact class of object that he want to create (ie: a Cthulhu::TpetraMultiVector or a Cthulhu::EpetraMultiVector).
// Each Build() method takes at least one Cthulhu object in argument (a Map or another MultiVector) and Build() methods return a MultiVector created by using the same underlying library (Epetra or Tpetra).

namespace Cthulhu {
  
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType>

  class MultiVectorFactory {
    
    typedef Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
    typedef TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMultiVector;
#endif

  private:
    //! Private constructor. This is a static class. 
    MultiVectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<MultiVector> Build(const Teuchos::RCP<const Map> &map, size_t NumVectors, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(map);
      if (tMap != null)
        return rcp( new TpetraMultiVector(map, NumVectors, zeroOut) );
#endif

      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map. The exact type of the Map 'map' is unknown.");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
    
  };

  template <>
  class MultiVectorFactory<int, int, Kokkos::DefaultNode::DefaultNodeType> {
  
    typedef Map<int, int> Map;
    typedef MultiVector<double, int, int> MultiVector;
#ifdef HAVE_CTHULHU_TPETRA
    typedef TpetraMap<int, int> TpetraMap;
    typedef TpetraMultiVector<double, int, int> TpetraMultiVector;
#endif

  private:
    //! Private constructor. This is a static class. 
    MultiVectorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<MultiVector> Build(const Teuchos::RCP<const Map> &map, size_t NumVectors, bool zeroOut=true) {
#ifdef HAVE_CTHULHU_TPETRA
      const RCP<const TpetraMap> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMap>(map);
      if (tMap != null)
        return rcp( new TpetraMultiVector(map, NumVectors, zeroOut) );
#endif
#ifdef HAVE_CTHULHU_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(map);
      if (eMap != null)
        return rcp( new EpetraMultiVector(map, NumVectors, zeroOut) );
#endif

      TEST_FOR_EXCEPTION(1,Cthulhu::Exceptions::BadCast,"Cannot dynamically cast Cthulhu::Map. The exact type of the Map 'map' is unknown.");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    // Other constructors here
#endif
   
  };

}

#define CTHULHU_MULTIVECTORFACTORY_SHORT
#endif
