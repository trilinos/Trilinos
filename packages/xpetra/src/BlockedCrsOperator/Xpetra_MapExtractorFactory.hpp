/*
 * Xpetra_MapExtractorFactory.hpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

#ifndef XPETRA_MAPEXTRACTORFACTORY_HPP_
#define XPETRA_MAPEXTRACTORFACTORY_HPP_

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMapExtractor.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMapExtractor.hpp>
#endif

namespace Xpetra
{
  // factory class
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MapExtractorFactory {
    typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
#ifdef HAVE_XPETRA_TPETRA
    typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMapClass;
    typedef TpetraMapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraMapExtractorClass;
#endif
  private:
    //! Private construtor, since this is a static class
    MapExtractorFactory() {}

  public:
    //! Constructor specifying the maps (and indirectly the used LINALG library (Tpetra vs Epetra))
    static Teuchos::RCP<MapExtractorClass> Build(const Teuchos::RCP<const MapClass>& fullmap, const std::vector<Teuchos::RCP<const MapClass> >& maps)
    {
#ifdef HAVE_XPETRA_EPETRA
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(fullmap);
      if (eMap != null)
        return Teuchos::rcp(new Xpetra::EpetraMapExtractor(eMap,maps));
#endif
#ifdef HAVE_XPETRA_TPETRA
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
      {
        Teuchos::RCP<TpetraMapExtractorClass> tMapExtractor = Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));
        return Teuchos::rcp_dynamic_cast<MapExtractorClass>(tMapExtractor);
      }
#endif

      TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"Cannot dynamically cast Xpetra::Map. The exact type of the Map 'map' is unknown.");
      return Teuchos::null;
    }

    //! Constructor for Tpetra
#ifdef HAVE_XPETRA_TPETRA
    static Teuchos::RCP<TpetraMapExtractorClass> Build(const Teuchos::RCP<const TpetraMapClass>& fullmap, const std::vector<Teuchos::RCP<const TpetraMapClass> >& maps)
    {
      const Teuchos::RCP<const TpetraMapClass> &tMap = Teuchos::rcp_dynamic_cast<const TpetraMapClass>(fullmap);
      if (tMap != null)
        //return Teuchos::rcp( new TpetraMultiVectorClass(map, NumVectors, zeroOut) );
        return Teuchos::rcp(new TpetraMapExtractorClass(tMap,maps));

      TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"Cannot dynamically cast Xpetra::Map. The exact type of the Map 'map' is unknown.");
      return Teuchos::null;
    }
#endif

    //! Constructor for Epetra
#ifdef HAVE_XPETRA_EPETRA
    static Teuchos::RCP<Xpetra::EpetraMapExtractor> Build(const Teuchos::RCP<const Xpetra::EpetraMap>& fullmap, const std::vector<Teuchos::RCP<const Xpetra::EpetraMap> >& maps)
    {
      const RCP<const EpetraMap> &eMap = Teuchos::rcp_dynamic_cast<const EpetraMap>(fullmap);
      if (eMap != null)
        return Teuchos::rcp(new Xpetra::EpetraMapExtractor(eMap,maps));

      TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"Cannot dynamically cast Xpetra::Map. The exact type of the Map 'map' is unknown.");
      return Teuchos::null;
    }
#endif
  };
}


#endif /* XPETRA_MAPEXTRACTORFACTORY_HPP_ */
