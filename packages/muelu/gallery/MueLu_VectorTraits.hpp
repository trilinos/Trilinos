/*
  Support for vectors.
*/

// TODO: rename variables (camelCase)

#ifndef MUELU_VECTORTRAITS_HPP
#define MUELU_VECTORTRAITS_HPP

#ifdef XPETRA_ENABLED
#  include "Xpetra_Map.hpp"  // needed for specialized traits
#endif

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {
  
  namespace Gallery {
    
    // Default traits
    template <class Map, class Vector>
    class VectorTraits 
    {
    public:
      static RCP<Vector> Build(const RCP<const Map> &map, size_t numVectors, bool zeroOut) {
        return rcp( new Vector(map, numVectors, zeroOut) );
      }
    };

#ifdef XPETRA_ENABLED

/*
    // Specialized traits for     Map = Xpetra::Map<...>, Vector = Xpetra::Vector<...>
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node >
    class VectorTraits <Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node > >
    {
    public:
      static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
        Build(const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &map, bool zeroOut)
      { return Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node>::Build(map,  zeroOut); };
    };
*/

    // Specialized traits for     Map = Xpetra::Map<...>, Vector = Xpetra::MultiVector<...>
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class VectorTraits <Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
    {
    public:
      static RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
        Build(const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &map, size_t numVectors, bool zeroOut)
      { return Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node>::Build(map, numVectors, zeroOut);};
    };

#endif

  } // namespace Gallery

} // namespace MueLu

#endif //ifndef MUELU_VECTORTRAITS_HPP
