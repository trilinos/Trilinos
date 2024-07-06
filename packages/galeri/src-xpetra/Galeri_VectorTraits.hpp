// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Support for vectors.
*/

// TODO: rename variables (camelCase)
#include "Galeri_config.h"

#ifndef GALERI_XPETRAVECTORTRAITS_HPP
#define GALERI_XPETRAVECTORTRAITS_HPP

#ifdef HAVE_GALERI_XPETRA
#  include "Xpetra_Map.hpp"  // needed for specialized traits
#  include "Xpetra_BlockedMultiVector.hpp"
#  include "Xpetra_BlockedVector.hpp"
#  include "Xpetra_MultiVectorFactory.hpp"
#endif

namespace Galeri {

  namespace Xpetra {

    using Teuchos::RCP;

    // Default traits
    template <class Map, class Vector>
    class VectorTraits
    {
    public:
      static RCP<Vector> Build(const RCP<const Map> &map, size_t numVectors, bool zeroOut) {
        return rcp( new Vector(map, numVectors, zeroOut) );
      }
    };

#ifdef HAVE_GALERI_XPETRA

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
    class VectorTraits < ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
    {
    public:
      static RCP< ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal, Node> >
      Build(const RCP<const ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &map, size_t numVectors, bool zeroOut)
      { return ::Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node>::Build(map, numVectors, zeroOut);};
    };

#endif

  } // namespace Xpetra

} // namespace Galeri

#endif //ifndef GALERI_XPETRAVECTORTRAITS_HPP
