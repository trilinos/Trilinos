// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_MULTIVECTORTRAITS_HPP
#define GALERI_MULTIVECTORTRAITS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Galeri_ConfigDefs.h"

#ifdef HAVE_GALERI_TPETRA
#include <Tpetra_Map.hpp>
#endif

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_MapFactory.hpp>
#endif

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapExtractor.hpp>
#endif // HAVE_GALERI_XPETRA

namespace Galeri {

  namespace Xpetra {

    typedef size_t global_size_t;

    // TODO: Epetra_Map trait not implemented

    template <typename T>
    struct UndefinedMultiVectorTraits
    {
      static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
    };

    /* Default traits (not implemented) */
    template <class Map, class MultiVector>
    class MultiVectorTraits {
    public:
      typedef void type;

      static Teuchos::RCP<MultiVector> Build(const Teuchos::RCP<const Map>& map, size_t num) { return UndefinedMultiVectorTraits<MultiVector>::notDefined(); }
    };

#ifdef HAVE_GALERI_TPETRA
    /* Specialized traits for Map = Tpetra::Map<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    class MultiVectorTraits<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>,Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
    public:
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> type;
      typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LocalOrdinal,GlobalOrdinal,Node> type_real;

      static Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Build(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& map, size_t num) {
        return Teuchos::rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map, num));
      }
    };
#endif // HAVE_GALERI_TPETRA

#ifdef HAVE_GALERI_XPETRA
    /* Specialized traits for Map = Xpetra::TpetraMap<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class Map>
    class MultiVectorTraits<Map,::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > {
    public:
      typedef ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> type;
      typedef ::Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LocalOrdinal,GlobalOrdinal,Node> type_real;

      static Teuchos::RCP< ::Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Build(const Teuchos::RCP<const Map>& map, size_t num) {
        return ::Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(map, num);
      }
    };
#endif // HAVE_GALERI_XPETRA

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_MULTIVECTORTRAITS_HPP
