// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_MAPTRAITS_HPP
#define GALERI_MAPTRAITS_HPP

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
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMap.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif
#endif // HAVE_GALERI_XPETRA

namespace Galeri {

  namespace Xpetra {

    typedef size_t global_size_t;

    // TODO: Epetra_Map trait not implemented

    template <typename T>
    struct UndefinedMapTraits
    {
      static inline T notDefined() { return T::this_type_is_missing_a_specialization(); }
    };

    /* Default traits (not implemented) */
    template <class GlobalOrdinal, class Map>
    class MapTraits
    {
    public:
      static Teuchos::RCP<Map> Build(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return UndefinedMapTraits<Map>::notDefined(); }

      static Teuchos::RCP<Map> Build(global_size_t numGlobalElements, global_size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return UndefinedMapTraits<Map>::notDefined(); }
    };

#ifdef HAVE_GALERI_TPETRA
    /* Specialized traits for Map = Tpetra::Map<...> */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    class MapTraits < GlobalOrdinal, Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    {
    public:
      static Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Build(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, elementList, indexBase, comm) ); }

      static Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > Build(global_size_t numGlobalElements, global_size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return rcp( new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElements, numLocalElements, indexBase, comm) ); }
    };
#endif // HAVE_GALERI_TPETRA

#ifdef HAVE_GALERI_XPETRA
#ifdef HAVE_XPETRA_TPETRA
    /* Specialized traits for Map = Xpetra::TpetraMap<...> */
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    class MapTraits <GlobalOrdinal, ::Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node> >
    {
    public:
      static Teuchos::RCP< ::Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node> > Build(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return Teuchos::rcp( new ::Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node>(numGlobalElements, elementList, indexBase, comm) ); }

      static Teuchos::RCP< ::Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node> > Build(global_size_t numGlobalElements, global_size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return Teuchos::rcp( new ::Xpetra::TpetraMap<LocalOrdinal,GlobalOrdinal, Node>(numGlobalElements, numLocalElements, indexBase, comm) ); }
    };
#endif

#ifdef HAVE_XPETRA_EPETRA
    /* Specialized traits for Map = Xpetra::EpetraMap<int,GlobalOrdinal,Node> */
    template <class GlobalOrdinal, class Node>
    class MapTraits <GlobalOrdinal, ::Xpetra::EpetraMapT<GlobalOrdinal,Node> >
    {
    public:
      static Teuchos::RCP< ::Xpetra::EpetraMapT<GlobalOrdinal,Node> > Build(global_size_t numGlobalElements, const Teuchos::ArrayView<const GlobalOrdinal> &elementList, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return Teuchos::rcp( new ::Xpetra::EpetraMapT<GlobalOrdinal,Node>(numGlobalElements, elementList, indexBase, comm) ); }

      static Teuchos::RCP< ::Xpetra::EpetraMapT<GlobalOrdinal,Node> > Build(global_size_t numGlobalElements, global_size_t numLocalElements, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
      { return Teuchos::rcp( new ::Xpetra::EpetraMapT<GlobalOrdinal,Node>(numGlobalElements, numLocalElements, indexBase, comm) ); }
    };
#endif

#endif // HAVE_GALERI_XPETRA

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_MAPTRAITS_HPP
