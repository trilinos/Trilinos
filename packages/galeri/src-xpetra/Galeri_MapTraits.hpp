// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
#include <Xpetra_TpetraMap.hpp>
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
