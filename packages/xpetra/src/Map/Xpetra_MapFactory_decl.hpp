// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MAPFACTORY_DECL_HPP
#define XPETRA_MAPFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Map_decl.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

/// \class MapFactory
/// \brief Create an Xpetra::Map instance.
///
/// Users must specify the exact class of the object that they want
/// to create (either an Xpetra::TpetraMap or an Xpetra::EpetraMap).
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = typename Map<LocalOrdinal, GlobalOrdinal>::node_type>
class MapFactory {
 private:
  //! Private constructor. This is a static class.
  MapFactory() {}

 public:
  //! Map constructor with Xpetra-defined contiguous uniform distribution.

  static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        GlobalOrdinal indexBase,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        LocalGlobal lg = Xpetra::GloballyDistributed);

  //! Map constructor with a user-defined contiguous distribution.

  static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        size_t numLocalElements,
        GlobalOrdinal indexBase,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Map constructor with user-defined non-contiguous (arbitrary) distribution.

  static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
        GlobalOrdinal indexBase,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  /*!
    @brief Map constructor generating degrees of freedom with numDofPerNode for given nodeMap

    @param[in] nodeMap Existing (node) map
    @param[in] numDofPerNode Number of DOFs per node for output map
    @param[in] gidOffset GID offset for output map
    @return Map

    \note This acts like a deep copy.
  */
  static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& nodeMap,
        const LocalOrdinal numDofPerNode,
        const GlobalOrdinal gidOffset = Teuchos::ScalarTraits<GlobalOrdinal>::zero());

#ifdef HAVE_XPETRA_TPETRA
  static Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node>>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        const Kokkos::View<const GlobalOrdinal*, typename Node::device_type>& indexList,
        GlobalOrdinal indexBase,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm);
#endif

  //! Create a locally replicated Map with the default node.
  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createLocalMap(UnderlyingLib lib,
                 size_t numElements,
                 const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a locally replicated Map with a specified node.

  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createLocalMapWithNode(UnderlyingLib lib,
                         size_t numElements,
                         const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a uniform, contiguous Map with a user-specified node.

  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createUniformContigMapWithNode(UnderlyingLib lib,
                                 global_size_t numElements,
                                 const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a uniform, contiguous Map with the default node.
  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createUniformContigMap(UnderlyingLib lib,
                         global_size_t numElements,
                         const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a (potentially) non-uniform, contiguous Map with the default node.
  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createContigMap(UnderlyingLib lib,
                  global_size_t numElements,
                  size_t localNumElements,
                  const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a (potentially) non-uniform, contiguous Map with a user-specified node.

  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  createContigMapWithNode(UnderlyingLib lib,
                          global_size_t numElements,
                          size_t localNumElements,
                          const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  //! Create a copy of the map, only using the new Comm object *if* the Comm would be valid
  // for this map.
  static Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
  copyMapWithNewComm(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& oldmap,
                     const Teuchos::RCP<const Teuchos::Comm<int>>& newComm);

};  // class MapFactory

//////////////////////////////////////////////////////////////
///  X P E T R A   E P E T R A   S P E C I A L I Z A T I O N
//////////////////////////////////////////////////////////////

}  // namespace Xpetra

#define XPETRA_MAPFACTORY_SHORT

#endif  // XPETRA_MAPFACTORY_DECL_HPP

// TODO: removed unused methods
