// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_STRIDEDMAPFACTORY_DECL_HPP
#define XPETRA_STRIDEDMAPFACTORY_DECL_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_StridedMap_decl.hpp"

namespace Xpetra {

/*!
  @class StridedMapFactory
  @brief This factory creates a \c Xpetra::StridedMap .

  @tparam LocalOrdinal
  @tparam GlobalOrdinal
  @tparam Node
*/
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class StridedMapFactory {
#undef XPETRA_STRIDEDMAPFACTORY_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"

 private:
  //! Private constructor. This is a static class.
  StridedMapFactory() = delete;

 public:
  //! Map constructor with Xpetra-defined contiguous uniform distribution.
  static RCP<Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node>>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        GlobalOrdinal indexBase,
        std::vector<size_t>& stridingInfo,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        LocalOrdinal stridedBlockId = -1,
        GlobalOrdinal offset        = 0,
        LocalGlobal lg              = Xpetra::GloballyDistributed);

  //! Map constructor with a user-defined contiguous distribution.
  static RCP<StridedMap>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        size_t numLocalElements,
        GlobalOrdinal indexBase,
        std::vector<size_t>& stridingInfo,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        LocalOrdinal stridedBlockId = -1,
        GlobalOrdinal offset        = 0);

  //! Build a strided map
  static RCP<StridedMap>
  Build(const RCP<const Map>& map, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0);

  //! special constructor for generating a given subblock of a strided map
  static RCP<StridedMap>
  Build(const RCP<const StridedMap>& map, LocalOrdinal stridedBlockId);

  //! Create copy of existing map (this just creates a copy of your map, it's not a clone in the sense of Tpetra)
  static RCP<StridedMap>
  Build(const StridedMap& map);

  /*!
    @brief Map constructor with a user-defined contiguous distribution.

    \warning For experts only. There is no special check whether the generated strided maps are valid.
  */
  static RCP<StridedMap>
  Build(UnderlyingLib lib,
        global_size_t numGlobalElements,
        const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
        GlobalOrdinal indexBase,
        std::vector<size_t>& stridingInfo,
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        LocalOrdinal stridedBlockId = -1,  // FIXME (mfh 03 Sep 2014) This breaks if LocalOrdinal is unsigned
        GlobalOrdinal /* offset */  = 0);

};  // class StridedMapFactory

}  // namespace Xpetra

#define XPETRA_STRIDEDMAPFACTORY_SHORT
#endif  // XPETRA_STRIDEDMAPFACTORY_DECL_HPP

// TODO: removed unused methods
