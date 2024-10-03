// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMAP_DEF_HPP
#define TPETRA_BLOCKEDMAP_DEF_HPP

#include <Kokkos_Core.hpp>

namespace Tpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const Teuchos::RCP<const map_type>& pointMap,
               const Teuchos::RCP<lo_vec_type>& blockSizes)
  : pointMap_(pointMap)
  , blockMap_(blockSizes->getMap())
  , blockSizes_(blockSizes) {
  auto lclBlockSizes          = blockSizes_->getLocalViewHost(Tpetra::Access::ReadOnly);
  LocalOrdinal minClusterSize = Teuchos::OrdinalTraits<LocalOrdinal>::max();
  LocalOrdinal maxClusterSize = 0;
  offsets_                    = Kokkos::View<size_t*>("offsets", blockMap_->getLocalNumElements() + 1);
  auto offsets_h              = Kokkos::create_mirror_view(offsets_);
  offsets_h(0)                = 0;
  for (size_t blockNum = 0; blockNum < blockMap_->getLocalNumElements(); ++blockNum) {
    offsets_h(blockNum + 1) = offsets_h(blockNum) + lclBlockSizes(blockNum, 0);
    minClusterSize          = std::min(minClusterSize, lclBlockSizes(blockNum, 0));
    maxClusterSize          = std::max(maxClusterSize, lclBlockSizes(blockNum, 0));
  }
  Kokkos::deep_copy(offsets_, offsets_h);
  TEUCHOS_ASSERT_EQUALITY(offsets_h(blockMap_->getLocalNumElements()), pointMap->getLocalNumElements());
  minClusterSize_ = minClusterSize;
  maxClusterSize_ = maxClusterSize;
}

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMAP_DEF_HPP
