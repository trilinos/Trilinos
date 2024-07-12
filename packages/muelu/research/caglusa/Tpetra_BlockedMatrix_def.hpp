// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMATRIX_DEF_HPP
#define TPETRA_BLOCKEDMATRIX_DEF_HPP

#include "Kokkos_DualView.hpp"
#include "Teuchos_Assert.hpp"
#include "Tpetra_Access.hpp"
namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMatrix(const Teuchos::RCP<matrix_type>& pointA,
                  const Teuchos::RCP<matrix_type>& blockA,
                  const Teuchos::RCP<blocked_map_type>& blockMap,
                  const Teuchos::RCP<blocked_map_type>& ghosted_blockMap)
  : pointA_(pointA)
  , blockA_(blockA)
  , blockMap_(blockMap)
  , ghosted_blockMap_(ghosted_blockMap) {
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRangeMap()));
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockA_->getRowMap()));
  TEUCHOS_ASSERT(blockA_->getDomainMap()->isSameAs(*blockMap_->blockMap_));

  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRangeMap()));
  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*pointA_->getRowMap()));
  TEUCHOS_ASSERT(pointA_->getDomainMap()->isSameAs(*blockMap_->pointMap_));

  {
    auto lcl_blockA           = blockA_->getLocalMatrixHost();
    auto lcl_pointA           = pointA_->getLocalMatrixHost();
    auto lcl_blockSizesRowMap = blockMap_->blockSizes_->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto lcl_offsets          = Kokkos::create_mirror_view(blockMap_->offsets_);
    Kokkos::deep_copy(lcl_offsets, blockMap_->offsets_);
    typename lo_vec_type::dual_view_type::t_host::const_type lcl_blockSizesColMap;
    Kokkos::View<size_t*, Kokkos::HostSpace> lcl_ghostedOffsets;
    if (!ghosted_blockMap.is_null()) {
      lcl_blockSizesColMap = ghosted_blockMap_->blockSizes_->getLocalViewHost(Tpetra::Access::ReadOnly);
      lcl_ghostedOffsets   = Kokkos::create_mirror_view(ghosted_blockMap_->offsets_);
      Kokkos::deep_copy(lcl_ghostedOffsets, ghosted_blockMap_->offsets_);
    } else {
      lcl_blockSizesColMap = lcl_blockSizesRowMap;
      lcl_ghostedOffsets   = lcl_offsets;
    }
    for (LocalOrdinal brlid = 0; brlid < lcl_blockA.numRows(); ++brlid) {
      size_t brsize = lcl_blockSizesRowMap(brlid, 0);
      auto brow     = lcl_blockA.row(brlid);
      for (LocalOrdinal k = 0; k < brow.length; ++k) {
        LocalOrdinal bclid = brow.colidx(k);
        size_t bcsize      = lcl_blockSizesColMap(bclid, 0);

        // We expect a dense block of size brsize*bcsize.
        const LocalOrdinal row_start = lcl_offsets(brlid);
        const LocalOrdinal row_end   = lcl_offsets(brlid + 1);
        const LocalOrdinal col_start = lcl_ghostedOffsets(bclid);
        const LocalOrdinal col_end   = lcl_ghostedOffsets(bclid + 1);

        TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(row_end - row_start), brsize);
        TEUCHOS_ASSERT_EQUALITY(Teuchos::as<size_t>(col_end - col_start), bcsize);

        size_t entries = 0;
        for (LocalOrdinal rlid = row_start; rlid < row_end; ++rlid) {
          auto row            = lcl_pointA.row(rlid);
          size_t entriesInRow = 0;
          for (LocalOrdinal n = 0; n < row.length; ++n) {
            auto clid = row.colidx(n);
            if ((col_start <= clid) && (clid < col_end)) {
              ++entriesInRow;
              ++entries;
            }
          }
          TEUCHOS_ASSERT_EQUALITY(entriesInRow, bcsize);
        }
        TEUCHOS_ASSERT_EQUALITY(entries, brsize * bcsize);
      }
    }
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
          Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const {
  pointA_->apply(X, Y, mode, alpha, beta);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    localApply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
               Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
               Teuchos::ETransp mode,
               Scalar alpha,
               Scalar beta) const {
  pointA_->localApply(X, Y, mode, alpha, beta);
}
}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMATRIX_DEF_HPP
