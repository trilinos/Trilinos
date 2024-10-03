// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMAP_DECL_HPP
#define TPETRA_BLOCKEDMAP_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Kokkos_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace Tpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class BlockedMap {
 public:
  using map_type    = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using lo_vec_type = Tpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

  BlockedMap(const Teuchos::RCP<const map_type>& pointMap,
             const Teuchos::RCP<lo_vec_type>& blockSizes);

  // private:
  Teuchos::RCP<const map_type> pointMap_;
  Teuchos::RCP<const map_type> blockMap_;
  Teuchos::RCP<lo_vec_type> blockSizes_;
  Kokkos::View<size_t*> offsets_;
  LocalOrdinal minClusterSize_;
  LocalOrdinal maxClusterSize_;
};

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMAP_DECL_HPP
