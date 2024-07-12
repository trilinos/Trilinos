// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMATRIX_DECL_HPP
#define TPETRA_BLOCKEDMATRIX_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

namespace Tpetra {

template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class BlockedMatrix {
 public:
  using matrix_type      = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using blocked_map_type = BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using lo_vec_type      = typename blocked_map_type::lo_vec_type;

  BlockedMatrix(const Teuchos::RCP<matrix_type>& pointA,
                const Teuchos::RCP<matrix_type>& blockA,
                const Teuchos::RCP<blocked_map_type>& blockMap,
                const Teuchos::RCP<blocked_map_type>& ghosted_blockMap = Teuchos::null);

  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  void localApply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                  Teuchos::ETransp mode = Teuchos::NO_TRANS,
                  Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                  Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const;

  // private:
  Teuchos::RCP<matrix_type> pointA_;
  Teuchos::RCP<matrix_type> blockA_;
  Teuchos::RCP<blocked_map_type> blockMap_;
  Teuchos::RCP<blocked_map_type> ghosted_blockMap_;
};

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMATRIX_DECL_HPP
