#ifndef XPETRA_TPETRABLOCKEDMATRIX_HPP
#define XPETRA_TPETRABLOCKEDMATRIX_HPP

#include <Xpetra_TpetraBlockedMap.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

#include <Tpetra_BlockedMatrix_decl.hpp>
#include <Tpetra_BlockedMatrix_def.hpp>

namespace Xpetra {

template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class TpetraBlockedMatrix {
 public:
  using matrix_type                = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using blocked_map_type           = TpetraBlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using tpetra_blocked_matrix_type = Tpetra::BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  TpetraBlockedMatrix(const Teuchos::RCP<matrix_type>& pointA,
                      const Teuchos::RCP<matrix_type>& blockA,
                      const Teuchos::RCP<blocked_map_type>& blockMap,
                      const Teuchos::RCP<blocked_map_type>& ghosted_blockMap = Teuchos::null) {
    using TpCrs   = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using CrsWrap = CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    RCP<typename blocked_map_type::tpetra_blocked_map_type> tp_ghosted_blockMap;
    if (!ghosted_blockMap.is_null())
      tp_ghosted_blockMap = ghosted_blockMap->getTpetra_BlockedMap();
    blockMatrix_ = Teuchos::rcp(new tpetra_blocked_matrix_type(Teuchos::rcp_dynamic_cast<TpCrs>(Teuchos::rcp_dynamic_cast<CrsWrap>(pointA)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                                                               Teuchos::rcp_dynamic_cast<TpCrs>(Teuchos::rcp_dynamic_cast<CrsWrap>(blockA)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                                                               blockMap->getTpetra_BlockedMap(),
                                                               tp_ghosted_blockMap));
  }

  void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {
    blockMatrix_->apply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
  }

  void localApply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                  Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
                  Teuchos::ETransp mode = Teuchos::NO_TRANS,
                  Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
                  Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {
    blockMatrix_->localApply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
  }

  Teuchos::RCP<tpetra_blocked_matrix_type> getTpetra_BlockedMatrix() const {
    return blockMatrix_;
  }

 private:
  Teuchos::RCP<tpetra_blocked_matrix_type> blockMatrix_;
};

}  // namespace Xpetra

#endif  // XPETRA_TPETRABLOCKEDMATRIX_HPP
