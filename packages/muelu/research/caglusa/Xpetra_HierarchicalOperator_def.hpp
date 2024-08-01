// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_HIERARCHICALOPERATOR_DEF_HPP
#define XPETRA_HIERARCHICALOPERATOR_DEF_HPP

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    HierarchicalOperator(const Teuchos::RCP<matrix_type>& nearField,
                         const Teuchos::RCP<blocked_matrix_type>& kernelApproximations,
                         const Teuchos::RCP<matrix_type>& basisMatrix,
                         std::vector<Teuchos::RCP<blocked_matrix_type> >& transferMatrices,
                         const Teuchos::RCP<Teuchos::ParameterList>& params) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using TpCrs   = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using TpGOVec = TpetraMultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  using CrsWrap = CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  std::vector<RCP<typename tHOp::blocked_matrix_type> > tTransferMatrices;
  for (size_t i = 0; i < transferMatrices.size(); i++) {
    auto transferT = transferMatrices[i]->getTpetra_BlockedMatrix();
    tTransferMatrices.push_back(transferT);
  }

  op_ = rcp(new tHOp(Teuchos::rcp_dynamic_cast<TpCrs>(Teuchos::rcp_dynamic_cast<CrsWrap>(nearField)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                     kernelApproximations->getTpetra_BlockedMatrix(),
                     Teuchos::rcp_dynamic_cast<TpCrs>(Teuchos::rcp_dynamic_cast<CrsWrap>(basisMatrix)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst(),
                     tTransferMatrices,
                     params));
  this->setTpetra_RowMatrix(op_);
}

}  // namespace Xpetra

#endif  // XPETRA_HIERARCHICALOPERATOR_DEF_HPP
