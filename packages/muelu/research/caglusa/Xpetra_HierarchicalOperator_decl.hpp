// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_HIERARCHICALOPERATOR_DECL_HPP
#define XPETRA_HIERARCHICALOPERATOR_DECL_HPP

#include <Tpetra_HierarchicalOperator_decl.hpp>
#include <Tpetra_HierarchicalOperator_def.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_TpetraBlockedMatrix.hpp>
#include "Xpetra_OperatorWithDiagonal.hpp"
#include "Xpetra_RowMatrix.hpp"

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class HierarchicalOperator : public OperatorWithDiagonal<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using tHOp                = Tpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using map_type            = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using mv_type             = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using matrix_type         = Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using blocked_matrix_type = Xpetra::TpetraBlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //! @name Constructor/Destructor
  //@{

  //! Constructor
  HierarchicalOperator(const RCP<tHOp>& op)
    : op_(op) {
    this->setTpetra_RowMatrix(op);
  }

  HierarchicalOperator(const RCP<matrix_type>& nearField,
                       const RCP<blocked_matrix_type>& kernelApproximations,
                       const RCP<matrix_type>& basisMatrix,
                       std::vector<RCP<blocked_matrix_type> >& transferMatrices,
                       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const map_type> getDomainMap() const {
    return toXpetra(op_->getDomainMap());
  }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const map_type> getRangeMap() const {
    return toXpetra(op_->getRangeMap());
  }

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
    vary according to the values of \c alpha and \c beta. Specifically
    - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
    - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
  */
  void apply(const mv_type& X, mv_type& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {
    op_->apply(Xpetra::toTpetra(X), Xpetra::toTpetra(Y), mode, alpha, beta);
  }

  //! Compute a residual R = B - (*this) * X
  void residual(const mv_type& X,
                const mv_type& B,
                mv_type& R) const {
    Tpetra::Details::residual(*op_, toTpetra(X), toTpetra(B), toTpetra(R));
  }

  RCP<HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > restrict(const RCP<matrix_type>& P) {
    using TpCrs   = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using CrsWrap = CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    return Teuchos::rcp(new HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op_->restrict(Teuchos::rcp_dynamic_cast<TpCrs>(Teuchos::rcp_dynamic_cast<CrsWrap>(P)->getCrsMatrix(), true)->getTpetra_CrsMatrixNonConst())));
  }

  RCP<matrix_type> toMatrix() {
    using TpCrs   = TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using CrsWrap = CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    auto tpMat    = Teuchos::rcp(new TpCrs(op_->toMatrix()));
    return Teuchos::rcp(new CrsWrap(Teuchos::rcp_dynamic_cast<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpMat)));
  }

  double getCompression() {
    return op_->getCompression();
  }

  void setDebugOutput(const bool debugOutput) {
    op_->setDebugOutput(debugOutput);
  }

  RCP<matrix_type> nearFieldMatrix() {
    auto tpMat = Teuchos::rcp(new TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(op_->nearFieldMatrix()));
    return Teuchos::rcp(new CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Teuchos::rcp_dynamic_cast<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpMat)));
  }

  //! Gets the operator out
  // RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_RowMatrix() const override { return op_; }

  // RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_RowMatrixNonConst() const override { return op_; }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
    describe(out, verbLevel, true);
  }

  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel, const bool printHeader) const {
    op_->describe(out, verbLevel, printHeader);
  }

  bool hasFarField() const {
    return op_->hasFarField();
  }

  bool denserThanDenseMatrix() const {
    return op_->denserThanDenseMatrix();
  }

 private:
  RCP<tHOp> op_;
};

}  // namespace Xpetra

#endif  // XPETRA_HIERARCHICALOPERATOR_DECL_HPP
