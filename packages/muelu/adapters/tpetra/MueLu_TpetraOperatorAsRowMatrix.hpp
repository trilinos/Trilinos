// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_RowMatrix.hpp>

namespace MueLu {

template <class Scalar        = Tpetra::Operator<>::scalar_type,
          class LocalOrdinal  = typename Tpetra::Operator<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class TpetraOperatorAsRowMatrix : public Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using op_type  = Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using vec_type = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //! The RowMatrix representing the base class of CrsMatrix
  using row_matrix_type = Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
  using mag_type         = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

  using local_inds_device_view_type =
      typename row_matrix_type::local_inds_device_view_type;
  using local_inds_host_view_type =
      typename row_matrix_type::local_inds_host_view_type;
  using nonconst_local_inds_host_view_type =
      typename row_matrix_type::nonconst_local_inds_host_view_type;

  using global_inds_device_view_type =
      typename row_matrix_type::global_inds_device_view_type;
  using global_inds_host_view_type =
      typename row_matrix_type::global_inds_host_view_type;
  using nonconst_global_inds_host_view_type =
      typename row_matrix_type::nonconst_global_inds_host_view_type;

  using values_device_view_type =
      typename row_matrix_type::values_device_view_type;
  using values_host_view_type =
      typename row_matrix_type::values_host_view_type;
  using nonconst_values_host_view_type =
      typename row_matrix_type::nonconst_values_host_view_type;

  //! @name Constructor/Destructor
  //@{

  //! Constructor
  TpetraOperatorAsRowMatrix(const RCP<op_type>& op)
    : op_(op)
    , diag_(Teuchos::null) {}

  TpetraOperatorAsRowMatrix(const RCP<op_type>& op,
                            const RCP<vec_type>& diag)
    : op_(op)
    , diag_(diag) {}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const {
    return op_->getDomainMap();
  }

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const {
    return op_->getRangeMap();
  }

  //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
  /*!
    \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::zero()) const {
    op_->apply(X, Y, mode, alpha, beta);
  }

  // Fake RowMatrix interface
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const {
    return op_->getRangeMap();
  }

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getColMap() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  typename row_matrix_type::local_ordinal_type getBlockSize() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
    return op_->getDomainMap()->getComm();
  }

  Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal, GlobalOrdinal, Node> > getGraph() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  Tpetra::global_size_t getGlobalNumRows() const {
    return getRowMap()->getGlobalNumElements();
  }

  Tpetra::global_size_t getGlobalNumCols() const {
    return getDomainMap()->getGlobalNumElements();
  }

  size_t getLocalNumRows() const {
    return getRowMap()->getLocalNumElements();
  }

  size_t getLocalNumCols() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  GlobalOrdinal getIndexBase() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  Tpetra::global_size_t getGlobalNumEntries() const {
    return 0;
  }

  size_t getLocalNumEntries() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  size_t getGlobalMaxNumRowEntries() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  size_t getLocalMaxNumRowEntries() const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  bool hasColMap() const {
    return false;
  }

  bool isLocallyIndexed() const {
    return true;
  }

  bool isGloballyIndexed() const {
    return true;
  }

  bool isFillComplete() const {
    return true;
  }

  bool supportsRowViews() const {
    return false;
  }

  void
  getGlobalRowCopy(GlobalOrdinal GlobalRow,
                   nonconst_global_inds_host_view_type& Indices,
                   nonconst_values_host_view_type& Values,
                   size_t& NumEntries) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  void
  getLocalRowCopy(LocalOrdinal LocalRow,
                  nonconst_local_inds_host_view_type& Indices,
                  nonconst_values_host_view_type& Values,
                  size_t& NumEntries) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  void
  getGlobalRowView(GlobalOrdinal GlobalRow,
                   global_inds_host_view_type& indices,
                   values_host_view_type& values) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  void
  getLocalRowView(LocalOrdinal LocalRow,
                  local_inds_host_view_type& indices,
                  values_host_view_type& values) const {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  void getLocalDiagCopy(Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& diag) const {
    if (diag_.is_null())
      throw MueLu::Exceptions::RuntimeError("No diagonal available.");
    else
      diag = *diag_;
  }

  void leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  void rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
    throw MueLu::Exceptions::RuntimeError("Not implemented.");
  }

  mag_type getFrobeniusNorm() const {
    return 0.;
  }

  // void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  //   using std::setw;
  //   using std::endl;
  //   const size_t numRows = nearField_->getRowMap()->getGlobalNumElements();
  //   const size_t nnzNearField = nearField_->getGlobalNumEntries();
  //   const double nnzNearPerRow = Teuchos::as<double>(nnzNearField)/numRows;
  //   const size_t nnzKernelApprox = kernelApproximations_->pointA_->getGlobalNumEntries();
  //   const size_t numClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
  //   const size_t nnzBasis = basisMatrix_->getGlobalNumEntries();
  //   size_t nnzTransfer = 0;
  //   for (size_t i = 0; i<transferMatrices_.size(); i++)
  //     nnzTransfer += transferMatrices_[i]->pointA_->getGlobalNumEntries();
  //   const size_t nnzTotal = nnzNearField+nnzKernelApprox+nnzBasis+nnzTransfer;
  //   const double nnzTotalPerRow = Teuchos::as<double>(nnzTotal)/numRows;
  //   std::ostringstream oss;
  //   oss << std::left;
  //   oss << setw(9) << "rows"  << setw(12) << "nnz(near)"  << setw(14) << "nnz(near)/row" << setw(12) << "nnz(basis)" << setw(15) << "#cluster pairs" << setw(12)<< "nnz(kernel)"    << setw(14) << "nnz(transfer)" << setw(12) << "nnz(total)" << setw(14) << "nnz(total)/row" << endl;
  //   oss << setw(9) << numRows << setw(12) << nnzNearField << setw(14) << nnzNearPerRow   << setw(12) << nnzBasis     << setw(15) << numClusterPairs  << setw(12) << nnzKernelApprox << setw(14) << nnzTransfer     << setw(12) << nnzTotal     << setw(14) << nnzTotalPerRow   << endl;
  //   out << oss.str();
  // }

 private:
  RCP<op_type> op_;
  RCP<vec_type> diag_;
};

}  // namespace MueLu
