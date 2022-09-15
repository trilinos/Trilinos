#ifndef TPETRA_HIERARCHICALOPERATOR_DECL_HPP
#define TPETRA_HIERARCHICALOPERATOR_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_BlockedMatrix_decl.hpp>


namespace Tpetra {

    /*
    The unknowns of the kernel approximations are collected in the clusterMap.
    For H-matrices, this is just a concatenation.
    For H2-matrices, the map also contains the intermediate clusters that might be needed in upward/downward pass.


    H = nearField
        + basisMatrix *
          ((I+transferMatrices[K-1]^T) * ... * (I+transferMatrices[0]^T)) *
          kernelApproximations *
          ((I+transferMatrices[0]) * ... * (I+transferMatrices[K-1])) *
          basisMatrix^T

   nearField and basisMatrix are CRS matrices.
   kernelApproximations and transferMatrices[.] are blocked CRS matrices

   Maps:
   map (standard): domain and range of H;
                   domain, range, row of nearField
   clusterMap (blocked map): domain, range, row, column of transferMatrices;
                             domain, range, row of kernelApproximations
   ghosted_clusterMap (blocked map): column of kernelApproximations


   For H-matrices:
   K = 0, i.e. no transfer matrices

   For H2-matrices:
   upward and downward pass in the cluster hierarchy are encoded in transfer matrices
   */

  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class HierarchicalOperator : public Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:
    using matrix_type = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using mv_type = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
    using map_type = Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>;

    //! The RowMatrix representing the base class of CrsMatrix
    using row_matrix_type = RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    using impl_scalar_type = typename row_matrix_type::impl_scalar_type;
    using mag_type = typename Kokkos::ArithTraits<impl_scalar_type>::mag_type;

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

    using blocked_matrix_type = BlockedMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using blocked_map_type = BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    HierarchicalOperator(const Teuchos::RCP<matrix_type>& nearField,
                         const Teuchos::RCP<blocked_matrix_type>& kernelApproximations,
                         const Teuchos::RCP<matrix_type>& basisMatrix,
                         std::vector<Teuchos::RCP<blocked_matrix_type> >& transferMatrices);

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      return nearField_->getDomainMap();
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      return nearField_->getRangeMap();
    }

    //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
    /*!
      \param[in]  X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
               Scalar beta  = Teuchos::ScalarTraits<Scalar>::zero()) const;


    Teuchos::RCP<HierarchicalOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > restrict(const Teuchos::RCP<matrix_type>& P);

    Teuchos::RCP<matrix_type> toMatrix();

    double getCompression() {
      size_t nnz = (nearField_->getGlobalNumEntries() +
                    kernelApproximations_->pointA_->getGlobalNumEntries() +
                    basisMatrix_->getGlobalNumEntries());
      for (size_t i = 0; i < transferMatrices_.size(); i++)
        nnz += transferMatrices_[i]->pointA_->getGlobalNumEntries();
      return Teuchos::as<double>(nnz) / (getDomainMap()->getGlobalNumElements()*getDomainMap()->getGlobalNumElements());
    }

    Teuchos::RCP<matrix_type> nearFieldMatrix() {
      return nearField_;
    }

    // Fake RowMatrix interface
    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const {
      return nearField_->getRowMap();
    }

    Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const {
      return nearField_->getColMap();
    }

    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
      return nearField_->getDomainMap()->getComm();
    }

    Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const {
      return nearField_->getCrsGraph();
    }

    global_size_t getGlobalNumRows() const {
      return nearField_->getGlobalNumRows();
    }

    global_size_t getGlobalNumCols() const {
      return nearField_->getGlobalNumCols();
    }

    size_t getLocalNumRows() const {
      return nearField_->getLocalNumRows();
    }

    size_t getLocalNumCols() const {
      return nearField_->getLocalNumCols();
    }

    GlobalOrdinal getIndexBase() const {
      return nearField_->getIndexBase();
    }

    global_size_t getGlobalNumEntries() const {
      return nearField_->getGlobalNumEntries();
    }

    size_t getLocalNumEntries() const {
      return nearField_->getLocalNumEntries();
    }

    size_t getNumEntriesInGlobalRow (GlobalOrdinal globalRow) const {
      throw std::runtime_error("Not implemented.");
    }

    size_t getNumEntriesInLocalRow (LocalOrdinal localRow) const {
      throw std::runtime_error("Not implemented.");
    }

    size_t getGlobalMaxNumRowEntries () const {
      throw std::runtime_error("Not implemented.");
    }

    size_t getLocalMaxNumRowEntries () const {
      throw std::runtime_error("Not implemented.");
    }

    bool hasColMap () const {
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
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      nonconst_global_inds_host_view_type &Indices,
                      nonconst_values_host_view_type &Values,
                      size_t& NumEntries) const {
      throw std::runtime_error("Not implemented.");
    }

    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     nonconst_local_inds_host_view_type &Indices,
                     nonconst_values_host_view_type &Values,
                     size_t& NumEntries) const {
      throw std::runtime_error("Not implemented.");
    }

    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      global_inds_host_view_type &indices,
                      values_host_view_type &values) const {
      throw std::runtime_error("Not implemented.");
    }

    void
    getLocalRowView (LocalOrdinal LocalRow,
                     local_inds_host_view_type & indices,
                     values_host_view_type & values) const {
      throw std::runtime_error("Not implemented.");
    }

    void getLocalDiagCopy (Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
      nearField_->getLocalDiagCopy(diag);
    }

    void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw std::runtime_error("Not implemented.");
    }

    void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
      throw std::runtime_error("Not implemented.");
    }

    mag_type getFrobeniusNorm() const {
      return 0.;
    }

    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
      using std::setw;
      using std::endl;
      const size_t numRows = nearField_->getRowMap()->getGlobalNumElements();
      const size_t nnzNearField = nearField_->getGlobalNumEntries();
      const double nnzNearPerRow = Teuchos::as<double>(nnzNearField)/numRows;
      const size_t nnzKernelApprox = kernelApproximations_->pointA_->getGlobalNumEntries();
      const size_t numClusterPairs = kernelApproximations_->blockA_->getGlobalNumEntries();
      const size_t nnzBasis = basisMatrix_->getGlobalNumEntries();
      size_t nnzTransfer = 0;
      for (size_t i = 0; i<transferMatrices_.size(); i++)
        nnzTransfer += transferMatrices_[i]->pointA_->getGlobalNumEntries();
      const size_t nnzTotal = nnzNearField+nnzKernelApprox+nnzBasis+nnzTransfer;
      const double nnzTotalPerRow = Teuchos::as<double>(nnzTotal)/numRows;
      std::ostringstream oss;
      oss << std::left;
      oss << setw(9) << "rows"  << setw(12) << "nnz(near)"  << setw(14) << "nnz(near)/row" << setw(12) << "nnz(basis)" << setw(15) << "#cluster pairs" << setw(12)<< "nnz(kernel)"    << setw(14) << "nnz(transfer)" << setw(12) << "nnz(total)" << setw(14) << "nnz(total)/row" << endl;
      oss << setw(9) << numRows << setw(12) << nnzNearField << setw(14) << nnzNearPerRow   << setw(12) << nnzBasis     << setw(15) << numClusterPairs  << setw(12) << nnzKernelApprox << setw(14) << nnzTransfer     << setw(12) << nnzTotal     << setw(14) << nnzTotalPerRow   << endl;
      out << oss.str();
    }

  private:

    void allocateMemory(size_t numVectors) const;

    Teuchos::RCP<matrix_type> nearField_;
    Teuchos::RCP<blocked_matrix_type> kernelApproximations_;
    Teuchos::RCP<matrix_type> basisMatrix_;
    std::vector<Teuchos::RCP<blocked_matrix_type> > transferMatrices_;
    Teuchos::RCP<const map_type> clusterCoeffMap_;
    mutable Teuchos::RCP<mv_type> coefficients_, coefficients2_;
    mutable Teuchos::RCP<mv_type> X_colmap_, coefficients_colmap_;
  };
}

#endif // TPETRA_HIERARCHICALOPERATOR_DECL_HPP
