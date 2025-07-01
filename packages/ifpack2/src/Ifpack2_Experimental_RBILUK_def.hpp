// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_EXPERIMENTAL_CRSRBILUK_DEF_HPP
#define IFPACK2_EXPERIMENTAL_CRSRBILUK_DEF_HPP

#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_BlockView.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers_decl.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"
#include "Ifpack2_Details_getCrsMatrix.hpp"
#include "Ifpack2_LocalFilter.hpp"
#include "Ifpack2_Utilities.hpp"
#include "Ifpack2_RILUK.hpp"
#include "KokkosSparse_sptrsv.hpp"

//#define IFPACK2_RBILUK_INITIAL
//#define IFPACK2_RBILUK_INITIAL_NOKK

#ifndef IFPACK2_RBILUK_INITIAL_NOKK
#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Util.hpp"
#endif

namespace Ifpack2 {

namespace Experimental {

namespace {
template <class MatrixType>
struct LocalRowHandler {
  using LocalOrdinal    = typename MatrixType::local_ordinal_type;
  using row_matrix_type = Tpetra::RowMatrix<
      typename MatrixType::scalar_type,
      LocalOrdinal,
      typename MatrixType::global_ordinal_type,
      typename MatrixType::node_type>;
  using inds_type = typename row_matrix_type::local_inds_host_view_type;
  using vals_type = typename row_matrix_type::values_host_view_type;

  LocalRowHandler(Teuchos::RCP<const row_matrix_type> A)
    : A_(std::move(A)) {
    if (!A_->supportsRowViews()) {
      const auto maxNumRowEntr = A_->getLocalMaxNumRowEntries();
      const auto blockSize     = A_->getBlockSize();
      ind_nc_                  = inds_type_nc("Ifpack2::RBILUK::LocalRowHandler::indices", maxNumRowEntr);
      val_nc_                  = vals_type_nc("Ifpack2::RBILUK::LocalRowHandler::values", maxNumRowEntr * blockSize * blockSize);
    }
  }

  void getLocalRow(LocalOrdinal local_row, inds_type& InI, vals_type& InV, LocalOrdinal& NumIn) {
    if (A_->supportsRowViews()) {
      A_->getLocalRowView(local_row, InI, InV);
      NumIn = (LocalOrdinal)InI.size();
    } else {
      size_t cnt;
      A_->getLocalRowCopy(local_row, ind_nc_, val_nc_, cnt);
      InI   = ind_nc_;
      InV   = val_nc_;
      NumIn = (LocalOrdinal)cnt;
    }
  }

 private:
  using inds_type_nc = typename row_matrix_type::nonconst_local_inds_host_view_type;
  using vals_type_nc = typename row_matrix_type::nonconst_values_host_view_type;

  Teuchos::RCP<const row_matrix_type> A_;
  inds_type_nc ind_nc_;
  vals_type_nc val_nc_;
};

template <typename BsrMatrixType, typename CrsMatrixType>
CrsMatrixType convertBsrToCrs(const BsrMatrixType& bsrMatrix) {
  using Ordinal     = typename BsrMatrixType::ordinal_type;
  using RowMapType  = typename BsrMatrixType::row_map_type::non_const_type;
  using EntriesType = typename BsrMatrixType::index_type::non_const_type;
  using ValuesType  = typename BsrMatrixType::values_type::non_const_type;
  using ExeSpace    = typename BsrMatrixType::execution_space;

  const Ordinal blockDim     = bsrMatrix.blockDim();
  const Ordinal blockSize    = blockDim * blockDim;
  const Ordinal numBlockRows = bsrMatrix.numRows();
  const Ordinal numBlockCols = bsrMatrix.numCols();

  const auto blockRowMap  = bsrMatrix.graph.row_map;
  const auto blockEntries = bsrMatrix.graph.entries;
  const auto blockValues  = bsrMatrix.values;

  const Ordinal numRows = numBlockRows * blockDim;
  const Ordinal numCols = numBlockCols * blockDim;

  // Allocate CrsMatrix row map
  RowMapType crsRowMap("crsRowMap", numRows + 1);
  EntriesType crsEntries("crsEntries", blockEntries.extent(0) * blockSize);
  ValuesType crsValues("crsValues", blockValues.extent(0) * blockSize);

  Kokkos::RangePolicy<ExeSpace> policy(0, numBlockRows);

  // Fill CrsMatrix row map, entries, and values
  Kokkos::parallel_for(
      "ConvertBsrToCrs", policy, KOKKOS_LAMBDA(const Ordinal blockRow) {
        const Ordinal blockRowStart = blockRowMap(blockRow);
        const Ordinal blockRowEnd   = blockRowMap(blockRow + 1);
        const Ordinal blockRowCount = blockRowEnd - blockRowStart;

        // Iterate over block entries in this row. This could use a teamRange parallel_for
        for (Ordinal blockNnz = blockRowStart; blockNnz < blockRowEnd; ++blockNnz) {
          const Ordinal blockCol = blockEntries(blockNnz);
          const Ordinal blockNum = blockNnz - blockRowStart;

          // Iterate over block dim to get the unblocked rows
          for (Ordinal blockRowOffset = 0; blockRowOffset < blockDim; ++blockRowOffset) {
            const Ordinal crsRow = blockRow * blockDim + blockRowOffset;
            // Each unblocked row has blockRowCount * blockDim items
            const Ordinal crsRowStart = blockRowStart * blockSize + blockRowCount * blockDim * blockRowOffset;
            crsRowMap(crsRow)         = crsRowStart;

            // Iterate over block dim to get the unblocked cols
            for (Ordinal blockColOffset = 0; blockColOffset < blockDim; ++blockColOffset) {
              const Ordinal crsCol = blockCol * blockDim + blockColOffset;
              const Ordinal crsNnz = crsRowStart + blockNum * blockDim + blockColOffset;
              crsEntries(crsNnz)   = crsCol;
              crsValues(crsNnz)    = blockValues(blockNnz * blockSize + blockRowOffset * blockDim + blockColOffset);
            }
          }
        }
      });

  // Finalize CrsMatrix row map
  Kokkos::RangePolicy<ExeSpace> policy2(0, 1);
  Kokkos::parallel_for(
      "FinalizeCrs", policy2, KOKKOS_LAMBDA(const Ordinal) {
        crsRowMap(numRows) = blockRowMap(numBlockRows) * blockSize;
      });

  // Construct CrsMatrix
  return CrsMatrixType("crsMatrix", numRows, numCols, crsEntries.extent(0), crsValues, crsRowMap, crsEntries);
}

}  // namespace

template <class MatrixType>
RBILUK<MatrixType>::RBILUK(const Teuchos::RCP<const row_matrix_type>& Matrix_in)
  : RILUK<row_matrix_type>(Teuchos::rcp_dynamic_cast<const row_matrix_type>(Matrix_in)) {}

template <class MatrixType>
RBILUK<MatrixType>::RBILUK(const Teuchos::RCP<const block_crs_matrix_type>& Matrix_in)
  : RILUK<row_matrix_type>(Teuchos::rcp_dynamic_cast<const row_matrix_type>(Matrix_in)) {}

template <class MatrixType>
RBILUK<MatrixType>::~RBILUK() {}

template <class MatrixType>
void RBILUK<MatrixType>::setMatrix(const Teuchos::RCP<const block_crs_matrix_type>& A) {
  // FIXME (mfh 04 Nov 2015) What about A_?  When does that get (re)set?

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  if (A.getRawPtr() != this->A_.getRawPtr()) {
    this->isAllocated_   = false;
    this->isInitialized_ = false;
    this->isComputed_    = false;
    this->Graph_         = Teuchos::null;
    L_block_             = Teuchos::null;
    U_block_             = Teuchos::null;
    D_block_             = Teuchos::null;
  }
}

template <class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getLBlock() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      L_block_.is_null(), std::runtime_error,
      "Ifpack2::RILUK::getL: The L factor "
      "is null.  Please call initialize() (and preferably also compute()) "
      "before calling this method.  If the input matrix has not yet been set, "
      "you must first call setMatrix() with a nonnull input matrix before you "
      "may call initialize() or compute().");
  return *L_block_;
}

template <class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getDBlock() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      D_block_.is_null(), std::runtime_error,
      "Ifpack2::RILUK::getD: The D factor "
      "(of diagonal entries) is null.  Please call initialize() (and "
      "preferably also compute()) before calling this method.  If the input "
      "matrix has not yet been set, you must first call setMatrix() with a "
      "nonnull input matrix before you may call initialize() or compute().");
  return *D_block_;
}

template <class MatrixType>
const typename RBILUK<MatrixType>::block_crs_matrix_type&
RBILUK<MatrixType>::getUBlock() const {
  TEUCHOS_TEST_FOR_EXCEPTION(
      U_block_.is_null(), std::runtime_error,
      "Ifpack2::RILUK::getU: The U factor "
      "is null.  Please call initialize() (and preferably also compute()) "
      "before calling this method.  If the input matrix has not yet been set, "
      "you must first call setMatrix() with a nonnull input matrix before you "
      "may call initialize() or compute().");
  return *U_block_;
}

template <class MatrixType>
void RBILUK<MatrixType>::allocate_L_and_U_blocks() {
  using Teuchos::null;
  using Teuchos::rcp;

  if (!this->isAllocated_) {
    if (!this->isKokkosKernelsStream_) {
      // Deallocate any existing storage.  This avoids storing 2x
      // memory, since RCP op= does not deallocate until after the
      // assignment.
      this->L_ = null;
      this->U_ = null;
      this->D_ = null;
      L_block_ = null;
      U_block_ = null;
      D_block_ = null;

      // Allocate Matrix using ILUK graphs
      L_block_ = rcp(new block_crs_matrix_type(*this->Graph_->getL_Graph(), blockSize_));
      U_block_ = rcp(new block_crs_matrix_type(*this->Graph_->getU_Graph(), blockSize_));
      D_block_ = rcp(new block_crs_matrix_type(*(Ifpack2::Details::computeDiagonalGraph(*(this->Graph_->getOverlapGraph()))),
                                               blockSize_));
      L_block_->setAllToScalar(STM::zero());  // Zero out L and U matrices
      U_block_->setAllToScalar(STM::zero());
      D_block_->setAllToScalar(STM::zero());

      // Allocate temp space for apply
      if (this->isKokkosKernelsSpiluk_) {
        const auto numRows = L_block_->getLocalNumRows();
        tmp_               = view1d("RBILUK::tmp_", numRows * blockSize_);
      }
    } else {
      this->D_ = null;
      D_block_ = null;

      L_block_v_ = std::vector<Teuchos::RCP<block_crs_matrix_type> >(this->num_streams_);
      U_block_v_ = std::vector<Teuchos::RCP<block_crs_matrix_type> >(this->num_streams_);
      for (int i = 0; i < this->num_streams_; i++) {
        L_block_v_[i] = rcp(new block_crs_matrix_type(*this->Graph_v_[i]->getL_Graph(), blockSize_));
        U_block_v_[i] = rcp(new block_crs_matrix_type(*this->Graph_v_[i]->getU_Graph(), blockSize_));
        L_block_v_[i]->setAllToScalar(STS::zero());  // Zero out L and U matrices
        U_block_v_[i]->setAllToScalar(STS::zero());
      }

      // Allocate temp space for apply
      auto lclMtx        = A_local_bcrs_->getLocalMatrixDevice();
      const auto numRows = lclMtx.numRows();
      tmp_               = view1d("RBILUK::tmp_", numRows * blockSize_);
    }
  }
  this->isAllocated_ = true;
}

namespace {

template <class MatrixType>
Teuchos::RCP<const typename RBILUK<MatrixType>::crs_graph_type>
getBlockCrsGraph(const Teuchos::RCP<const typename RBILUK<MatrixType>::row_matrix_type>& A) {
  using local_ordinal_type = typename MatrixType::local_ordinal_type;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;
  using row_matrix_type       = typename RBILUK<MatrixType>::row_matrix_type;
  using crs_graph_type        = typename RBILUK<MatrixType>::crs_graph_type;
  using block_crs_matrix_type = typename RBILUK<MatrixType>::block_crs_matrix_type;

  auto A_local = RBILUK<MatrixType>::makeLocalFilter(A);

  {
    RCP<const LocalFilter<row_matrix_type> > filteredA =
        rcp_dynamic_cast<const LocalFilter<row_matrix_type> >(A_local);
    RCP<const OverlappingRowMatrix<row_matrix_type> > overlappedA = Teuchos::null;
    RCP<const block_crs_matrix_type> A_block                      = Teuchos::null;
    if (!filteredA.is_null()) {
      overlappedA = rcp_dynamic_cast<const OverlappingRowMatrix<row_matrix_type> >(filteredA->getUnderlyingMatrix());
    }

    if (!overlappedA.is_null()) {
      A_block = rcp_dynamic_cast<const block_crs_matrix_type>(overlappedA->getUnderlyingMatrix());
    } else if (!filteredA.is_null()) {
      // If there is no overlap, filteredA could be the block CRS matrix
      A_block = rcp_dynamic_cast<const block_crs_matrix_type>(filteredA->getUnderlyingMatrix());
    } else {
      A_block = rcp_dynamic_cast<const block_crs_matrix_type>(A_local);
    }

    if (!A_block.is_null()) {
      return rcpFromRef(A_block->getCrsGraph());
    }
  }

  // Could not extract block crs, make graph manually
  {
    local_ordinal_type numRows = A_local->getLocalNumRows();
    Teuchos::Array<size_t> entriesPerRow(numRows);
    for (local_ordinal_type i = 0; i < numRows; i++) {
      entriesPerRow[i] = A_local->getNumEntriesInLocalRow(i);
    }
    RCP<crs_graph_type> A_local_crs_nc =
        rcp(new crs_graph_type(A_local->getRowMap(),
                               A_local->getColMap(),
                               entriesPerRow()));

    {
      using LocalRowHandler_t = LocalRowHandler<MatrixType>;
      LocalRowHandler_t localRowHandler(A_local);
      typename LocalRowHandler_t::inds_type indices;
      typename LocalRowHandler_t::vals_type values;
      for (local_ordinal_type i = 0; i < numRows; i++) {
        local_ordinal_type numEntries = 0;
        localRowHandler.getLocalRow(i, indices, values, numEntries);
        A_local_crs_nc->insertLocalIndices(i, numEntries, indices.data());
      }
    }

    A_local_crs_nc->fillComplete(A_local->getDomainMap(), A_local->getRangeMap());
    return rcp_const_cast<const crs_graph_type>(A_local_crs_nc);
  }
}

}  // namespace

template <class MatrixType>
void RBILUK<MatrixType>::initialize() {
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  using crs_map_type = Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>;

  const char prefix[] = "Ifpack2::Experimental::RBILUK::initialize: ";

  TEUCHOS_TEST_FOR_EXCEPTION(this->A_.is_null(), std::runtime_error, prefix << "The matrix (A_, the "
                                                                               "RowMatrix) is null.  Please call setMatrix() with a nonnull input "
                                                                               "before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(!this->A_->isFillComplete(), std::runtime_error, prefix << "The matrix "
                                                                                        "(A_, the BlockCrsMatrix) is not fill complete.  You may not invoke "
                                                                                        "initialize() or compute() with this matrix until the matrix is fill "
                                                                                        "complete.  Note: BlockCrsMatrix is fill complete if and only if its "
                                                                                        "underlying graph is fill complete.");

  blockSize_     = this->A_->getBlockSize();
  this->A_local_ = this->makeLocalFilter(this->A_);

  Teuchos::Time timer("RBILUK::initialize");
  double startTime = timer.wallTime();
  {  // Start timing
    Teuchos::TimeMonitor timeMon(timer);

    // Calling initialize() means that the user asserts that the graph
    // of the sparse matrix may have changed.  We must not just reuse
    // the previous graph in that case.
    //
    // Regarding setting isAllocated_ to false: Eventually, we may want
    // some kind of clever memory reuse strategy, but it's always
    // correct just to blow everything away and start over.
    this->isInitialized_ = false;
    this->isAllocated_   = false;
    this->isComputed_    = false;
    this->Graph_         = Teuchos::null;

    if (this->isKokkosKernelsSpiluk_) {
      A_local_bcrs_                          = Details::getBcrsMatrix(this->A_local_);
      RCP<const crs_matrix_type> A_local_crs = Details::getCrsMatrix(this->A_local_);
      if (A_local_bcrs_.is_null()) {
        if (A_local_crs.is_null()) {
          local_ordinal_type numRows = this->A_local_->getLocalNumRows();
          Array<size_t> entriesPerRow(numRows);
          for (local_ordinal_type i = 0; i < numRows; i++) {
            entriesPerRow[i] = this->A_local_->getNumEntriesInLocalRow(i);
          }
          A_local_crs_nc_ =
              rcp(new crs_matrix_type(this->A_local_->getRowMap(),
                                      this->A_local_->getColMap(),
                                      entriesPerRow()));
          // copy entries into A_local_crs
          nonconst_local_inds_host_view_type indices("indices", this->A_local_->getLocalMaxNumRowEntries());
          nonconst_values_host_view_type values("values", this->A_local_->getLocalMaxNumRowEntries());
          for (local_ordinal_type i = 0; i < numRows; i++) {
            size_t numEntries = 0;
            this->A_local_->getLocalRowCopy(i, indices, values, numEntries);
            A_local_crs_nc_->insertLocalValues(i, numEntries, reinterpret_cast<scalar_type*>(values.data()), indices.data());
          }
          A_local_crs_nc_->fillComplete(this->A_local_->getDomainMap(), this->A_local_->getRangeMap());
          A_local_crs = Teuchos::rcp_const_cast<const crs_matrix_type>(A_local_crs_nc_);
        }
        // Create bcrs from crs
        // We can skip fillLogicalBlocks if we know the input is already filled
        if (blockSize_ > 1) {
          auto crs_matrix_block_filled = Tpetra::fillLogicalBlocks(*A_local_crs, blockSize_);
          A_local_bcrs_                = Tpetra::convertToBlockCrsMatrix(*crs_matrix_block_filled, blockSize_);
        } else {
          A_local_bcrs_ = Tpetra::convertToBlockCrsMatrix(*A_local_crs, blockSize_);
        }
      }
    }

    if (this->isKokkosKernelsStream_) {
      std::vector<int> weights(this->num_streams_);
      std::fill(weights.begin(), weights.end(), 1);
      this->exec_space_instances_ = Kokkos::Experimental::partition_space(execution_space(), weights);

      auto lclMtx               = A_local_bcrs_->getLocalMatrixDevice();
      this->Graph_v_            = std::vector<Teuchos::RCP<Ifpack2::IlukGraph<crs_graph_type, kk_handle_type> > >(this->num_streams_);
      A_block_local_diagblks_v_ = std::vector<local_matrix_device_type>(this->num_streams_);
      std::vector<local_crs_matrix_device_type> A_crs_local_diagblks_v(this->num_streams_);

      if (!this->hasStreamReordered_) {
        auto lclCrsMtx = convertBsrToCrs<local_matrix_device_type, local_crs_matrix_device_type>(lclMtx);
        KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_sequential(lclCrsMtx, A_crs_local_diagblks_v);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, prefix << "Stream reordering is not currently supported for RBILUK");
      }

      for (int i = 0; i < this->num_streams_; i++) {
        A_block_local_diagblks_v_[i] = local_matrix_device_type(A_crs_local_diagblks_v[i], blockSize_);
      }

      A_block_local_diagblks_rowmap_v_  = std::vector<lno_row_view_t>(this->num_streams_);
      A_block_local_diagblks_entries_v_ = std::vector<lno_nonzero_view_t>(this->num_streams_);
      A_block_local_diagblks_values_v_  = std::vector<scalar_nonzero_view_t>(this->num_streams_);

      for (int i = 0; i < this->num_streams_; i++) {
        A_block_local_diagblks_rowmap_v_[i]  = A_block_local_diagblks_v_[i].graph.row_map;
        A_block_local_diagblks_entries_v_[i] = A_block_local_diagblks_v_[i].graph.entries;
        A_block_local_diagblks_values_v_[i]  = A_block_local_diagblks_v_[i].values;

        Teuchos::RCP<const crs_map_type> A_local_diagblks_RowMap = rcp(new crs_map_type(A_block_local_diagblks_v_[i].numRows(),
                                                                                        A_block_local_diagblks_v_[i].numRows(),
                                                                                        A_local_bcrs_->getRowMap()->getComm()));
        Teuchos::RCP<const crs_map_type> A_local_diagblks_ColMap = rcp(new crs_map_type(A_block_local_diagblks_v_[i].numCols(),
                                                                                        A_block_local_diagblks_v_[i].numCols(),
                                                                                        A_local_bcrs_->getColMap()->getComm()));

        Teuchos::RCP<const crs_graph_type> graph = rcp(new crs_graph_type(A_local_diagblks_RowMap, A_local_diagblks_ColMap, A_block_local_diagblks_v_[i].graph));
        this->Graph_v_[i]                        = rcp(new Ifpack2::IlukGraph<crs_graph_type, kk_handle_type>(graph,
                                                                                       this->LevelOfFill_, 0));
      }
    } else {
      RCP<const crs_graph_type> matrixCrsGraph = getBlockCrsGraph<MatrixType>(this->A_);
      this->Graph_                             = rcp(new Ifpack2::IlukGraph<crs_graph_type, kk_handle_type>(matrixCrsGraph,
                                                                                this->LevelOfFill_, 0));
    }

    if (this->isKokkosKernelsSpiluk_) {
      if (!this->isKokkosKernelsStream_) {
        KernelHandle_block_ = Teuchos::rcp(new kk_handle_type());
        const auto numRows  = this->A_local_->getLocalNumRows();
        KernelHandle_block_->create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1,
                                                  numRows,
                                                  2 * this->A_local_->getLocalNumEntries() * (this->LevelOfFill_ + 1),
                                                  2 * this->A_local_->getLocalNumEntries() * (this->LevelOfFill_ + 1),
                                                  blockSize_);
        this->Graph_->initialize(KernelHandle_block_);  // this calls spiluk_symbolic

        L_Sptrsv_KernelHandle_ = Teuchos::rcp(new kk_handle_type());
        U_Sptrsv_KernelHandle_ = Teuchos::rcp(new kk_handle_type());

        KokkosSparse::Experimental::SPTRSVAlgorithm alg = KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1;

        L_Sptrsv_KernelHandle_->create_sptrsv_handle(alg, numRows, true /*lower*/, blockSize_);
        U_Sptrsv_KernelHandle_->create_sptrsv_handle(alg, numRows, false /*upper*/, blockSize_);
      } else {
        KokkosSparse::Experimental::SPTRSVAlgorithm alg = KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1;
        KernelHandle_block_v_                           = std::vector<Teuchos::RCP<kk_handle_type> >(this->num_streams_);
        L_Sptrsv_KernelHandle_v_                        = std::vector<Teuchos::RCP<kk_handle_type> >(this->num_streams_);
        U_Sptrsv_KernelHandle_v_                        = std::vector<Teuchos::RCP<kk_handle_type> >(this->num_streams_);
        for (int i = 0; i < this->num_streams_; i++) {
          const auto numRows       = A_block_local_diagblks_v_[i].numRows();
          KernelHandle_block_v_[i] = Teuchos::rcp(new kk_handle_type());
          KernelHandle_block_v_[i]->create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm::SEQLVLSCHD_TP1,
                                                         numRows,
                                                         2 * A_block_local_diagblks_v_[i].nnz() * (this->LevelOfFill_ + 1),
                                                         2 * A_block_local_diagblks_v_[i].nnz() * (this->LevelOfFill_ + 1),
                                                         blockSize_);
          this->Graph_v_[i]->initialize(KernelHandle_block_v_[i]);  // this calls spiluk_symbolic

          L_Sptrsv_KernelHandle_v_[i] = Teuchos::rcp(new kk_handle_type());
          U_Sptrsv_KernelHandle_v_[i] = Teuchos::rcp(new kk_handle_type());

          L_Sptrsv_KernelHandle_v_[i]->create_sptrsv_handle(alg, numRows, true /*lower*/, blockSize_);
          U_Sptrsv_KernelHandle_v_[i]->create_sptrsv_handle(alg, numRows, false /*upper*/, blockSize_);
        }
      }
    } else {
      this->Graph_->initialize();
    }

    allocate_L_and_U_blocks();

#ifdef IFPACK2_RBILUK_INITIAL
    initAllValues();
#endif
  }  // Stop timing

  this->isInitialized_ = true;
  this->numInitialize_ += 1;
  this->initializeTime_ += (timer.wallTime() - startTime);
}

template <class MatrixType>
void RBILUK<MatrixType>::
    initAllValues() {
  using Teuchos::RCP;
  typedef Tpetra::Map<LO, GO, node_type> map_type;

  LO NumIn = 0, NumL = 0, NumU = 0;
  bool DiagFound         = false;
  size_t NumNonzeroDiags = 0;
  size_t MaxNumEntries   = this->A_->getLocalMaxNumRowEntries();
  LO blockMatSize        = blockSize_ * blockSize_;

  // First check that the local row map ordering is the same as the local portion of the column map.
  // The extraction of the strictly lower/upper parts of A, as well as the factorization,
  // implicitly assume that this is the case.
  Teuchos::ArrayView<const GO> rowGIDs = this->A_->getRowMap()->getLocalElementList();
  Teuchos::ArrayView<const GO> colGIDs = this->A_->getColMap()->getLocalElementList();
  bool gidsAreConsistentlyOrdered      = true;
  GO indexOfInconsistentGID            = 0;
  for (GO i = 0; i < rowGIDs.size(); ++i) {
    if (rowGIDs[i] != colGIDs[i]) {
      gidsAreConsistentlyOrdered = false;
      indexOfInconsistentGID     = i;
      break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(gidsAreConsistentlyOrdered == false, std::runtime_error,
                             "The ordering of the local GIDs in the row and column maps is not the same"
                                 << std::endl
                                 << "at index " << indexOfInconsistentGID
                                 << ".  Consistency is required, as all calculations are done with"
                                 << std::endl
                                 << "local indexing.");

  // Allocate temporary space for extracting the strictly
  // lower and upper parts of the matrix A.
  Teuchos::Array<LO> LI(MaxNumEntries);
  Teuchos::Array<LO> UI(MaxNumEntries);
  Teuchos::Array<scalar_type> LV(MaxNumEntries * blockMatSize);
  Teuchos::Array<scalar_type> UV(MaxNumEntries * blockMatSize);

  Teuchos::Array<scalar_type> diagValues(blockMatSize);

  L_block_->setAllToScalar(STM::zero());  // Zero out L and U matrices
  U_block_->setAllToScalar(STM::zero());
  D_block_->setAllToScalar(STM::zero());  // Set diagonal values to zero

  // NOTE (mfh 27 May 2016) The factorization below occurs entirely on
  // host, so sync to host first.  The const_cast is unfortunate but
  // is our only option to make this correct.

  /*
  const_cast<block_crs_matrix_type&> (A).sync_host ();
  L_block_->sync_host ();
  U_block_->sync_host ();
  D_block_->sync_host ();
  // NOTE (mfh 27 May 2016) We're modifying L, U, and D on host.
  L_block_->modify_host ();
  U_block_->modify_host ();
  D_block_->modify_host ();
  */

  RCP<const map_type> rowMap = L_block_->getRowMap();

  // First we copy the user's matrix into L and U, regardless of fill level.
  // It is important to note that L and U are populated using local indices.
  // This means that if the row map GIDs are not monotonically increasing
  // (i.e., permuted or gappy), then the strictly lower (upper) part of the
  // matrix is not the one that you would get if you based L (U) on GIDs.
  // This is ok, as the *order* of the GIDs in the rowmap is a better
  // expression of the user's intent than the GIDs themselves.

  // TODO BMK: Revisit this fence when BlockCrsMatrix is refactored.
  Kokkos::fence();
  using LocalRowHandler_t = LocalRowHandler<MatrixType>;
  LocalRowHandler_t localRowHandler(this->A_);
  typename LocalRowHandler_t::inds_type InI;
  typename LocalRowHandler_t::vals_type InV;
  for (size_t myRow = 0; myRow < this->A_->getLocalNumRows(); ++myRow) {
    LO local_row = myRow;

    localRowHandler.getLocalRow(local_row, InI, InV, NumIn);

    // Split into L and U (we don't assume that indices are ordered).
    NumL      = 0;
    NumU      = 0;
    DiagFound = false;

    for (LO j = 0; j < NumIn; ++j) {
      const LO k           = InI[j];
      const LO blockOffset = blockMatSize * j;

      if (k == local_row) {
        DiagFound = true;
        // Store perturbed diagonal in Tpetra::Vector D_
        for (LO jj = 0; jj < blockMatSize; ++jj)
          diagValues[jj] = this->Rthresh_ * InV[blockOffset + jj] + IFPACK2_SGN(InV[blockOffset + jj]) * this->Athresh_;
        D_block_->replaceLocalValues(local_row, &InI[j], diagValues.getRawPtr(), 1);
      } else if (k < 0) {  // Out of range
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::runtime_error,
            "Ifpack2::RILUK::initAllValues: current "
            "GID k = "
                << k << " < 0.  I'm not sure why this is an error; it is "
                        "probably an artifact of the undocumented assumptions of the "
                        "original implementation (likely copied and pasted from Ifpack).  "
                        "Nevertheless, the code I found here insisted on this being an error "
                        "state, so I will throw an exception here.");
      } else if (k < local_row) {
        LI[NumL]              = k;
        const LO LBlockOffset = NumL * blockMatSize;
        for (LO jj = 0; jj < blockMatSize; ++jj)
          LV[LBlockOffset + jj] = InV[blockOffset + jj];
        NumL++;
      } else if (Teuchos::as<size_t>(k) <= rowMap->getLocalNumElements()) {
        UI[NumU]              = k;
        const LO UBlockOffset = NumU * blockMatSize;
        for (LO jj = 0; jj < blockMatSize; ++jj)
          UV[UBlockOffset + jj] = InV[blockOffset + jj];
        NumU++;
      }
    }

    // Check in things for this row of L and U

    if (DiagFound) {
      ++NumNonzeroDiags;
    } else {
      for (LO jj = 0; jj < blockSize_; ++jj)
        diagValues[jj * (blockSize_ + 1)] = this->Athresh_;
      D_block_->replaceLocalValues(local_row, &local_row, diagValues.getRawPtr(), 1);
    }

    if (NumL) {
      L_block_->replaceLocalValues(local_row, &LI[0], &LV[0], NumL);
    }

    if (NumU) {
      U_block_->replaceLocalValues(local_row, &UI[0], &UV[0], NumU);
    }
  }

  // NOTE (mfh 27 May 2016) Sync back to device, in case compute()
  // ever gets a device implementation.
  /*
  {
    typedef typename block_crs_matrix_type::device_type device_type;
    const_cast<block_crs_matrix_type&> (A).template sync<device_type> ();
    L_block_->template sync<device_type> ();
    U_block_->template sync<device_type> ();
    D_block_->template sync<device_type> ();
  }
  */
  this->isInitialized_ = true;
}

namespace {  // (anonymous)

// For a given Kokkos::View type, possibly unmanaged, get the
// corresponding managed Kokkos::View type.  This is handy for
// translating from little_block_type or little_host_vec_type (both
// possibly unmanaged) to their managed versions.
template <class LittleBlockType>
struct GetManagedView {
  static_assert(Kokkos::is_view<LittleBlockType>::value,
                "LittleBlockType must be a Kokkos::View.");
  typedef Kokkos::View<typename LittleBlockType::non_const_data_type,
                       typename LittleBlockType::array_layout,
                       typename LittleBlockType::device_type>
      managed_non_const_type;
  static_assert(static_cast<int>(managed_non_const_type::rank) ==
                    static_cast<int>(LittleBlockType::rank),
                "managed_non_const_type::rank != LittleBlock::rank.  "
                "Please report this bug to the Ifpack2 developers.");
};

}  // namespace

template <class MatrixType>
void RBILUK<MatrixType>::compute() {
  using Teuchos::Array;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef impl_scalar_type IST;
  const char prefix[] = "Ifpack2::Experimental::RBILUK::compute: ";

  // initialize() checks this too, but it's easier for users if the
  // error shows them the name of the method that they actually
  // called, rather than the name of some internally called method.
  TEUCHOS_TEST_FOR_EXCEPTION(this->A_.is_null(), std::runtime_error, prefix << "The matrix (A_, "
                                                                               "the BlockCrsMatrix) is null.  Please call setMatrix() with a nonnull "
                                                                               "input before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(!this->A_->isFillComplete(), std::runtime_error, prefix << "The matrix "
                                                                                        "(A_, the BlockCrsMatrix) is not fill complete.  You may not invoke "
                                                                                        "initialize() or compute() with this matrix until the matrix is fill "
                                                                                        "complete.  Note: BlockCrsMatrix is fill complete if and only if its "
                                                                                        "underlying graph is fill complete.");

  if (!this->isInitialized()) {
    initialize();  // Don't count this in the compute() time
  }

  // // NOTE (mfh 27 May 2016) The factorization below occurs entirely on
  // // host, so sync to host first.  The const_cast is unfortunate but
  // // is our only option to make this correct.
  // if (! A_block_.is_null ()) {
  //   Teuchos::RCP<block_crs_matrix_type> A_nc =
  //     Teuchos::rcp_const_cast<block_crs_matrix_type> (A_block_);
  //   //    A_nc->sync_host ();
  // }
  /*
  L_block_->sync_host ();
  U_block_->sync_host ();
  D_block_->sync_host ();
  // NOTE (mfh 27 May 2016) We're modifying L, U, and D on host.
  L_block_->modify_host ();
  U_block_->modify_host ();
  D_block_->modify_host ();
  */

  Teuchos::Time timer("RBILUK::compute");
  double startTime = timer.wallTime();
  {  // Start timing
    Teuchos::TimeMonitor timeMon(timer);
    this->isComputed_ = false;

    // MinMachNum should be officially defined, for now pick something a little
    // bigger than IEEE underflow value

    //    const scalar_type MinDiagonalValue = STS::rmin ();
    //    const scalar_type MaxDiagonalValue = STS::one () / MinDiagonalValue;
    if (!this->isKokkosKernelsSpiluk_) {
      initAllValues();
      size_t NumIn;
      LO NumL, NumU, NumURead;

      TEUCHOS_TEST_FOR_EXCEPTION(
          this->isKokkosKernelsStream_, std::runtime_error,
          "Ifpack2::RBILUK::compute: "
          "streams are not yet supported without KokkosKernels.");

      // Get Maximum Row length
      const size_t MaxNumEntries =
          L_block_->getLocalMaxNumRowEntries() + U_block_->getLocalMaxNumRowEntries() + 1;

      const LO blockMatSize = blockSize_ * blockSize_;

      // FIXME (mfh 08 Nov 2015, 24 May 2016) We need to move away from
      // expressing these strides explicitly, in order to finish #177
      // (complete Kokkos-ization of BlockCrsMatrix) thoroughly.
      const LO rowStride = blockSize_;

      Teuchos::Array<int> ipiv_teuchos(blockSize_);
      Kokkos::View<int*, Kokkos::HostSpace,
                   Kokkos::MemoryUnmanaged>
          ipiv(ipiv_teuchos.getRawPtr(), blockSize_);
      Teuchos::Array<IST> work_teuchos(blockSize_);
      Kokkos::View<IST*, Kokkos::HostSpace,
                   Kokkos::MemoryUnmanaged>
          work(work_teuchos.getRawPtr(), blockSize_);

      size_t num_cols = U_block_->getColMap()->getLocalNumElements();
      Teuchos::Array<int> colflag(num_cols);

      typename GetManagedView<little_block_host_type>::managed_non_const_type diagModBlock("diagModBlock", blockSize_, blockSize_);
      typename GetManagedView<little_block_host_type>::managed_non_const_type matTmp("matTmp", blockSize_, blockSize_);
      typename GetManagedView<little_block_host_type>::managed_non_const_type multiplier("multiplier", blockSize_, blockSize_);

      //    Teuchos::ArrayRCP<scalar_type> DV = D_->get1dViewNonConst(); // Get view of diagonal

      // Now start the factorization.

      // Need some integer workspace and pointers
      LO NumUU;
      for (size_t j = 0; j < num_cols; ++j) {
        colflag[j] = -1;
      }
      Teuchos::Array<LO> InI(MaxNumEntries, 0);
      Teuchos::Array<scalar_type> InV(MaxNumEntries * blockMatSize, STM::zero());

      const LO numLocalRows = L_block_->getLocalNumRows();
      for (LO local_row = 0; local_row < numLocalRows; ++local_row) {
        // Fill InV, InI with current row of L, D and U combined

        NumIn = MaxNumEntries;
        local_inds_host_view_type colValsL;
        values_host_view_type valsL;

        L_block_->getLocalRowView(local_row, colValsL, valsL);
        NumL = (LO)colValsL.size();
        for (LO j = 0; j < NumL; ++j) {
          const LO matOffset = blockMatSize * j;
          little_block_host_type lmat((typename little_block_host_type::value_type*)&valsL[matOffset], blockSize_, rowStride);
          little_block_host_type lmatV((typename little_block_host_type::value_type*)&InV[matOffset], blockSize_, rowStride);
          // lmatV.assign(lmat);
          Tpetra::COPY(lmat, lmatV);
          InI[j] = colValsL[j];
        }

        little_block_host_type dmat = D_block_->getLocalBlockHostNonConst(local_row, local_row);
        little_block_host_type dmatV((typename little_block_host_type::value_type*)&InV[NumL * blockMatSize], blockSize_, rowStride);
        // dmatV.assign(dmat);
        Tpetra::COPY(dmat, dmatV);
        InI[NumL] = local_row;

        local_inds_host_view_type colValsU;
        values_host_view_type valsU;
        U_block_->getLocalRowView(local_row, colValsU, valsU);
        NumURead = (LO)colValsU.size();
        NumU     = 0;
        for (LO j = 0; j < NumURead; ++j) {
          if (!(colValsU[j] < numLocalRows)) continue;
          InI[NumL + 1 + j]  = colValsU[j];
          const LO matOffset = blockMatSize * (NumL + 1 + j);
          little_block_host_type umat((typename little_block_host_type::value_type*)&valsU[blockMatSize * j], blockSize_, rowStride);
          little_block_host_type umatV((typename little_block_host_type::value_type*)&InV[matOffset], blockSize_, rowStride);
          // umatV.assign(umat);
          Tpetra::COPY(umat, umatV);
          NumU += 1;
        }
        NumIn = NumL + NumU + 1;

        // Set column flags
        for (size_t j = 0; j < NumIn; ++j) {
          colflag[InI[j]] = j;
        }

#ifndef IFPACK2_RBILUK_INITIAL
        for (LO i = 0; i < blockSize_; ++i)
          for (LO j = 0; j < blockSize_; ++j) {
            {
              diagModBlock(i, j) = 0;
            }
          }
#else
        scalar_type diagmod = STM::zero();  // Off-diagonal accumulator
        Kokkos::deep_copy(diagModBlock, diagmod);
#endif

        for (LO jj = 0; jj < NumL; ++jj) {
          LO j = InI[jj];
          little_block_host_type currentVal((typename little_block_host_type::value_type*)&InV[jj * blockMatSize], blockSize_, rowStride);  // current_mults++;
          // multiplier.assign(currentVal);
          Tpetra::COPY(currentVal, multiplier);

          const little_block_host_type dmatInverse = D_block_->getLocalBlockHostNonConst(j, j);
          // alpha = 1, beta = 0
#ifndef IFPACK2_RBILUK_INITIAL_NOKK
          KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Algo::Gemm::Blocked>::invoke(STS::one(), currentVal, dmatInverse, STS::zero(), matTmp);
#else
          Tpetra::GEMM("N", "N", STS::one(), currentVal, dmatInverse,
                       STS::zero(), matTmp);
#endif
          // blockMatOpts.square_matrix_matrix_multiply(reinterpret_cast<impl_scalar_type*> (currentVal.data ()), reinterpret_cast<impl_scalar_type*> (dmatInverse.data ()), reinterpret_cast<impl_scalar_type*> (matTmp.data ()), blockSize_);
          // currentVal.assign(matTmp);
          Tpetra::COPY(matTmp, currentVal);
          local_inds_host_view_type UUI;
          values_host_view_type UUV;

          U_block_->getLocalRowView(j, UUI, UUV);
          NumUU = (LO)UUI.size();

          if (this->RelaxValue_ == STM::zero()) {
            for (LO k = 0; k < NumUU; ++k) {
              if (!(UUI[k] < numLocalRows)) continue;
              const int kk = colflag[UUI[k]];
              if (kk > -1) {
                little_block_host_type kkval((typename little_block_host_type::value_type*)&InV[kk * blockMatSize], blockSize_, rowStride);
                little_block_host_type uumat((typename little_block_host_type::value_type*)&UUV[k * blockMatSize], blockSize_, rowStride);
#ifndef IFPACK2_RBILUK_INITIAL_NOKK
                KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Algo::Gemm::Blocked>::invoke(magnitude_type(-STM::one()), multiplier, uumat, STM::one(), kkval);
#else
                Tpetra::GEMM("N", "N", magnitude_type(-STM::one()), multiplier, uumat,
                             STM::one(), kkval);
#endif
                // blockMatOpts.square_matrix_matrix_multiply(reinterpret_cast<impl_scalar_type*> (multiplier.data ()), reinterpret_cast<impl_scalar_type*> (uumat.data ()), reinterpret_cast<impl_scalar_type*> (kkval.data ()), blockSize_, -STM::one(), STM::one());
              }
            }
          } else {
            for (LO k = 0; k < NumUU; ++k) {
              if (!(UUI[k] < numLocalRows)) continue;
              const int kk = colflag[UUI[k]];
              little_block_host_type uumat((typename little_block_host_type::value_type*)&UUV[k * blockMatSize], blockSize_, rowStride);
              if (kk > -1) {
                little_block_host_type kkval((typename little_block_host_type::value_type*)&InV[kk * blockMatSize], blockSize_, rowStride);
#ifndef IFPACK2_RBILUK_INITIAL_NOKK
                KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Algo::Gemm::Blocked>::invoke(magnitude_type(-STM::one()), multiplier, uumat, STM::one(), kkval);
#else
                Tpetra::GEMM("N", "N", magnitude_type(-STM::one()), multiplier, uumat,
                             STM::one(), kkval);
#endif
                // blockMatOpts.square_matrix_matrix_multiply(reinterpret_cast<impl_scalar_type*>(multiplier.data ()), reinterpret_cast<impl_scalar_type*>(uumat.data ()), reinterpret_cast<impl_scalar_type*>(kkval.data ()), blockSize_, -STM::one(), STM::one());
              } else {
#ifndef IFPACK2_RBILUK_INITIAL_NOKK
                KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Trans::NoTranspose,
                                          KokkosBatched::Algo::Gemm::Blocked>::invoke(magnitude_type(-STM::one()), multiplier, uumat, STM::one(), diagModBlock);
#else
                Tpetra::GEMM("N", "N", magnitude_type(-STM::one()), multiplier, uumat,
                             STM::one(), diagModBlock);
#endif
                // blockMatOpts.square_matrix_matrix_multiply(reinterpret_cast<impl_scalar_type*>(multiplier.data ()), reinterpret_cast<impl_scalar_type*>(uumat.data ()), reinterpret_cast<impl_scalar_type*>(diagModBlock.data ()), blockSize_, -STM::one(), STM::one());
              }
            }
          }
        }
        if (NumL) {
          // Replace current row of L
          L_block_->replaceLocalValues(local_row, InI.getRawPtr(), InV.getRawPtr(), NumL);
        }

        // dmat.assign(dmatV);
        Tpetra::COPY(dmatV, dmat);

        if (this->RelaxValue_ != STM::zero()) {
          // dmat.update(this->RelaxValue_, diagModBlock);
          Tpetra::AXPY(this->RelaxValue_, diagModBlock, dmat);
        }

        //      if (STS::magnitude (DV[i]) > STS::magnitude (MaxDiagonalValue)) {
        //        if (STS::real (DV[i]) < STM::zero ()) {
        //          DV[i] = -MinDiagonalValue;
        //        }
        //        else {
        //          DV[i] = MinDiagonalValue;
        //        }
        //      }
        //      else
        {
          int lapackInfo = 0;
          for (int k = 0; k < blockSize_; ++k) {
            ipiv[k] = 0;
          }

          Tpetra::GETF2(dmat, ipiv, lapackInfo);
          // lapack.GETRF(blockSize_, blockSize_, d_raw, blockSize_, ipiv.getRawPtr(), &lapackInfo);
          TEUCHOS_TEST_FOR_EXCEPTION(
              lapackInfo != 0, std::runtime_error,
              "Ifpack2::Experimental::RBILUK::compute: "
              "lapackInfo = "
                  << lapackInfo << " which indicates an error in the factorization GETRF.");

          Tpetra::GETRI(dmat, ipiv, work, lapackInfo);
          // lapack.GETRI(blockSize_, d_raw, blockSize_, ipiv.getRawPtr(), work.getRawPtr(), lwork, &lapackInfo);
          TEUCHOS_TEST_FOR_EXCEPTION(
              lapackInfo != 0, std::runtime_error,
              "Ifpack2::Experimental::RBILUK::compute: "
              "lapackInfo = "
                  << lapackInfo << " which indicates an error in the matrix inverse GETRI.");
        }

        for (LO j = 0; j < NumU; ++j) {
          little_block_host_type currentVal((typename little_block_host_type::value_type*)&InV[(NumL + 1 + j) * blockMatSize], blockSize_, rowStride);  // current_mults++;
                                                                                                                                                        // scale U by the diagonal inverse
#ifndef IFPACK2_RBILUK_INITIAL_NOKK
          KokkosBatched::SerialGemm<KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Trans::NoTranspose,
                                    KokkosBatched::Algo::Gemm::Blocked>::invoke(STS::one(), dmat, currentVal, STS::zero(), matTmp);
#else
          Tpetra::GEMM("N", "N", STS::one(), dmat, currentVal,
                       STS::zero(), matTmp);
#endif
          // blockMatOpts.square_matrix_matrix_multiply(reinterpret_cast<impl_scalar_type*>(dmat.data ()), reinterpret_cast<impl_scalar_type*>(currentVal.data ()), reinterpret_cast<impl_scalar_type*>(matTmp.data ()), blockSize_);
          // currentVal.assign(matTmp);
          Tpetra::COPY(matTmp, currentVal);
        }

        if (NumU) {
          // Replace current row of L and U
          U_block_->replaceLocalValues(local_row, &InI[NumL + 1], &InV[blockMatSize * (NumL + 1)], NumU);
        }

#ifndef IFPACK2_RBILUK_INITIAL
        // Reset column flags
        for (size_t j = 0; j < NumIn; ++j) {
          colflag[InI[j]] = -1;
        }
#else
        // Reset column flags
        for (size_t j = 0; j < num_cols; ++j) {
          colflag[j] = -1;
        }
#endif
      }
    }  // !this->isKokkosKernelsSpiluk_
    else {
      // If we are not using A_local_ directly, we must copy values in case of pattern reuse. Due
      // to having to deal with filling logical blocks and converting back and forth between CRS and
      // BSR, this is not very performant and a lot of computation has to be repeated. If you are going
      // to do patten reuse, I recommend having A_local be a BSR.
      auto A_local_bcrs_temp = Details::getBcrsMatrix(this->A_local_);
      auto A_local_crs_temp  = Details::getCrsMatrix(this->A_local_);
      if (A_local_bcrs_temp.is_null()) {
        if (A_local_crs_temp.is_null()) {
          // copy entries into A_local_crs_nc
          local_ordinal_type numRows = this->A_local_->getLocalNumRows();
          A_local_crs_nc_->resumeFill();
          nonconst_local_inds_host_view_type indices("indices", this->A_local_->getLocalMaxNumRowEntries());
          nonconst_values_host_view_type values("values", this->A_local_->getLocalMaxNumRowEntries());
          for (local_ordinal_type i = 0; i < numRows; i++) {
            size_t numEntries = 0;
            this->A_local_->getLocalRowCopy(i, indices, values, numEntries);
            A_local_crs_nc_->insertLocalValues(i, numEntries, reinterpret_cast<scalar_type*>(values.data()), indices.data());
          }
          A_local_crs_nc_->fillComplete(this->A_local_->getDomainMap(), this->A_local_->getRangeMap());
        }

        if (blockSize_ > 1) {
          auto crs_matrix_block_filled = Tpetra::fillLogicalBlocks(*A_local_crs_nc_, blockSize_);
          Kokkos::deep_copy(A_local_bcrs_->getLocalMatrixDevice().values,
                            crs_matrix_block_filled->getLocalMatrixDevice().values);
        } else {
          Kokkos::deep_copy(A_local_bcrs_->getLocalMatrixDevice().values,
                            A_local_crs_nc_->getLocalMatrixDevice().values);
        }
      }

      if (!this->isKokkosKernelsStream_) {
        auto lclMtx          = A_local_bcrs_->getLocalMatrixDevice();
        auto A_local_rowmap  = lclMtx.graph.row_map;
        auto A_local_entries = lclMtx.graph.entries;
        auto A_local_values  = lclMtx.values;

        // L_block_->resumeFill ();
        // U_block_->resumeFill ();

        if (L_block_->isLocallyIndexed()) {
          L_block_->setAllToScalar(STS::zero());  // Zero out L and U matrices
          U_block_->setAllToScalar(STS::zero());
        }

        using row_map_type = typename local_matrix_device_type::row_map_type;

        auto lclL             = L_block_->getLocalMatrixDevice();
        row_map_type L_rowmap = lclL.graph.row_map;
        auto L_entries        = lclL.graph.entries;
        auto L_values         = lclL.values;

        auto lclU             = U_block_->getLocalMatrixDevice();
        row_map_type U_rowmap = lclU.graph.row_map;
        auto U_entries        = lclU.graph.entries;
        auto U_values         = lclU.values;

        KokkosSparse::spiluk_numeric(KernelHandle_block_.getRawPtr(), this->LevelOfFill_,
                                     A_local_rowmap, A_local_entries, A_local_values,
                                     L_rowmap, L_entries, L_values, U_rowmap, U_entries, U_values);

        // Now call symbolic for sptrsvs
        KokkosSparse::sptrsv_symbolic(L_Sptrsv_KernelHandle_.getRawPtr(), L_rowmap, L_entries, L_values);
        KokkosSparse::sptrsv_symbolic(U_Sptrsv_KernelHandle_.getRawPtr(), U_rowmap, U_entries, U_values);
      } else {
        // Due to A_block_local_diagblks_values_v_, we must always refresh when streams are on
        assert(!this->hasStreamReordered_);
        auto lclMtx    = A_local_bcrs_->getLocalMatrixDevice();
        auto lclCrsMtx = convertBsrToCrs<local_matrix_device_type, local_crs_matrix_device_type>(lclMtx);
        std::vector<local_crs_matrix_device_type> A_crs_local_diagblks_v(this->num_streams_);
        KokkosSparse::Impl::kk_extract_diagonal_blocks_crsmatrix_sequential(lclCrsMtx, A_crs_local_diagblks_v);

        for (int i = 0; i < this->num_streams_; i++) {
          A_block_local_diagblks_v_[i] = local_matrix_device_type(A_crs_local_diagblks_v[i], blockSize_);
          // rowmap and entries should be the same, but not values
          A_block_local_diagblks_values_v_[i] = A_block_local_diagblks_v_[i].values;
        }

        std::vector<lno_row_view_t> L_rowmap_v(this->num_streams_);
        std::vector<lno_nonzero_view_t> L_entries_v(this->num_streams_);
        std::vector<scalar_nonzero_view_t> L_values_v(this->num_streams_);
        std::vector<lno_row_view_t> U_rowmap_v(this->num_streams_);
        std::vector<lno_nonzero_view_t> U_entries_v(this->num_streams_);
        std::vector<scalar_nonzero_view_t> U_values_v(this->num_streams_);
        std::vector<kk_handle_type*> KernelHandle_rawptr_v_(this->num_streams_);

        for (int i = 0; i < this->num_streams_; i++) {
          L_block_v_[i]->setAllToScalar(STS::zero());  // Zero out L and U matrices
          U_block_v_[i]->setAllToScalar(STS::zero());

          auto lclL      = L_block_v_[i]->getLocalMatrixDevice();
          L_rowmap_v[i]  = lclL.graph.row_map;
          L_entries_v[i] = lclL.graph.entries;
          L_values_v[i]  = lclL.values;

          auto lclU                 = U_block_v_[i]->getLocalMatrixDevice();
          U_rowmap_v[i]             = lclU.graph.row_map;
          U_entries_v[i]            = lclU.graph.entries;
          U_values_v[i]             = lclU.values;
          KernelHandle_rawptr_v_[i] = KernelHandle_block_v_[i].getRawPtr();
        }

        // L_block_->resumeFill ();
        // U_block_->resumeFill ();

        KokkosSparse::Experimental::spiluk_numeric_streams(
            this->exec_space_instances_, KernelHandle_rawptr_v_, this->LevelOfFill_,
            A_block_local_diagblks_rowmap_v_, A_block_local_diagblks_entries_v_, A_block_local_diagblks_values_v_,
            L_rowmap_v, L_entries_v, L_values_v,
            U_rowmap_v, U_entries_v, U_values_v);

        // Now call symbolic for sptrsvs
        for (int i = 0; i < this->num_streams_; i++) {
          KokkosSparse::sptrsv_symbolic(L_Sptrsv_KernelHandle_v_[i].getRawPtr(), L_rowmap_v[i], L_entries_v[i], L_values_v[i]);
          KokkosSparse::sptrsv_symbolic(U_Sptrsv_KernelHandle_v_[i].getRawPtr(), U_rowmap_v[i], U_entries_v[i], U_values_v[i]);
        }
      }
    }
  }  // Stop timing

  // Sync everything back to device, for efficient solves.
  /*
  {
    typedef typename block_crs_matrix_type::device_type device_type;
    if (! A_block_.is_null ()) {
      Teuchos::RCP<block_crs_matrix_type> A_nc =
        Teuchos::rcp_const_cast<block_crs_matrix_type> (A_block_);
      A_nc->template sync<device_type> ();
    }
    L_block_->template sync<device_type> ();
    U_block_->template sync<device_type> ();
    D_block_->template sync<device_type> ();
  }
  */

  this->isComputed_ = true;
  this->numCompute_ += 1;
  this->computeTime_ += (timer.wallTime() - startTime);
}

template <class MatrixType>
template <typename X, typename Y>
void RBILUK<MatrixType>::
    stream_apply_impl(const X& X_views, Y& Y_views, const bool use_temp_x, const bool use_temp_y, const std::vector<Teuchos::RCP<block_crs_matrix_type> >& LU_block_v, const std::vector<Teuchos::RCP<kk_handle_type> >& LU_sptrsv_handle_v, const LO numVecs) const {
  std::vector<lno_row_view_t> ptr_v(this->num_streams_);
  std::vector<lno_nonzero_view_t> ind_v(this->num_streams_);
  std::vector<scalar_nonzero_view_t> val_v(this->num_streams_);
  std::vector<kk_handle_type*> KernelHandle_rawptr_v(this->num_streams_);

  for (int j = 0; j < numVecs; ++j) {
    auto X_view = Kokkos::subview(X_views, Kokkos::ALL(), j);
    auto Y_view = Kokkos::subview(Y_views, Kokkos::ALL(), j);
    std::vector<decltype(X_view)> x_v(this->num_streams_);
    std::vector<decltype(Y_view)> y_v(this->num_streams_);
    local_ordinal_type stream_begin = 0;
    local_ordinal_type stream_end;
    for (int i = 0; i < this->num_streams_; i++) {
      auto LU_bcrs_i = LU_block_v[i];
      auto LUlocal_i = LU_bcrs_i->getLocalMatrixDevice();
      ptr_v[i]       = LUlocal_i.graph.row_map;
      ind_v[i]       = LUlocal_i.graph.entries;
      val_v[i]       = LUlocal_i.values;
      stream_end     = stream_begin + LUlocal_i.numRows() * blockSize_;
      if (use_temp_x) {
        x_v[i] = Kokkos::subview(tmp_, Kokkos::make_pair(stream_begin, stream_end));
      } else {
        x_v[i] = Kokkos::subview(X_view, Kokkos::make_pair(stream_begin, stream_end));
      }
      if (use_temp_y) {
        y_v[i] = Kokkos::subview(tmp_, Kokkos::make_pair(stream_begin, stream_end));
      } else {
        y_v[i] = Kokkos::subview(Y_view, Kokkos::make_pair(stream_begin, stream_end));
      }
      KernelHandle_rawptr_v[i] = LU_sptrsv_handle_v[i].getRawPtr();
      stream_begin             = stream_end;
    }

    Kokkos::fence();
    KokkosSparse::Experimental::sptrsv_solve_streams(this->exec_space_instances_, KernelHandle_rawptr_v, ptr_v, ind_v, val_v, x_v, y_v);
    Kokkos::fence();
  }
}

template <class MatrixType>
void RBILUK<MatrixType>::
    apply(const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
          Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
          Teuchos::ETransp mode,
          scalar_type alpha,
          scalar_type beta) const {
  using Teuchos::RCP;
  typedef Tpetra::BlockMultiVector<scalar_type,
                                   local_ordinal_type, global_ordinal_type, node_type>
      BMV;

  TEUCHOS_TEST_FOR_EXCEPTION(
      this->A_.is_null(), std::runtime_error,
      "Ifpack2::Experimental::RBILUK::apply: The matrix is "
      "null.  Please call setMatrix() with a nonnull input, then initialize() "
      "and compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      !this->isComputed(), std::runtime_error,
      "Ifpack2::Experimental::RBILUK::apply: If you have not yet called compute(), "
      "you must call compute() before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
      "Ifpack2::Experimental::RBILUK::apply: X and Y do not have the same number of columns.  "
      "X.getNumVectors() = "
          << X.getNumVectors()
          << " != Y.getNumVectors() = " << Y.getNumVectors() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isComplex && mode == Teuchos::CONJ_TRANS, std::logic_error,
      "Ifpack2::Experimental::RBILUK::apply: mode = Teuchos::CONJ_TRANS is not implemented for "
      "complex Scalar type.  Please talk to the Ifpack2 developers to get this "
      "fixed.  There is a FIXME in this file about this very issue.");

  const LO blockMatSize = blockSize_ * blockSize_;

  const LO rowStride = blockSize_;

  Teuchos::Array<scalar_type> lclarray(blockSize_);
  little_host_vec_type lclvec((typename little_host_vec_type::value_type*)&lclarray[0], blockSize_);
  const scalar_type one  = STM::one();
  const scalar_type zero = STM::zero();

  Teuchos::Time timer("RBILUK::apply");
  double startTime = timer.wallTime();
  {  // Start timing
    Teuchos::TimeMonitor timeMon(timer);
    if (!this->isKokkosKernelsSpiluk_) {
      BMV yBlock(Y, *(this->Graph_->getA_Graph()->getDomainMap()), blockSize_);
      const BMV xBlock(X, *(this->A_->getColMap()), blockSize_);

      if (alpha == one && beta == zero) {
        if (mode == Teuchos::NO_TRANS) {  // Solve L (D (U Y)) = X for Y.
          // Start by solving L C = X for C.  C must have the same Map
          // as D.  We have to use a temp multivector, since our
          // implementation of triangular solves does not allow its
          // input and output to alias one another.
          //
          // FIXME (mfh 24 Jan 2014) Cache this temp multivector.
          const LO numVectors = xBlock.getNumVectors();
          BMV cBlock(*(this->Graph_->getA_Graph()->getDomainMap()), blockSize_, numVectors);
          BMV rBlock(*(this->Graph_->getA_Graph()->getDomainMap()), blockSize_, numVectors);
          for (LO imv = 0; imv < numVectors; ++imv) {
            for (size_t i = 0; i < D_block_->getLocalNumRows(); ++i) {
              LO local_row = i;
              const_host_little_vec_type xval =
                  xBlock.getLocalBlockHost(local_row, imv, Tpetra::Access::ReadOnly);
              little_host_vec_type cval =
                  cBlock.getLocalBlockHost(local_row, imv, Tpetra::Access::OverwriteAll);
              // cval.assign(xval);
              Tpetra::COPY(xval, cval);

              local_inds_host_view_type colValsL;
              values_host_view_type valsL;
              L_block_->getLocalRowView(local_row, colValsL, valsL);
              LO NumL = (LO)colValsL.size();

              for (LO j = 0; j < NumL; ++j) {
                LO col = colValsL[j];
                const_host_little_vec_type prevVal =
                    cBlock.getLocalBlockHost(col, imv, Tpetra::Access::ReadOnly);

                const LO matOffset = blockMatSize * j;
                little_block_host_type lij((typename little_block_host_type::value_type*)&valsL[matOffset], blockSize_, rowStride);

                // cval.matvecUpdate(-one, lij, prevVal);
                Tpetra::GEMV(-one, lij, prevVal, cval);
              }
            }
          }

          // Solve D R = C. Note that D has been replaced by D^{-1} at this point.
          D_block_->applyBlock(cBlock, rBlock);

          // Solve U Y = R.
          for (LO imv = 0; imv < numVectors; ++imv) {
            const LO numRows = D_block_->getLocalNumRows();
            for (LO i = 0; i < numRows; ++i) {
              LO local_row = (numRows - 1) - i;
              const_host_little_vec_type rval =
                  rBlock.getLocalBlockHost(local_row, imv, Tpetra::Access::ReadOnly);
              little_host_vec_type yval =
                  yBlock.getLocalBlockHost(local_row, imv, Tpetra::Access::OverwriteAll);
              // yval.assign(rval);
              Tpetra::COPY(rval, yval);

              local_inds_host_view_type colValsU;
              values_host_view_type valsU;
              U_block_->getLocalRowView(local_row, colValsU, valsU);
              LO NumU = (LO)colValsU.size();

              for (LO j = 0; j < NumU; ++j) {
                LO col = colValsU[NumU - 1 - j];
                const_host_little_vec_type prevVal =
                    yBlock.getLocalBlockHost(col, imv, Tpetra::Access::ReadOnly);

                const LO matOffset = blockMatSize * (NumU - 1 - j);
                little_block_host_type uij((typename little_block_host_type::value_type*)&valsU[matOffset], blockSize_, rowStride);

                // yval.matvecUpdate(-one, uij, prevVal);
                Tpetra::GEMV(-one, uij, prevVal, yval);
              }
            }
          }
        } else {  // Solve U^P (D^P (L^P Y)) = X for Y (where P is * or T).
          TEUCHOS_TEST_FOR_EXCEPTION(
              true, std::runtime_error,
              "Ifpack2::Experimental::RBILUK::apply: transpose apply is not implemented for the block algorithm");
        }
      } else {  // alpha != 1 or beta != 0
        if (alpha == zero) {
          if (beta == zero) {
            Y.putScalar(zero);
          } else {
            Y.scale(beta);
          }
        } else {  // alpha != zero
          MV Y_tmp(Y.getMap(), Y.getNumVectors());
          apply(X, Y_tmp, mode);
          Y.update(alpha, Y_tmp, beta);
        }
      }
    } else {
      // Kokkos kernels impl.
      auto X_views     = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
      auto Y_views     = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
      const LO numVecs = X.getNumVectors();

      TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS, std::runtime_error,
                                 "Ifpack2::Experimental::RBILUK::apply: transpose apply is not implemented for the block algorithm");

      if (!this->isKokkosKernelsStream_) {
        auto lclL      = L_block_->getLocalMatrixDevice();
        auto L_rowmap  = lclL.graph.row_map;
        auto L_entries = lclL.graph.entries;
        auto L_values  = lclL.values;

        auto lclU      = U_block_->getLocalMatrixDevice();
        auto U_rowmap  = lclU.graph.row_map;
        auto U_entries = lclU.graph.entries;
        auto U_values  = lclU.values;

        for (LO vec = 0; vec < numVecs; ++vec) {
          auto X_view = Kokkos::subview(X_views, Kokkos::ALL(), vec);
          KokkosSparse::sptrsv_solve(L_Sptrsv_KernelHandle_.getRawPtr(), L_rowmap, L_entries, L_values, X_view, tmp_);
        }

        for (LO vec = 0; vec < numVecs; ++vec) {
          auto Y_view = Kokkos::subview(Y_views, Kokkos::ALL(), vec);
          KokkosSparse::sptrsv_solve(U_Sptrsv_KernelHandle_.getRawPtr(), U_rowmap, U_entries, U_values, tmp_, Y_view);
        }
      } else {
        stream_apply_impl(X_views, Y_views, false, true, L_block_v_, L_Sptrsv_KernelHandle_v_, numVecs);
        stream_apply_impl(X_views, Y_views, true, false, U_block_v_, U_Sptrsv_KernelHandle_v_, numVecs);
      }

      KokkosBlas::axpby(alpha, Y_views, beta, Y_views);
      // Y.getWrappedDualView().sync();
    }
  }  // Stop timing

  this->numApply_ += 1;
  this->applyTime_ += (timer.wallTime() - startTime);
}

template <class MatrixType>
std::string RBILUK<MatrixType>::description() const {
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Experimental::RBILUK\": {";
  os << "Initialized: " << (this->isInitialized() ? "true" : "false") << ", "
     << "Computed: " << (this->isComputed() ? "true" : "false") << ", ";

  os << "Level-of-fill: " << this->getLevelOfFill() << ", ";

  if (this->A_.is_null()) {
    os << "Matrix: null";
  } else {
    os << "Global matrix dimensions: ["
       << this->A_->getGlobalNumRows() << ", " << this->A_->getGlobalNumCols() << "]"
       << ", Global nnz: " << this->A_->getGlobalNumEntries();
  }

  os << "}";
  return os.str();
}

}  // namespace Experimental

}  // namespace Ifpack2

// FIXME (mfh 26 Aug 2015) We only need to do instantiation for
// MatrixType = Tpetra::RowMatrix.  Conversions to BlockCrsMatrix are
// handled internally via dynamic cast.

#define IFPACK2_EXPERIMENTAL_RBILUK_INSTANT(S, LO, GO, N)                              \
  template class Ifpack2::Experimental::RBILUK<Tpetra::BlockCrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::Experimental::RBILUK<Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
