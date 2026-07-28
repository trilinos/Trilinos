// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef AMESOS2_CUSOLVER_DEF_HPP
#define AMESOS2_CUSOLVER_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_cuSOLVER_decl.hpp"

namespace Amesos2 {
namespace Impl {

// Standalone functor for scattering CSR entries into a dense column-major matrix.
// Must live outside the private member function to satisfy CUDA extended-lambda rules.
template<class MatrixView, class RowPtrView, class ColIndView, class NzView>
struct CsrToDenseFunctor {
  MatrixView matrix;
  RowPtrView row_ptr;
  ColIndView col_ind;
  NzView     nzvals;

  CsrToDenseFunctor(MatrixView m, RowPtrView r, ColIndView c, NzView v)
    : matrix(m), row_ptr(r), col_ind(c), nzvals(v) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int row) const {
    for(int j = row_ptr(row); j < row_ptr(row+1); ++j)
      matrix(row, col_ind(j)) = nzvals(j);
  }
};

// Functor to fill a square matrix with the identity.
template<class MatrixView>
struct SetIdentityFunctor {
  MatrixView matrix;
  SetIdentityFunctor(MatrixView m) : matrix(m) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    typedef typename MatrixView::value_type val_t;
    const int n = static_cast<int>(matrix.extent(1));
    for(int j = 0; j < n; ++j)
      matrix(i, j) = (i == j) ? val_t(1) : val_t(0);
  }
};

} // namespace Impl

template <class Matrix, class Vector>
cuSOLVER<Matrix,Vector>::cuSOLVER(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::cuSOLVER,Matrix,Vector>(A, X, B)
{
  auto status = cusolverDnCreate(&data_.dn_handle);
  TEUCHOS_TEST_FOR_EXCEPTION( status != CUSOLVER_STATUS_SUCCESS,
    std::runtime_error, "cusolverDnCreate failed");
  status = cusolverSpCreate(&data_.sp_handle);
  TEUCHOS_TEST_FOR_EXCEPTION( status != CUSOLVER_STATUS_SUCCESS,
    std::runtime_error, "cusolverSpCreate failed");
  status = cusolverSpCreateCsrcholInfo(&data_.chol_info);
  TEUCHOS_TEST_FOR_EXCEPTION( status != CUSOLVER_STATUS_SUCCESS,
    std::runtime_error, "cusolverSpCreateCsrcholInfo failed");
  auto sparse_status = cusparseCreateMatDescr(&data_.desc);
  TEUCHOS_TEST_FOR_EXCEPTION( sparse_status != CUSPARSE_STATUS_SUCCESS,
    std::runtime_error, "cusparseCreateMatDescr failed");
  auto blas_status = cublasCreate(&data_.blas_handle);
  TEUCHOS_TEST_FOR_EXCEPTION( blas_status != CUBLAS_STATUS_SUCCESS,
    std::runtime_error, "cublasCreate failed");
#ifdef KOKKOS_ENABLE_CUDA
  // Run GEMM on the same stream as Kokkos so copies and multiply are
  // serialized without cross-stream synchronization gaps.
  cublasSetStream(data_.blas_handle,
    Kokkos::DefaultExecutionSpace().cuda_stream());
#endif
}

template <class Matrix, class Vector>
cuSOLVER<Matrix,Vector>::~cuSOLVER( )
{
  cublasDestroy(data_.blas_handle);
  cusparseDestroyMatDescr(data_.desc);
  cusolverSpDestroyCsrcholInfo(data_.chol_info);
  cusolverSpDestroy(data_.sp_handle);
  cusolverDnDestroy(data_.dn_handle);
}

template<class Matrix, class Vector>
int
cuSOLVER<Matrix,Vector>::preOrdering_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif
  if(do_optimization() &&
     this->globalNumRows_ > data_.small_matrix_threshold) {
    this->matrixA_->returnRowPtr_kokkos_view(device_row_ptr_view_);
    this->matrixA_->returnColInd_kokkos_view(device_cols_view_);

    // reorder to optimize cuSolver
    if(data_.bReorder) {
      Amesos2::Util::reorder(
        device_row_ptr_view_, device_cols_view_,
        device_perm_, device_peri_, sorted_nnz,
        true);
    }
  }

  return 0;
}

template <class Matrix, class Vector>
int
cuSOLVER<Matrix,Vector>::symbolicFactorization_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor symFactTimer(this->timers_.symFactTime_);
#endif

  int err = 0;
  if ( this->root_ ) {
    const int n = this->globalNumRows_;

    if(n > data_.small_matrix_threshold) {
      if(device_row_ptr_view_.extent(0) == 0) {
        this->matrixA_->returnRowPtr_kokkos_view(device_row_ptr_view_);
        this->matrixA_->returnColInd_kokkos_view(device_cols_view_);
      }
      const int nnz = device_cols_view_.size(); // reorder may have changed this
      const int * colIdx = device_cols_view_.data();
      const int * rowPtr = device_row_ptr_view_.data();
      auto status = cusolverSpXcsrcholAnalysis(
        data_.sp_handle, n, nnz, data_.desc, rowPtr, colIdx, data_.chol_info);
      err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
    } else {
      // Allocate dense matrix, pivot array, info scalar, and inverse.
      if((int)device_matrix_.extent(0) != n) {
        device_matrix_ = device_value_type_matrix(
          Kokkos::ViewAllocateWithoutInitializing("cusolver_dense"), n, n);
        device_ipiv_ = Kokkos::View<int*, device_type>(
          Kokkos::ViewAllocateWithoutInitializing("cusolver_ipiv"), n);
        device_info_ = Kokkos::View<int, device_type>("cusolver_info");
      }
      if((int)device_inverse_.extent(0) != n) {
        device_inverse_ = device_value_type_matrix(
          Kokkos::ViewAllocateWithoutInitializing("cusolver_inverse"), n, n);
      }

      // Query factorization workspace size
      int lwork = 0;
      auto status = function_map::bufferInfo(
        data_.dn_handle, n, device_matrix_.data(), n, &lwork);
      if(status == CUSOLVER_STATUS_SUCCESS) {
        if((size_t)lwork > buffer_.extent(0)) {
          buffer_ = device_value_type_array(
            Kokkos::ViewAllocateWithoutInitializing("cusolver_buf"), lwork);
        }
      }
      err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
    }
  }

  Teuchos::broadcast(*(this->getComm()), 0, &err);
  TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
    std::runtime_error, "Amesos2 cuSolver symbolic failed.");

  return err;
}

template <class Matrix, class Vector>
int
cuSOLVER<Matrix,Vector>::numericFactorization_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

  int err = 0;
  if(do_optimization()) {
    const int n = this->globalNumRows_;

    if(n > data_.small_matrix_threshold) {
      this->matrixA_->returnValues_kokkos_view(device_nzvals_view_);

      // reorder to optimize cuSolver
      if(data_.bReorder) {
        // must have original row and cols - maybe cache this from 1st symbolic setup
        // this setup exists to support the refactor option
        device_size_type_array orig_device_row_ptr_view;
        device_ordinal_type_array orig_device_cols_view;
        this->matrixA_->returnRowPtr_kokkos_view(orig_device_row_ptr_view);
        this->matrixA_->returnColInd_kokkos_view(orig_device_cols_view);
        Amesos2::Util::reorder_values(
          device_nzvals_view_, orig_device_row_ptr_view, device_row_ptr_view_,
          orig_device_cols_view, device_perm_, device_peri_, sorted_nnz);
      }

      const int nnz = device_cols_view_.size(); // reorder may have changed this
      const cusolver_type * values = device_nzvals_view_.data();
      const int * colIdx = device_cols_view_.data();
      const int * rowPtr = device_row_ptr_view_.data();

      size_t internalDataInBytes, workspaceInBytes;
      auto status = function_map::sparseBufferInfo(
        data_.sp_handle, n, nnz, data_.desc,
        values, rowPtr, colIdx, data_.chol_info,
        &internalDataInBytes, &workspaceInBytes);

      if(status == CUSOLVER_STATUS_SUCCESS) {
        const size_t buffer_size = workspaceInBytes / sizeof(cusolver_type);
        if(buffer_size > buffer_.extent(0)) {
          buffer_ = device_value_type_array(
            Kokkos::ViewAllocateWithoutInitializing("cusolver_buf"), buffer_size);
        }
        status = function_map::sparseNumeric(
          data_.sp_handle, n, nnz, data_.desc,
          values, rowPtr, colIdx, data_.chol_info, buffer_.data());
      }
      err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
      Teuchos::broadcast(*(this->getComm()), 0, &err);
      TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
        std::runtime_error, "Amesos2 cuSolver numeric failed.");

      return err;
    }

    // Extract CSR views from matrix (views into matrix internal storage)
    device_value_type_array  nzvals;
    device_size_type_array   row_ptr;
    device_ordinal_type_array col_ind;
    this->matrixA_->returnValues_kokkos_view(nzvals);
    this->matrixA_->returnRowPtr_kokkos_view(row_ptr);
    this->matrixA_->returnColInd_kokkos_view(col_ind);

    // Zero-fill dense matrix then scatter sparse entries
    Kokkos::deep_copy(device_matrix_,
      Teuchos::ScalarTraits<cusolver_type>::zero());
    Impl::CsrToDenseFunctor<device_value_type_matrix,
                            device_size_type_array,
                            device_ordinal_type_array,
                            device_value_type_array>
      scatter(device_matrix_, row_ptr, col_ind, nzvals);
    Kokkos::parallel_for("Amesos2_cuSOLVER_csr_to_dense",
      Kokkos::RangePolicy<typename device_type::execution_space>(0, n),
      scatter);

    // LU factorization in-place
    auto status = function_map::numeric(
      data_.dn_handle, n, device_matrix_.data(), n,
      buffer_.data(), device_ipiv_.data(), device_info_.data());

    if(status == CUSOLVER_STATUS_SUCCESS) {
      // Fill device_inverse_ with identity, then solve LU * inv = I, inv = A^{-1}.
      Impl::SetIdentityFunctor<device_value_type_matrix> set_id(device_inverse_);
      Kokkos::parallel_for("Amesos2_cuSOLVER_set_identity",
        Kokkos::RangePolicy<typename device_type::execution_space>(0, n),
        set_id);

      status = function_map::solveLU(
        data_.dn_handle, CUBLAS_OP_N, n, n,
        device_matrix_.data(), n, device_ipiv_.data(),
        device_inverse_.data(), n, device_info_.data());
    }

    if(status == CUSOLVER_STATUS_SUCCESS) {
      auto host_info = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), device_info_);
      err = (host_info() != 0) ? 1 : 0;
    } else {
      err = 1;
    }
  }

  Teuchos::broadcast(*(this->getComm()), 0, &err);
  TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
    std::runtime_error, "Amesos2 cuSolver numeric failed.");

  return err;
}

template <class Matrix, class Vector>
int
cuSOLVER<Matrix,Vector>::solve_impl(
  const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
  const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  int err = 0;
  const int n = this->globalNumRows_;

  if(n > data_.small_matrix_threshold) {
    TEUCHOS_TEST_FOR_EXCEPTION(this->control_.useTranspose_,
      std::runtime_error,
      "Amesos2 cuSolver sparse Cholesky path does not support transpose solves.");

    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    const ordinal_type nrhs = X->getGlobalNumVectors();

    bool bAssignedX;
    {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
#endif

      const bool initialize_data = true;
      const bool do_not_initialize_data = false;
      Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        device_solve_array_t>::do_get(initialize_data, B, this->bValues_,
        Teuchos::as<size_t>(ld_rhs), ROOTED, this->rowIndexBase_);

      bAssignedX = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        device_solve_array_t>::do_get(do_not_initialize_data, X, this->xValues_,
        Teuchos::as<size_t>(ld_rhs), ROOTED, this->rowIndexBase_);
    }

    if ( this->root_ ) {  // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

      if(data_.bReorder) {
        Amesos2::Util::apply_reorder_permutation(
          this->bValues_, this->permute_result_, this->device_perm_);
      }
      else {
        this->permute_result_ = this->bValues_; // no permutation
      }

      for(ordinal_type rhs = 0; rhs < nrhs; ++rhs) {
        const cusolver_type * b = this->permute_result_.data() + rhs * n;
        cusolver_type * x = this->xValues_.data() + rhs * n;
        auto status = function_map::sparseSolve(
          data_.sp_handle, n, b, x, data_.chol_info, buffer_.data());
        err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
        if(err != 0) {
          break;
        }
      }

      if(data_.bReorder && err == 0) {
        Amesos2::Util::apply_reorder_permutation(
          this->xValues_, this->permute_result_, this->device_peri_);
        Kokkos::deep_copy(this->xValues_, this->permute_result_);
      }
    }

    if(!bAssignedX) {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      Util::template put_1d_data_helper_kokkos_view<
        MultiVecAdapter<Vector>,device_solve_array_t>::do_put(X, xValues_,
        Teuchos::as<size_t>(ld_rhs), ROOTED, this->rowIndexBase_);
    }

    Teuchos::broadcast(*(this->getComm()), 0, &err);
    TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
      std::runtime_error, "Amesos2 cuSolver solve failed.");

    return err;
  }

  if(this->root_) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
    const int nrhs = static_cast<int>(X->getGlobalNumVectors());

    // Zero-copy: read B and write X directly through their device views.
    // Both B and X are device-resident (do_optimization() guarantees single-process).
    auto b_view = B->getLocalDeviceView2d_ReadOnly();
    auto x_view = X->getLocalDeviceView2d_OverwriteAll();

    // stride(1) = distance between consecutive columns (includes any MV padding).
    const int ldb = std::max(n, static_cast<int>(b_view.stride(1)));
    const int ldx = std::max(n, static_cast<int>(x_view.stride(1)));
    const cublasOperation_t trans =
      this->control_.useTranspose_ ? CUBLAS_OP_C : CUBLAS_OP_N;

    auto blas_status = function_map::solve(
      data_.blas_handle, trans, n, nrhs,
      device_inverse_.data(), n,
      reinterpret_cast<const cusolver_type*>(b_view.data()), ldb,
      reinterpret_cast<cusolver_type*>(x_view.data()), ldx);

    err = (blas_status != CUBLAS_STATUS_SUCCESS) ? 1 : 0;
  }

  Teuchos::broadcast(*(this->getComm()), 0, &err);
  TEUCHOS_TEST_FOR_EXCEPTION(err != 0,
    std::runtime_error, "Amesos2 cuSolver solve failed.");

  return err;
}

template <class Matrix, class Vector>
bool
cuSOLVER<Matrix,Vector>::matrixShapeOK_impl() const
{
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}

template <class Matrix, class Vector>
void
cuSOLVER<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  if( parameterList->isParameter("Reorder") ){
    RCP<const ParameterEntryValidator> reorder_validator = valid_params->getEntry("Reorder").validator();
    parameterList->getEntry("Reorder").setValidator(reorder_validator);
  }
  if( parameterList->isParameter("small matrix threshold") ){
    RCP<const ParameterEntryValidator> threshold_validator =
      valid_params->getEntry("small matrix threshold").validator();
    parameterList->getEntry("small matrix threshold").setValidator(threshold_validator);
  }
  data_.bReorder = parameterList->get<bool>("Reorder",
#ifdef HAVE_AMESOS2_METIS
    true
#else
    false
#endif
  );
  data_.small_matrix_threshold = parameterList->get<int>(
    "small matrix threshold", 2500);
}

template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
cuSOLVER<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if (is_null(valid_params)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

#ifdef HAVE_AMESOS2_METIS
    pl->set("Reorder", true, "Whether GIDs contiguous");
#else
    pl->set("Reorder", false, "Whether GIDs contiguous");
#endif
    pl->set("small matrix threshold", 2500,
      "Use explicit inverse for matrices with rows less than or equal to this value");
    valid_params = pl;
  }
  return valid_params;
}

template <class Matrix, class Vector>
bool
cuSOLVER<Matrix,Vector>::do_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1));
}

template <class Matrix, class Vector>
bool
cuSOLVER<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  if(current_phase == SOLVE) {
    return(false);
  }

  if(!do_optimization()) {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,
      "cuSolver is only implemented for serial.");
  }

  return true;
}

template <class Matrix, class Vector>
void
cuSOLVER<Matrix,Vector>::describe_impl(Teuchos::FancyOStream &out,
                                       const Teuchos::EVerbosityLevel verbLevel) const
{
  out << " cuSOLVER: sparse Cholesky with dense inverse for small matrices" << std::endl;
  out << "  > Reorder = " << (data_.bReorder ? "YES" : "NO") << std::endl;
  out << "  > small matrix threshold = "
      << data_.small_matrix_threshold << std::endl;
}

template<class Matrix, class Vector>
const char* cuSOLVER<Matrix,Vector>::name = "cuSOLVER";

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_DEF_HPP
