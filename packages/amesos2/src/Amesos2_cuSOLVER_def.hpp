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

template <class Matrix, class Vector>
cuSOLVER<Matrix,Vector>::cuSOLVER(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::cuSOLVER,Matrix,Vector>(A, X, B)
{
  auto status = cusolverSpCreate(&data_.handle);
  TEUCHOS_TEST_FOR_EXCEPTION( status != CUSOLVER_STATUS_SUCCESS,
    std::runtime_error, "cusolverSpCreate failed");

  status = cusolverSpCreateCsrcholInfo(&data_.chol_info);
  TEUCHOS_TEST_FOR_EXCEPTION( status != CUSOLVER_STATUS_SUCCESS,
    std::runtime_error, "cusolverSpCreateCsrcholInfo failed");

  auto sparse_status = cusparseCreateMatDescr(&data_.desc);
  TEUCHOS_TEST_FOR_EXCEPTION( sparse_status != CUSPARSE_STATUS_SUCCESS,
    std::runtime_error, "cusparseCreateMatDescr failed");
}

template <class Matrix, class Vector>
cuSOLVER<Matrix,Vector>::~cuSOLVER( )
{
  cusparseDestroyMatDescr(data_.desc);
  cusolverSpDestroyCsrcholInfo(data_.chol_info);
  cusolverSpDestroy(data_.handle);
}

template<class Matrix, class Vector>
int
cuSOLVER<Matrix,Vector>::preOrdering_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif
  if(do_optimization()) {
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

  return(0);
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
    const int size = this->globalNumRows_;
    const int nnz = device_cols_view_.size(); // reorder may have changed this
    const int * colIdx = device_cols_view_.data();
    const int * rowPtr = device_row_ptr_view_.data();
    auto status = cusolverSpXcsrcholAnalysis(
      data_.handle, size, nnz, data_.desc, rowPtr, colIdx, data_.chol_info);
    err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
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
  int err = 0;
  if(do_optimization()) { // just supporting one rank right now
    this->matrixA_->returnValues_kokkos_view(device_nzvals_view_);

    // reorder to optimize cuSolver
    if(data_.bReorder) {
      // must have original row and cols - maybe cache this from 1st symbiolic setup
      // this setup exists to support the refactor option
      device_size_type_array orig_device_row_ptr_view;
      device_ordinal_type_array orig_device_cols_view;
      this->matrixA_->returnRowPtr_kokkos_view(orig_device_row_ptr_view);
      this->matrixA_->returnColInd_kokkos_view(orig_device_cols_view);
      Amesos2::Util::reorder_values(
        device_nzvals_view_, orig_device_row_ptr_view, device_row_ptr_view_, orig_device_cols_view,
        device_perm_, device_peri_, sorted_nnz);
    }

    const int size = this->globalNumRows_;
    const int nnz = device_cols_view_.size(); // reorder may have changed this
    const cusolver_type * values = device_nzvals_view_.data();
    const int * colIdx = device_cols_view_.data();
    const int * rowPtr = device_row_ptr_view_.data();

    size_t internalDataInBytes, workspaceInBytes;
    auto status = function_map::bufferInfo(data_.handle, size, nnz, data_.desc,
      values, rowPtr, colIdx, data_.chol_info,
      &internalDataInBytes, &workspaceInBytes);

    if(status == CUSOLVER_STATUS_SUCCESS) {
      const size_t buffer_size = workspaceInBytes / sizeof(cusolver_type);
      if(buffer_size > buffer_.extent(0)) {
        buffer_ = device_value_type_array(
          Kokkos::ViewAllocateWithoutInitializing("cusolver buf"), buffer_size);
      }
      status = function_map::numeric(data_.handle, size, nnz, data_.desc,
        values, rowPtr, colIdx, data_.chol_info, buffer_.data());
    }
    err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
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
  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const ordinal_type nrhs = X->getGlobalNumVectors();

  bool bAssignedX;
  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    const bool initialize_data = true;
    const bool do_not_initialize_data = false;
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      device_solve_array_t>::do_get(initialize_data, B, this->bValues_, Teuchos::as<size_t>(ld_rhs),
      ROOTED, this->rowIndexBase_);

    // In general we may want to write directly to the x space without a copy.
    // So we 'get' x which may be a direct view assignment to the MV.
    bAssignedX = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      device_solve_array_t>::do_get(do_not_initialize_data, X, this->xValues_, Teuchos::as<size_t>(ld_rhs),
      ROOTED, this->rowIndexBase_);
  }

  int err = 0;

  if ( this->root_ ) {  // Do solve!
#ifdef HAVE_AMESOS2_TIMER
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

    const int size = this->globalNumRows_;

    if(data_.bReorder) {
      Amesos2::Util::apply_reorder_permutation(
        this->bValues_, this->permute_result_, this->device_perm_);
    }
    else {
      this->permute_result_ = this->bValues_; // no permutation
    }

    for(ordinal_type n = 0; n < nrhs; ++n) {
      const cusolver_type * b = this->permute_result_.data() + n * size;
      cusolver_type * x = this->xValues_.data() + n * size;
      auto status = function_map::solve(
        data_.handle, size, b, x, data_.chol_info, buffer_.data());
      err = (status != CUSOLVER_STATUS_SUCCESS) ? 1 : 0;
      if(err != 0) {
        break;
      }
    }

    if(data_.bReorder && err == 0) {
      Amesos2::Util::apply_reorder_permutation(
        this->xValues_, this->permute_result_, this->device_peri_);
      Kokkos::deep_copy(this->xValues_, this->permute_result_); // full copy since permute_result_ is reused
    }
  }

  /* Update X's global values */

  // if bDidAssignX, then we solved straight to the adapter's X memory space without
  // requiring additional memory allocation, so the x data is already in place.
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

  data_.bReorder = parameterList->get<bool>("Reorder", true);
}

template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
cuSOLVER<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set("Reorder", true, "Whether GIDs contiguous");

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

  if(!do_optimization()) { // we're only doing serial right now for cuSolver
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,
      "cuSolver is only implemented for serial.");
  }

  return true;
}

template<class Matrix, class Vector>
const char* cuSOLVER<Matrix,Vector>::name = "cuSOLVER";

} // end namespace Amesos2

#endif  // AMESOS2_CUSOLVER_DEF_HPP
