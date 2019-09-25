// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Sivasankaran Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// This file is used for TachoSolver and TachoHostSolver
// So include guards exist twice
#if (defined(TACHO_BUILD_SOLVER) && !defined(TACHO_BUILT_SOLVER_DEF_HPP)) || \
    (defined(TACHOHOST_BUILD_SOLVER) && !defined(TACHOHOST_BUILT_SOLVER_DEF_HPP))

#ifdef TACHO_BUILD_SOLVER
  #define TACHO_BUILT_SOLVER_DEF_HPP
  #define TACHO_SOLVER_CHAR_NAME "Tacho"
#endif

#ifdef TACHOHOST_BUILD_SOLVER
  #define TACHOHOST_BUILT_SOLVER_DEF_HPP
  #define TACHO_SOLVER_CHAR_NAME "TachoHost"
#endif

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Tacho_decl.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

template <class Matrix, class Vector>
TACHO_SOLVER_NAME<Matrix,Vector>::TACHO_SOLVER_NAME(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::TACHO_SOLVER_NAME,Matrix,Vector>(A, X, B)
{
}


template <class Matrix, class Vector>
TACHO_SOLVER_NAME<Matrix,Vector>::~TACHO_SOLVER_NAME( )
{
}

template <class Matrix, class Vector>
std::string
TACHO_SOLVER_NAME<Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << "Tacho solver interface";
  return oss.str();
}

template<class Matrix, class Vector>
int
TACHO_SOLVER_NAME<Matrix,Vector>::preOrdering_impl()
{
  return(0);
}

template <class Matrix, class Vector>
int
TACHO_SOLVER_NAME<Matrix,Vector>::symbolicFactorization_impl()
{
  int status = 0;

  if ( this->root_ ) {
    if(do_optimization()) {
      // MDM-TODO some of these comments are about to become resolved so return
      // and remove as necessary.
      /*
          A couple of things to resolve/discuss here:

            (1) Do we want to try and do a direct view when the types match?
                If we do, then we need something that compiles for non-matching.
                I think that forces us to have a force cast - something like below.

            (2) Is this a problem for complex? Basker prevents optimization for
                complex but I wasn't sure if that necessary.

            (3) Can I safely read the size and index types directly without using
                std::remove_pointer to determine the regular type.

      */

      typedef typename Matrix::execution_space m_exec_space;

      // first create a View in the matrix mem space as efficiently as possible,
      // either directly, or by copy if types mismatch
      typedef typename MatrixAdapter<Matrix>::spmtx_ptr_t row_ptr_t;
      row_ptr_t sp_rowptr = this->matrixA_->returnRowPtr();
      TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null");
      Kokkos::View<size_type*, m_exec_space> row_ptr_view_src;
      if(std::is_same<size_type, typename std::remove_pointer<row_ptr_t>::type>::value) {
        // they match so get the view in matrix space and then copy or assign to us
        Kokkos::View<size_type*, m_exec_space> src_row_ptr((size_type*)sp_rowptr, this->globalNumRows_ + 1);
        deep_copy_or_assign_view(host_row_ptr_view_, src_row_ptr);
      }
      else {
        // they don't match so we'll need to loop - first make the view in matrix space
        // but with the dst type, copy the types, then send to the copy manager
        Kokkos::View<size_type*, m_exec_space> src_row_ptr(
          Kokkos::ViewAllocateWithoutInitializing("src_row_ptr"), this->globalNumRows_ + 1);

        // now copy in the matrix exec space
        Kokkos::parallel_for(Kokkos::RangePolicy<m_exec_space, global_size_type>
          (0, this->globalNumRows_ + 1), KOKKOS_LAMBDA (global_size_type n) {
          src_row_ptr(n) = static_cast<size_type>(sp_rowptr[n]);
        });
        deep_copy_or_assign_view(host_row_ptr_view_, src_row_ptr);
      }

      // do the same procedure as above, but for cols
      typedef typename MatrixAdapter<Matrix>::spmtx_idx_t col_ind_t;
      col_ind_t sp_colind = this->matrixA_->returnColInd();
      TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null");
      if(std::is_same<ordinal_type, typename std::remove_pointer<col_ind_t>::type>::value) {
        // they match so get the view in matrix space and then copy or assign to us
        Kokkos::View<ordinal_type*, m_exec_space> src_cols((ordinal_type*)sp_colind, this->globalNumNonZeros_);
        deep_copy_or_assign_view(host_cols_view_, src_cols);
      }
      else {
        // they don't match so we'll need to loop - first make the view in matrix space
        // but with the dst type, copy the types, then send to the copy manager
        Kokkos::View<ordinal_type*, m_exec_space> src_cols(
          Kokkos::ViewAllocateWithoutInitializing("src_cols"), this->globalNumNonZeros_);

        // now copy in the matrix exec space
        Kokkos::parallel_for(Kokkos::RangePolicy<m_exec_space, global_size_type>
          (0, this->globalNumNonZeros_), KOKKOS_LAMBDA (global_size_type n) {
          src_cols(n) = static_cast<ordinal_type>(sp_colind[n]);
        });
        deep_copy_or_assign_view(host_cols_view_, src_cols);
      }
    }
    else {
      // if both are host, then this will act like a host mirror and just point
      // to itself. But if the src was on Cuda (UVM off) then this will do a
      // a deep copy to prepare for the symbolic
      deep_copy_or_assign_view(host_row_ptr_view_, device_row_ptr_view_);
      deep_copy_or_assign_view(host_cols_view_, device_cols_view_);
    }

    // TODO: Confirm param options
    // data_.solver.setMaxNumberOfSuperblocks(data_.max_num_superblocks);

    // Symbolic factorization currently must be done on host
    data_.solver.analyze(this->globalNumCols_, host_row_ptr_view_, host_cols_view_);
  }

  return status;
}


template <class Matrix, class Vector>
int
TACHO_SOLVER_NAME<Matrix,Vector>::numericFactorization_impl()
{
  int status = 0;

  if ( this->root_ ) {
    device_value_type_array values;

    if(do_optimization()) {
      // in the optimized case we read the values directly from the matrix
      typename MatrixAdapter<Matrix>::spmtx_vals_t sp_values = this->matrixA_->returnValues();
      TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
        std::runtime_error, "Amesos2 Runtime Error: sp_values returned null");
      // this type should always be ok for float/double because Tacho solver
      // is templated to value_type to guarantee a match here.
      values = device_value_type_array(sp_values, this->globalNumNonZeros_);
    }
    else
    {
      // Non optimized case used the arrays set up in loadA_impl
      values = this->device_nzvals_view_;
    }

    data_.solver.factorize(values);
  }

  return status;
}

template <class Matrix, class Vector>
int
TACHO_SOLVER_NAME<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  // if we did match the source we don't allocate
  if(xValues_.extent(0) < ld_rhs || xValues_.extent(1) < nrhs) {
    // allocate x because the solve expects it
    xValues_ = device_solve_array_t(
      Kokkos::ViewAllocateWithoutInitializing("xValues"), ld_rhs, nrhs);
  }

  // don't allocate b since it's handled by the copy manager and might just be
  // be assigned, not copied anyways.

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
                             device_solve_array_t>::do_get(B, this->bValues_,
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
  }

  int ierr = 0; // returned error code

  if ( this->root_ ) {  // Do solve!
#ifdef HAVE_AMESOS2_TIMER
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
    // Bump up the workspace size if needed
    // MDM-TODO This was assuming extra size is ok but now with views it seems
    // a little dangerous - would we prefer to have this be == checks so the
    // workspace always matches exactly? Need to learn about internal handling.
    if (workspace_.extent(0) < this->globalNumRows_ || workspace_.extent(1) < nrhs) {
      workspace_ = device_solve_array_t(
        Kokkos::ViewAllocateWithoutInitializing("t"), this->globalNumRows_, nrhs);
    }

    data_.solver.solve(xValues_, bValues_, workspace_);

    int status = 0; // TODO: determine what error handling will be
    if(status != 0) {
      ierr = status;
    }
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  TEUCHOS_TEST_FOR_EXCEPTION( ierr != 0, std::runtime_error,
    "tacho_solve has error code: " << ierr );

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    Util::template put_1d_data_helper_kokkos_view<
      MultiVecAdapter<Vector>,device_solve_array_t>::do_put(X, xValues_,
                                        as<size_t>(ld_rhs),
                                        ROOTED, this->rowIndexBase_);
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
TACHO_SOLVER_NAME<Matrix,Vector>::matrixShapeOK_impl() const
{
  // Tacho can only apply the solve routines to square matrices
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
TACHO_SOLVER_NAME<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  // TODO: Confirm param options
  // data_.num_kokkos_threads = parameterList->get<int>("kokkos-threads", 1);
  // data_.max_num_superblocks = parameterList->get<int>("max-num-superblocks", 4);
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
TACHO_SOLVER_NAME<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    // TODO: Confirm param options
    // pl->set("kokkos-threads", 1, "Number of threads");
    // pl->set("max-num-superblocks", 4, "Max number of superblocks");

    valid_params = pl;
  }

  return valid_params;
}

template <class Matrix, class Vector>
bool
TACHO_SOLVER_NAME<Matrix,Vector>::do_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1));
}

template <class Matrix, class Vector>
bool
TACHO_SOLVER_NAME<Matrix,Vector>::loadA_impl(EPhase current_phase)
{

  if(current_phase == SOLVE) {
    return(false);
  }

  if(!do_optimization()) {
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Note views are allocated but eventually we should remove this.
    // The internal copy manager will decide if we can assign or deep_copy
    // and then allocate if necessary. However I don't have working Tacho GPU
    // right now for multiple ranks. So this code won't use the copy manager
    // yet. Therefore allocate here to have it working at least for non-GPU
    // parallel builds with more than 1 rank/
    // Only the root image needs storage allocated
    if( this->root_ ) {
      device_nzvals_view_ = device_value_type_array(
        Kokkos::ViewAllocateWithoutInitializing("nzvals"), this->globalNumNonZeros_);
      device_cols_view_ = device_ordinal_type_array(
        Kokkos::ViewAllocateWithoutInitializing("colind"), this->globalNumNonZeros_);
      device_row_ptr_view_ = device_size_type_array(
        Kokkos::ViewAllocateWithoutInitializing("rowptr"), this->globalNumRows_ + 1);
    }

    size_type nnz_ret = 0;
    {
  #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
  #endif

      TEUCHOS_TEST_FOR_EXCEPTION( this->rowIndexBase_ != this->columnIndexBase_,
                          std::runtime_error,
                          "Row and column maps have different indexbase ");

      Util::get_crs_helper_kokkos_view<
      MatrixAdapter<Matrix>,device_value_type_array,device_ordinal_type_array,device_size_type_array>::do_get(
                                                      this->matrixA_.ptr(),
                                                      device_nzvals_view_,
                                                      device_cols_view_,
                                                      device_row_ptr_view_,
                                                      nnz_ret,
                                                      ROOTED, ARBITRARY,
                                                      this->columnIndexBase_);
    }
  }

  return true;
}

template<class Matrix, class Vector>
const char* TACHO_SOLVER_NAME<Matrix,Vector>::name = TACHO_SOLVER_CHAR_NAME;

} // end namespace Amesos2

#undef TACHO_SOLVER_CHAR_NAME

#endif  // AMESOS2_TACHO_DEF_HPP or AMESOS2_TACHOHOST_DEF_HPP
