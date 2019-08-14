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

#ifndef AMESOS2_TACHO_DEF_HPP
#define AMESOS2_TACHO_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Tacho_decl.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

template <class Matrix, class Vector>
TachoSolver<Matrix,Vector>::TachoSolver(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::TachoSolver,Matrix,Vector>(A, X, B)
  , nzvals_()                   // initialize to empty arrays
  , colind_()
  , rowptr_()
{
}


template <class Matrix, class Vector>
TachoSolver<Matrix,Vector>::~TachoSolver( )
{
}

template <class Matrix, class Vector>
std::string
TachoSolver<Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << "Tacho solver interface";
  return oss.str();
}

template<class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::preOrdering_impl()
{
  return(0);
}

template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::symbolicFactorization_impl()
{
  int status = 0;

  device_size_type_array row_ptr;
  device_ordinal_type_array cols;

  if ( this->root_ ) {
    if(do_optimization()) {
      // in the optimized case we read the values directly from the matrix
      typedef typename MatrixAdapter<Matrix>::spmtx_ptr_t row_ptr_t;
      typedef typename MatrixAdapter<Matrix>::spmtx_idx_t col_ind_t;
      row_ptr_t sp_rowptr = this->matrixA_->returnRowPtr();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null");
      col_ind_t sp_colind = this->matrixA_->returnColInd();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null");

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

      if(std::is_same<size_type, typename std::remove_pointer<row_ptr_t>::type>::value) {
        // This line is used when types match exactly, but we need compilation when they don't
        // match, hence the cast - not sure if there is a better way to do this.
        row_ptr = device_size_type_array((size_type*)sp_rowptr, this->globalNumRows_ + 1);
      }
      else {
        row_ptr = device_size_type_array("r", this->globalNumRows_ + 1);

        // Keep this on device? We could copy to host first, then loop.
        // Not sure if symbolic factorization might have a device mode later
        Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceSpaceType, global_size_type>
              (0, this->globalNumRows_ + 1), KOKKOS_LAMBDA (int n) {
          row_ptr(n) = static_cast<size_type>(sp_rowptr[n]);
        });
      }

      if(std::is_same<ordinal_type, typename std::remove_pointer<col_ind_t>::type>::value) {
        // cast is so that things compile when this line is not used - same issues as above.
        cols = device_ordinal_type_array((ordinal_type*)sp_colind, this->globalNumNonZeros_);
      }
      else {
        cols = device_ordinal_type_array("c", this->globalNumNonZeros_);

        // Keep this on device? We could copy to host first, then loop.
        // Not sure if symbolic factorization might have a device mode later
        Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceSpaceType, global_size_type>
              (0, this->globalNumNonZeros_), KOKKOS_LAMBDA (int n) {
          cols(n) = static_cast<ordinal_type>(sp_colind[n]);
        });
      }
    }
    else
    {
      // Non optimized case used the arrays set up in loadA_impl
      row_ptr = this->rowptr_;
      cols = this->colind_;
    }

    // TODO: Confirm param options
    // data_.solver.setMaxNumberOfSuperblocks(data_.max_num_superblocks);

    // Symbolic factorization currently must be done on host
    host_size_type_array host_row_ptr("host_row_ptr", this->globalNumRows_ + 1);
    host_ordinal_type_array host_cols("host_cols", this->globalNumNonZeros_);
    Kokkos::deep_copy(host_row_ptr, row_ptr);
    Kokkos::deep_copy(host_cols, cols);

    data_.solver.analyze(this->globalNumCols_, host_row_ptr, host_cols);
  }

  return status;
}


template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::numericFactorization_impl()
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
      values = this->nzvals_;
    }

    data_.solver.factorize(values);
  }

  return status;
}

template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  device_solve_array_t x("x", this->globalNumRows_, nrhs);
  device_solve_array_t b("b", this->globalNumRows_, nrhs);

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
                             device_solve_array_t>::do_get(B, b,
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
  }

  int ierr = 0; // returned error code

  if ( this->root_ ) {  // Do solve!
#ifdef HAVE_AMESOS2_TIMER
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
    // Bump up the workspace size if needed
    if (workspace_.extent(0) < this->globalNumRows_ || workspace_.extent(1) < nrhs) {
      workspace_ = device_solve_array_t("t", this->globalNumRows_, nrhs);
    }

    data_.solver.solve(x, b, workspace_);

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

    Util::put_1d_data_helper_kokkos_view<
      MultiVecAdapter<Vector>,device_solve_array_t>::do_put(X, x,
                                         as<size_t>(ld_rhs),
                                         ROOTED, this->rowIndexBase_);
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::matrixShapeOK_impl() const
{
  // Tacho can only apply the solve routines to square matrices
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
TachoSolver<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  // TODO: Confirm param options
  // data_.num_kokkos_threads = parameterList->get<int>("kokkos-threads", 1);
  // data_.max_num_superblocks = parameterList->get<int>("max-num-superblocks", 4);
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
TachoSolver<Matrix,Vector>::getValidParameters_impl() const
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
TachoSolver<Matrix,Vector>::do_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1));
}

template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  if(current_phase == SOLVE) {
    return(false);
  }

  if(!do_optimization()) {
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Only the root image needs storage allocated
    if( this->root_ ) {
      nzvals_ = device_value_type_array("nzvals", this->globalNumNonZeros_);
      colind_ = device_ordinal_type_array("colind", this->globalNumNonZeros_);
      rowptr_ = device_size_type_array("rowptr", this->globalNumRows_ + 1);
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
                                                      nzvals_,
                                                      colind_,
                                                      rowptr_,
                                                      nnz_ret,
                                                      ROOTED, ARBITRARY,
                                                      this->columnIndexBase_);
    }
  }

  return true;
}


template<class Matrix, class Vector>
const char* TachoSolver<Matrix,Vector>::name = "Tacho";


} // end namespace Amesos2

#endif  // AMESOS2_TACHO_DEF_HPP
