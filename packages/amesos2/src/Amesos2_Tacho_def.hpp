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
  data_.Symbolic = NULL;
  data_.Numeric = NULL;
}


template <class Matrix, class Vector>
TachoSolver<Matrix,Vector>::~TachoSolver( )
{
  delete data_.Symbolic;
  delete data_.Numeric;
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
  if ( this->root_ ) {
    if(data_.Symbolic) {
      delete data_.Symbolic;
      data_.Symbolic = NULL;
    }

    size_type_array row_ptr;
    ordinal_type_array cols;

#ifndef HAVE_TEUCHOS_COMPLEX
    if(single_process_optim_check()) {
      // in the optimized case we read the values directly from the matrix
      // without converting through the Teuchos::Array setup. Note that in
      // this case nzvals_, colind_, and rowptr_ are never set in loadA_impl.
      // For HAVE_TEUCHOS_COMPLEX is_optimized_case() return false because we
      // need nzvals_.
      auto sp_rowptr = this->matrixA_->returnRowPtr();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null");
      auto sp_colind = this->matrixA_->returnColInd();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null");

      // Tacho uses size_type size_t which matches Tpetra but not Epetra (int)
      // So we need a converter
      row_ptr = size_type_array("r", this->globalNumRows_ + 1);
      for(global_size_type n = 0; n < this->globalNumRows_ + 1; ++n) {
        row_ptr(n) = static_cast<size_type>(sp_rowptr[n]);
      }

      // TODO - For Tpetra, we could have a direct view like this...
      // row_ptr = size_type_array(sp_rowptr, this->globalNumRows_ + 1);

      // cols type is ok for both Tpetra and Epetra already
      cols = ordinal_type_array(sp_colind, this->globalNumNonZeros_);
    }
    else
#endif
    {
      // Non optimized case used the arrays set up in loadA_impl
      row_ptr = size_type_array(this->rowptr_.getRawPtr(), this->globalNumRows_ + 1);
      cols = ordinal_type_array(this->colind_.getRawPtr(), this->globalNumNonZeros_);
    }

    data_.Symbolic = new Tacho::Experimental::SymbolicTools(
      this->globalNumCols_, row_ptr, cols, idx_, idx_);
    data_.Symbolic->symbolicFactorize();
  }

  return status;
}


template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::numericFactorization_impl()
{
  int status = 0;
  if ( this->root_ ) {
    if(!data_.Symbolic) {
      symbolicFactorization_impl();
    }
    if(data_.Numeric) {
      delete data_.Numeric;
      data_.Numeric = NULL;
    }
    size_type_array row_ptr;
    ordinal_type_array cols;
    value_type_array values;

#ifndef HAVE_TEUCHOS_COMPLEX
    if(single_process_optim_check()) {
      // in the optimized case we read the values directly from the matrix
      // without converting through the Teuchos::Array setup. Note that in
      // this case nzvals_, colind_, and rowptr_ are never set in loadA_impl.
      // For HAVE_TEUCHOS_COMPLEX is_optimized_case() return false because we
      // need nzvals_.
      auto sp_rowptr = this->matrixA_->returnRowPtr();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
         std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null");
      auto sp_colind = this->matrixA_->returnColInd();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null");
      auto sp_values = this->matrixA_->returnValues();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_values returned null");

      // Tacho uses size_type size_t which matches Tpetra but not Epetra (int)
      // So we need a converter
      row_ptr = size_type_array("r", this->globalNumRows_ + 1);
      for(global_size_type n = 0; n < this->globalNumRows_ + 1; ++n) {
        row_ptr(n) = static_cast<size_type>(sp_rowptr[n]);
      }

      // TODO - For Tpetra, we could have a direct view like this...
      // row_ptr = size_type_array(sp_rowptr, this->globalNumRows_ + 1);

      // cols and values are ok for both Tpetra and Epetra already
      cols = ordinal_type_array(sp_colind, this->globalNumNonZeros_);
      values = value_type_array(sp_values, this->globalNumNonZeros_);
    }
    else
#endif
    {
      // Non optimized case used the arrays set up in loadA_impl
      row_ptr = size_type_array(this->rowptr_.getRawPtr(), this->globalNumRows_ + 1);
      cols = ordinal_type_array(this->colind_.getRawPtr(), this->globalNumNonZeros_);
      values = value_type_array(this->nzvals_.getRawPtr(), this->globalNumNonZeros_);
    }

    data_.Numeric = new Tacho::Experimental::NumericTools<scalar_type,DeviceSpaceType>(
      this->globalNumRows_,
      row_ptr,
      cols,
      idx_,
      idx_,
      data_.Symbolic->NumSupernodes(),
      data_.Symbolic->Supernodes(),
      data_.Symbolic->gidSuperPanelPtr(),
      data_.Symbolic->gidSuperPanelColIdx(),
      data_.Symbolic->sidSuperPanelPtr(),
      data_.Symbolic->sidSuperPanelColIdx(),
      data_.Symbolic->blkSuperPanelColIdx(),
      data_.Symbolic->SupernodesTreeParent(),
      data_.Symbolic->SupernodesTreePtr(),
      data_.Symbolic->SupernodesTreeChildren(),
      data_.Symbolic->SupernodesTreeRoots());

    data_.Numeric->factorizeCholesky_Serial(values);
  }

  return status;
}

template <class Matrix, class Vector>
int
TachoSolver<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  const size_t val_store_size = as<size_t>(ld_rhs * nrhs);
  Teuchos::Array<tacho_type> xValues(val_store_size);
  Teuchos::Array<tacho_type> bValues(val_store_size);

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif
    Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
                             tacho_type>::do_get(B, bValues(),
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
  }

  // TODO: Decide if Transpose should be tested?
  // Currently this is turned on in the test but does nothing - to fix or delete
  // int TachoRequest = this->control_.useTranspose_ ?
  //   Tacho::Experimental::Trans::Transpose :
  //   Tacho::Experimental::Trans::NoTranspose;

  int ierr = 0; // returned error code

  if ( this->root_ ) {
    {                           // Do solve!
  #ifdef HAVE_AMESOS2_TIMER
      Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
  #endif
      if (data_.Symbolic) {
        delete data_.Symbolic;
        data_.Symbolic = NULL;
      }

      // validate
      int i_ld_rhs = as<int>(ld_rhs);

      for(size_t j = 0 ; j < nrhs; j++) {
        typedef typename Tacho::Experimental::NumericTools
          <tacho_type,DeviceSpaceType>::value_type_matrix_host solve_array_t;

        // TODO Decide how to handle t for Tacho and verify this setup is ok
        const int n = 1;
        solve_array_t x(&xValues.getRawPtr()[j*i_ld_rhs], this->globalNumRows_, n);
        solve_array_t b(&bValues.getRawPtr()[j*i_ld_rhs], this->globalNumRows_, n);
        solve_array_t t("t", this->globalNumRows_, n);

        // Still need to implement parallel, decided about thread handling, etc.
        data_.Numeric->solveCholesky_Serial(x, b, t);

        int status = 0; // TODO - determine what error handling will be

        if(status != 0) {
          ierr = status;
          break;
        }
      }
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

    Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,tacho_type>::do_put(X, xValues(),
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
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
TachoSolver<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    valid_params = pl;
  }

  return valid_params;
}

template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::single_process_optim_check() const {
  return ( (this->root_) &&
           (this->matrixA_->getComm()->getRank() == 0) &&
           (this->matrixA_->getComm()->getSize() == 1) );
}

template <class Matrix, class Vector>
bool
TachoSolver<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  if(current_phase == SOLVE) {
    return(false);
  }

  // TODO: How should we handle generating/storing this perm index?
  if( this->root_ ) {
    idx_ = ordinal_type_array("idx", this->globalNumNonZeros_);
    for (size_t i=0; i<this->idx_.size(); ++i) {
      idx_(i) = i;
    }
  }

#ifndef HAVE_TEUCHOS_COMPLEX
  if(single_process_optim_check()) {
    // Do nothing
  }
  else
#endif
  {
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // Only the root image needs storage allocated
    if( this->root_ ) {
      nzvals_.resize(this->globalNumNonZeros_);
      colind_.resize(this->globalNumNonZeros_);
      rowptr_.resize(this->globalNumRows_ + 1);
    }

    size_type nnz_ret = 0;
    {
  #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
  #endif

      TEUCHOS_TEST_FOR_EXCEPTION( this->rowIndexBase_ != this->columnIndexBase_,
                          std::runtime_error,
                          "Row and column maps have different indexbase ");

      Util::get_crs_helper<
      MatrixAdapter<Matrix>,tacho_type,ordinal_type,size_type>::do_get(
                                                      this->matrixA_.ptr(),
                                                      nzvals_(), colind_(),
                                                      rowptr_(), nnz_ret,
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
