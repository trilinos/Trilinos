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

#ifndef AMESOS2_UMFPACK_DEF_HPP
#define AMESOS2_UMFPACK_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Umfpack_decl.hpp"
#include "Amesos2_Util.hpp"

namespace Amesos2 {

template <class Matrix, class Vector>
Umfpack<Matrix,Vector>::Umfpack(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Umfpack,Matrix,Vector>(A, X, B)
  , nzvals_()                   // initialize to empty arrays
  , rowind_()
  , colptr_()
  , is_contiguous_(true)
{
  data_.Symbolic = NULL;
  data_.Numeric = NULL;
}


template <class Matrix, class Vector>
Umfpack<Matrix,Vector>::~Umfpack( )
{
  if (data_.Symbolic) function_map::umfpack_free_symbolic (&data_.Symbolic);
  if (data_.Numeric) function_map::umfpack_free_numeric (&data_.Numeric);
}

template <class Matrix, class Vector>
std::string
Umfpack<Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << "Umfpack solver interface";
  return oss.str();
}

template<class Matrix, class Vector>
int
Umfpack<Matrix,Vector>::preOrdering_impl()
{
  return(0);
}


template <class Matrix, class Vector>
int
Umfpack<Matrix,Vector>::symbolicFactorization_impl()
{
  int status = 0;
  if ( this->root_ ) {
    if (data_.Symbolic) {
      function_map::umfpack_free_symbolic(&(data_.Symbolic));
    }

    function_map::umfpack_defaults(data_.Control);

    status = function_map::umfpack_symbolic(
      this->globalNumRows_,this->globalNumCols_,
      &(this->colptr_[0]), &(this->rowind_[0]),
      &(this->nzvals_[0]), &(data_.Symbolic), data_.Control, data_.Info);
  }

  return status;
}


template <class Matrix, class Vector>
int
Umfpack<Matrix,Vector>::numericFactorization_impl()
{
  int status = 0;
  if ( this->root_ ) {
    if(!data_.Symbolic) {
      symbolicFactorization_impl();
    }

    function_map::umfpack_defaults(data_.Control);

    if (data_.Numeric) {
      function_map::umfpack_free_numeric(&(data_.Numeric));
    }

    status = function_map::umfpack_numeric(
      &(this->colptr_[0]),
      &(this->rowind_[0]), &(this->nzvals_[0]), data_.Symbolic,
      &(data_.Numeric), data_.Control, data_.Info);
  }
  return status;
}

template <class Matrix, class Vector>
int
Umfpack<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  const size_t val_store_size = as<size_t>(ld_rhs * nrhs);
  Teuchos::Array<umfpack_type> xValues(val_store_size);
  Teuchos::Array<umfpack_type> bValues(val_store_size);

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif
    if ( is_contiguous_ == true ) {
      Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
                               umfpack_type>::do_get(B, bValues(),
                                                     as<size_t>(ld_rhs),
                                                     ROOTED, 
                                                     this->rowIndexBase_);
    } 
    else {
      Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
                               umfpack_type>::do_get(B, bValues(),
                                                     as<size_t>(ld_rhs),
                                                     CONTIGUOUS_AND_ROOTED, 
                                                     this->rowIndexBase_);
    }
  }

  int UmfpackRequest = this->control_.useTranspose_ ? UMFPACK_At : UMFPACK_A;

  int ierr = 0; // returned error code

  if ( this->root_ ) {
    {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMER
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
      if (data_.Symbolic) {
        function_map::umfpack_free_symbolic(&(data_.Symbolic));
      }

      // validate
      int i_ld_rhs = as<int>(ld_rhs);

      for(size_t j = 0 ; j < nrhs; j++) {
        int status = function_map::umfpack_solve(
          UmfpackRequest,
          &(this->colptr_[0]), &(this->rowind_[0]), &(this->nzvals_[0]),
          &xValues.getRawPtr()[j*i_ld_rhs],
          &bValues.getRawPtr()[j*i_ld_rhs],
          data_.Numeric, data_.Control, data_.Info);

        if(status != 0) {
          ierr = status;
          break; // need to verify best ways to handle error in this loop
        }
      }
    }
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  TEUCHOS_TEST_FOR_EXCEPTION( ierr != 0, std::runtime_error,
    "umfpack_solve has error code: " << ierr );

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    if ( is_contiguous_ == true ) {
      Util::put_1d_data_helper<
        MultiVecAdapter<Vector>,umfpack_type>::do_put(X, xValues(),
                                                      as<size_t>(ld_rhs),
                                                      ROOTED, 
                                                      this->rowIndexBase_);
    }
    else {
      Util::put_1d_data_helper<
        MultiVecAdapter<Vector>,umfpack_type>::do_put(X, xValues(),
                                                      as<size_t>(ld_rhs),
                                                      CONTIGUOUS_AND_ROOTED, 
                                                      this->rowIndexBase_);
    }
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
Umfpack<Matrix,Vector>::matrixShapeOK_impl() const
{
  // The Umfpack factorization routines can handle square as well as
  // rectangular matrices, but Umfpack can only apply the solve routines to
  // square matrices, so we check the matrix for squareness.
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
Umfpack<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  if( parameterList->isParameter("IsContiguous") ){
    is_contiguous_ = parameterList->get<bool>("IsContiguous");
  }
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Umfpack<Matrix,Vector>::getValidParameters_impl() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set("IsContiguous", true, "Whether GIDs contiguous");

    valid_params = pl;
  }

  return valid_params;
}


template <class Matrix, class Vector>
bool
Umfpack<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  if(current_phase == SOLVE) {
    return(false);
  }

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

  // Only the root image needs storage allocated
  if( this->root_ ){
    nzvals_.resize(this->globalNumNonZeros_);
    rowind_.resize(this->globalNumNonZeros_);
    colptr_.resize(this->globalNumCols_ + 1);
  }

  int nnz_ret = 0;
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    TEUCHOS_TEST_FOR_EXCEPTION( this->rowIndexBase_ != this->columnIndexBase_,
                        std::runtime_error,
                        "Row and column maps have different indexbase ");
    if ( is_contiguous_ == true ) {
      Util::get_ccs_helper<
      MatrixAdapter<Matrix>,umfpack_type,int,int>::do_get(this->matrixA_.ptr(),
                                                          nzvals_(), rowind_(),
                                                          colptr_(), nnz_ret, 
                                                          ROOTED,
                                                          ARBITRARY,
                                                          this->rowIndexBase_);
    } else {
      Util::get_ccs_helper<
      MatrixAdapter<Matrix>,umfpack_type,int,int>::do_get(this->matrixA_.ptr(),
                                                          nzvals_(), rowind_(),
                                                          colptr_(), nnz_ret, 
                                                          CONTIGUOUS_AND_ROOTED,
                                                          ARBITRARY,
                                                          this->rowIndexBase_);
    }
  }

  return true;
}


template<class Matrix, class Vector>
const char* Umfpack<Matrix,Vector>::name = "Umfpack";


} // end namespace Amesos2

#endif  // AMESOS2_UMFPACK_DEF_HPP
