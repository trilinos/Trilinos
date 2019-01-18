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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
   \file   Amesos2_Basker_def.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>

   \brief  Definitions for the Amesos2 Basker solver interface
*/


#ifndef AMESOS2_BASKER_DEF_HPP
#define AMESOS2_BASKER_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Basker_decl.hpp"

namespace Amesos2 {


template <class Matrix, class Vector>
Basker<Matrix,Vector>::Basker(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Basker,Matrix,Vector>(A, X, B)
  , nzvals_()                   // initialize to empty arrays
  , rowind_()
  , colptr_()
  , is_contiguous_(true)
//  , basker()
{
  //Nothing
}


template <class Matrix, class Vector>
Basker<Matrix,Vector>::~Basker( )
{}

template <class Matrix, class Vector>
bool
Basker<Matrix,Vector>::single_proc_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1) && is_contiguous_);
}

template<class Matrix, class Vector>
int
Basker<Matrix,Vector>::preOrdering_impl()
{
  /* TODO: Define what it means for Basker
   */
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif

  return(0);
}


template <class Matrix, class Vector>
int
Basker<Matrix,Vector>::symbolicFactorization_impl()
{
  /*No symbolic factoriztion*/
  return(0);
}


template <class Matrix, class Vector>
int
Basker<Matrix,Vector>::numericFactorization_impl()
{
  using Teuchos::as;

  int info = 0;
  if ( this->root_ ){
    { // Do factorization
  #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
  #endif

  #ifdef HAVE_AMESOS2_VERBOSE_DEBUG
      std::cout << "Basker:: Before numeric factorization" << std::endl;
      std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
      std::cout << "rowind_ : " << rowind_.toString() << std::endl;
      std::cout << "colptr_ : " << colptr_.toString() << std::endl;
  #endif
     
      info = basker.factor(this->globalNumRows_, this->globalNumCols_, this->globalNumNonZeros_, colptr_.getRawPtr(), rowind_.getRawPtr(), nzvals_.getRawPtr());

      // This is set after numeric factorization complete as pivoting can be used;
      // In this case, a discrepancy between symbolic and numeric nnz total can occur.
      this->setNnzLU( as<size_t>(basker.get_NnzLU() ) ) ;
    }

  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);

  //global_size_type info_st = as<global_size_type>(info);
  /* TODO : Proper error messages*/
  TEUCHOS_TEST_FOR_EXCEPTION( (info == -1) ,
    std::runtime_error,
    "Basker: Could not alloc space for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION( (info ==  -2),
    std::runtime_error,
    "Basker: Could not alloc needed work space");
  TEUCHOS_TEST_FOR_EXCEPTION( (info == -3) ,
    std::runtime_error,
    "Basker: Could not alloc additional memory needed for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION( (info > 0) ,
    std::runtime_error,
    "Basker: Zero pivot found at: " << info );

  return(info);
}


template <class Matrix, class Vector>
int
Basker<Matrix,Vector>::solve_impl(
 const Teuchos::Ptr<MultiVecAdapter<Vector> >  X,
 const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  int ierr = 0; // returned error code

  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  if ( single_proc_optimization() && nrhs == 1 ) {

#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

#ifndef HAVE_TEUCHOS_COMPLEX
    auto b_vector = Util::vector_pointer_helper< MultiVecAdapter<Vector>, Vector >::get_pointer_to_vector( B );
    auto x_vector = Util::vector_pointer_helper< MultiVecAdapter<Vector>, Vector >::get_pointer_to_vector( X );
#else
    // NDE: 09/25/2017
    // Cannot convert Kokkos::complex<T>* to std::complex<T>*; in this case, use reinterpret_cast
    using complex_type = typename Util::getStdCplxType< magnitude_type, typename matrix_adapter_type::spmtx_vals_t >::type;
    complex_type * b_vector = reinterpret_cast< complex_type * >( Util::vector_pointer_helper< MultiVecAdapter<Vector>, Vector >::get_pointer_to_vector( B ) );
    complex_type * x_vector = reinterpret_cast< complex_type * >( Util::vector_pointer_helper< MultiVecAdapter<Vector>, Vector >::get_pointer_to_vector( X ) );
#endif
    TEUCHOS_TEST_FOR_EXCEPTION(b_vector == nullptr,
        std::runtime_error, "Amesos2 Runtime Error: b_vector returned null ");

    TEUCHOS_TEST_FOR_EXCEPTION(x_vector  == nullptr,
        std::runtime_error, "Amesos2 Runtime Error: x_vector returned null ");

    if ( this->root_ ) {
      {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif
        ierr = basker.solveMultiple(nrhs, b_vector, x_vector);
      }

      /* All processes should have the same error code */
      Teuchos::broadcast(*(this->getComm()), 0, &ierr);

      TEUCHOS_TEST_FOR_EXCEPTION( ierr  > 0,
          std::runtime_error,
          "Encountered zero diag element at: " << ierr);
      TEUCHOS_TEST_FOR_EXCEPTION( ierr == -1,
          std::runtime_error,
          "Could not alloc needed working memory for solve" );
    }
  }
  else 
  {
    const size_t val_store_size = as<size_t>(ld_rhs * nrhs);

    xvals_.resize(val_store_size);
    bvals_.resize(val_store_size);

    {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      if ( is_contiguous_ == true ) {
        Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
          slu_type>::do_get(B, bvals_(), as<size_t>(ld_rhs), ROOTED, this->rowIndexBase_);
      }
      else {
        Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
          slu_type>::do_get(B, bvals_(), as<size_t>(ld_rhs), CONTIGUOUS_AND_ROOTED, this->rowIndexBase_);
      }
    }

    if ( this->root_ ) {
      {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

        ierr = basker.solveMultiple(nrhs, bvals_.getRawPtr(),xvals_.getRawPtr());
      }

    }

    /* All processes should have the same error code */
    Teuchos::broadcast(*(this->getComm()), 0, &ierr);

    TEUCHOS_TEST_FOR_EXCEPTION( ierr  > 0,
        std::runtime_error,
        "Encountered zero diag element at: " << ierr);
    TEUCHOS_TEST_FOR_EXCEPTION( ierr == -1,
        std::runtime_error,
        "Could not alloc needed working memory for solve" );

    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      if ( is_contiguous_ == true ) {
        Util::put_1d_data_helper<
          MultiVecAdapter<Vector>,slu_type>::do_put(X, xvals_(),
              as<size_t>(ld_rhs),
              ROOTED);
      }
      else {
        Util::put_1d_data_helper<
          MultiVecAdapter<Vector>,slu_type>::do_put(X, xvals_(),
              as<size_t>(ld_rhs),
              CONTIGUOUS_AND_ROOTED);
      }
    }
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
Basker<Matrix,Vector>::matrixShapeOK_impl() const
{
  // The Basker can only handle square for right now
  return( this->globalNumRows_ == this->globalNumCols_ );
}


template <class Matrix, class Vector>
void
Basker<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  if(parameterList->isParameter("IsContiguous"))
    {
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }

}

template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Basker<Matrix,Vector>::getValidParameters_impl() const
{
  using Teuchos::ParameterList;

  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set("IsContiguous", true, "Are GIDs contiguous");
    pl->set("alnnz",  2, "Approx number of nonzeros in L, default is 2*nnz(A)");
    pl->set("aunnx",  2, "Approx number of nonzeros in I, default is 2*nnz(U)");
    valid_params = pl;
  }
  return valid_params;
}


template <class Matrix, class Vector>
bool
Basker<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  using Teuchos::as;
  if(current_phase == SOLVE) return (false);

  #ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
  #endif

  // Only the root image needs storage allocated
  if( this->root_ ){
    nzvals_.resize(this->globalNumNonZeros_);
    rowind_.resize(this->globalNumNonZeros_);
    colptr_.resize(this->globalNumCols_ + 1);
  }

  local_ordinal_type nnz_ret = 0;
  {
  #ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
  #endif

    if ( is_contiguous_ == true ) {
      Util::get_ccs_helper<
        MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
        ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
            nnz_ret, ROOTED, ARBITRARY, this->rowIndexBase_);
    }
    else {
      Util::get_ccs_helper<
        MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
        ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
            nnz_ret, CONTIGUOUS_AND_ROOTED, ARBITRARY, this->rowIndexBase_);
    }
  }

  if( this->root_ ){
    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
                        std::runtime_error,
                        "Amesos2_Basker loadA_impl: Did not get the expected number of non-zero vals");
  }
  return true;
}


template<class Matrix, class Vector>
const char* Basker<Matrix,Vector>::name = "Basker";


} // end namespace Amesos2

#endif  // AMESOS2_Basker_DEF_HPP
