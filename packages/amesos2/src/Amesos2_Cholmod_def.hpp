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
   \file   Amesos2_Cholmod_def.hpp
   \author Kevin Deweese <kdeweese@sandia.gov>
   \date   Wed Jul  24 15::48:51 2013

   \brief  Definitions for the Amesos2 Cholmod solver interface
*/


#ifndef AMESOS2_CHOLMOD_DEF_HPP
#define AMESOS2_CHOLMOD_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Cholmod_decl.hpp"


namespace Amesos2 {


template <class Matrix, class Vector>
Cholmod<Matrix,Vector>::Cholmod(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Cholmod,Matrix,Vector>(A, X, B)
  , firstsolve(true)
  , skip_symfact(false)
  , map()
  , is_contiguous_(true)
{
  CHOL::cholmod_l_start(&data_.c); // long form required for CUDA

  // Always use GPU is possible - note this is called after cholmod_l_start.
  // If we don't call this, it remains the default value of -1.
  // Then GPU is only used if env variable CHOLMOD_USE_GPU is set to 1.
  data_.c.useGPU = 1;
  data_.c.nmethods = 9; // try all methods - did not explore this yet
  data_.c.quick_return_if_not_posdef = 1;
  data_.c.itype = CHOLMOD_LONG;

  data_.L = NULL;
  data_.Y = NULL;
  data_.E = NULL;
}


template <class Matrix, class Vector>
Cholmod<Matrix,Vector>::~Cholmod( )
{
  CHOL::cholmod_l_gpu_stats(&(data_.c));

  CHOL::cholmod_l_free_factor(&(data_.L), &(data_.c));
  CHOL::cholmod_l_free_dense(&(data_.Y), &data_.c);
  CHOL::cholmod_l_free_dense(&(data_.E), &data_.c);

  CHOL::cholmod_l_finish(&(data_.c));
}

template<class Matrix, class Vector>
int
Cholmod<Matrix,Vector>::preOrdering_impl()
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif

  int info = 0;

  data_.L = CHOL::cholmod_l_analyze(&data_.A, &(data_.c));
  info = data_.c.status;
  skip_symfact = true;

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);
  
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
    std::runtime_error,
    "Amesos2 cholmod_l_analyze failure in Cholmod preOrdering_impl");

  return(0);
}


template <class Matrix, class Vector>
int
Cholmod<Matrix,Vector>::symbolicFactorization_impl()
{
  int info = 0;

  if(!skip_symfact) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor symFactTimer(this->timers_.symFactTime_);
#endif
    CHOL::cholmod_l_resymbol (&data_.A, NULL, 0, true, data_.L, &(data_.c));
    info = data_.c.status;
  } else {
    /*
     * Symbolic factorization has already occured in preOrdering_impl,
     * but if the user calls this routine directly later, we need to
     * redo the symbolic factorization.
     */
    skip_symfact = false;
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);
  
  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
    std::runtime_error,
    "Amesos2 cholmod_l_resymbol failure in Cholmod symbolicFactorization_impl");

  return(0);
}


template <class Matrix, class Vector>
int
Cholmod<Matrix,Vector>::numericFactorization_impl()
{
  int info = 0;

#ifdef HAVE_AMESOS2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(data_.A.ncol != Teuchos::as<size_t>(this->globalNumCols_),
    std::runtime_error,
    "Error in converting to cholmod_sparse: wrong number of global columns." );
  TEUCHOS_TEST_FOR_EXCEPTION(data_.A.nrow != Teuchos::as<size_t>(this->globalNumRows_),
    std::runtime_error,
    "Error in converting to cholmod_sparse: wrong number of global rows." );
#endif

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
  // Add print out of views here - need host conversion first
#endif

  CHOL::cholmod_l_factorize(&data_.A, data_.L, &(data_.c));
  info = data_.c.status;

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);

  TEUCHOS_TEST_FOR_EXCEPTION(info == 2,
    std::runtime_error,
    "Memory allocation failure in Cholmod factorization");

  TEUCHOS_TEST_FOR_EXCEPTION(info == 1,
    std::runtime_error,
    "Amesos2 cholmod_l_factorize is attempting to use Cholmod features that are not available yet.");

  TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
    std::runtime_error,
    "Amesos2 cholmod_l_factorize failure in Cholmod factorization");

  return(info);
}


template <class Matrix, class Vector>
int
Cholmod<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  const global_size_type ld_rhs = X->getGlobalLength();
  const size_t nrhs = X->getGlobalNumVectors();

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(B, bValues_,
          Teuchos::as<size_t>(ld_rhs),
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          this->rowIndexBase_);

    // In general we may want to write directly to the x space without a copy.
    // So we 'get' x which may be a direct view assignment to the MV.
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(X, xValues_,
          Teuchos::as<size_t>(ld_rhs),
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          this->rowIndexBase_);
  }

  int ierr = 0; // returned error code

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
#endif

  function_map::cholmod_init_dense(Teuchos::as<long>(this->globalNumRows_),
    Teuchos::as<int>(nrhs), Teuchos::as<int>(ld_rhs), bValues_.data(), &data_.b);
  function_map::cholmod_init_dense(Teuchos::as<long>(this->globalNumRows_),
    Teuchos::as<int>(nrhs), Teuchos::as<int>(ld_rhs), xValues_.data(), &data_.x);

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

  CHOL::cholmod_dense *xtemp = &(data_.x);
  CHOL::cholmod_l_solve2(CHOLMOD_A, data_.L, &data_.b, NULL,
    &xtemp, NULL, &data_.Y, &data_.E, &data_.c);

  ierr = data_.c.status;

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  TEUCHOS_TEST_FOR_EXCEPTION(ierr == -2, std::runtime_error, "Ran out of memory" );

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    Util::put_1d_data_helper_kokkos_view<
      MultiVecAdapter<Vector>,host_solve_array_t>::do_put(X, xValues_,
          Teuchos::as<size_t>(ld_rhs),
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          this->rowIndexBase_);
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
Cholmod<Matrix,Vector>::matrixShapeOK_impl() const
{
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
Cholmod<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  if( parameterList->isParameter("IsContiguous") ){
    is_contiguous_ = parameterList->get<bool>("IsContiguous");
  }

  data_.c.dbound = parameterList->get<double>("dbound", 0.0);
  data_.c.prefer_upper = (parameterList->get<bool>("PreferUpper", true)) ? 1 : 0;
  data_.c.print = parameterList->get<int>("print",3);
  data_.c.nmethods = parameterList->get<int>("nmethods",0);
  data_.c.useGPU = parameterList->get<int>("useGPU",1);
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Cholmod<Matrix,Vector>::getValidParameters_impl() const
{
  using std::string;
  using Teuchos::tuple;
  using Teuchos::ParameterList;
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::stringToIntegralParameterEntryValidator;

  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();


    Teuchos::RCP<EnhancedNumberValidator<int> > print_validator
      = Teuchos::rcp( new EnhancedNumberValidator<int>(0,5));

    Teuchos::RCP<EnhancedNumberValidator<int> > nmethods_validator
      = Teuchos::rcp( new EnhancedNumberValidator<int>(0,9));

    pl->set("nmethods", 0, "Specifies the number of different ordering methods to try", nmethods_validator);

    pl->set("print", 3, "Specifies the verbosity of the print statements", print_validator);

    pl->set("dbound", 0.0,
            "Specifies the smallest absolute value on the diagonal D for the LDL' factorization");


    pl->set("Equil", true, "Whether to equilibrate the system before solve");

    pl->set("PreferUpper", true,
            "Specifies whether the matrix will be " 
            "stored in upper triangular form.");

    pl->set("useGPU", -1, "1: Use GPU is 1, 0: Do not use GPU, -1: ENV CHOLMOD_USE_GPU set GPU usage.");

    pl->set("IsContiguous", true, "Whether GIDs contiguous");

    valid_params = pl;
  }

  return valid_params;
}


template <class Matrix, class Vector>
bool
Cholmod<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

  // Only the root image needs storage allocated

  Kokkos::resize(host_nzvals_view_, this->globalNumNonZeros_);
  Kokkos::resize(host_rows_view_, this->globalNumNonZeros_);
  Kokkos::resize(host_col_ptr_view_, this->globalNumRows_ + 1);

  long nnz_ret = 0;
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    TEUCHOS_TEST_FOR_EXCEPTION(this->rowIndexBase_ != this->columnIndexBase_,
              std::runtime_error,
              "Row and column maps have different indexbase ");

    if ( is_contiguous_ == true ) {
      Util::get_ccs_helper_kokkos_view<
        MatrixAdapter<Matrix>,host_value_type_array,host_ordinal_type_array,
          host_size_type_array>::do_get(this->matrixA_.ptr(),
            host_nzvals_view_, host_rows_view_,
            host_col_ptr_view_, nnz_ret, ROOTED,
            ARBITRARY,
            this->rowIndexBase_);
    }
    else {
      Util::get_ccs_helper_kokkos_view<
        MatrixAdapter<Matrix>,host_value_type_array,host_ordinal_type_array,
          host_size_type_array>::do_get(this->matrixA_.ptr(),
            host_nzvals_view_, host_rows_view_,
            host_col_ptr_view_, nnz_ret, CONTIGUOUS_AND_ROOTED,
            ARBITRARY,
            this->rowIndexBase_);
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(nnz_ret != Teuchos::as<long>(this->globalNumNonZeros_),
           std::runtime_error,
           "Did not get the expected number of non-zero vals");

  function_map::cholmod_init_sparse(Teuchos::as<size_t>(this->globalNumRows_),
            Teuchos::as<size_t>(this->globalNumCols_),
            Teuchos::as<size_t>(this->globalNumNonZeros_),
            0,
            host_col_ptr_view_.data(),
            host_nzvals_view_.data(),
            host_rows_view_.data(),
            &(data_.A));

  TEUCHOS_TEST_FOR_EXCEPTION(data_.A.stype == 0, std::runtime_error,
    "CHOLMOD loadA_impl loaded matrix but it is not symmetric.");

  return true;
}


template<class Matrix, class Vector>
const char* Cholmod<Matrix,Vector>::name = "Cholmod";
  

} // end namespace Amesos2

#endif  // AMESOS2_CHOLMOD_DEF_HPP
