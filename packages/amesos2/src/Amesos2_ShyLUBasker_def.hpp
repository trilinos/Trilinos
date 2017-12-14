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
   \file   Amesos2_ShyLUBasker_def.hpp
   \author Joshua Dennis Booth <jdbooth@sandia.gov>
           Siva Rajamanickam <srajama@sandia.gov>
           Nathan Ellingwood <ndellin@sandia.gov>

   \brief  Definitions for the Amesos2 ShyLUBasker solver interface
*/


#ifndef AMESOS2_SHYLUBASKER_DEF_HPP
#define AMESOS2_SHYLUBASKER_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_ShyLUBasker_decl.hpp"

namespace Amesos2 {


template <class Matrix, class Vector>
ShyLUBasker<Matrix,Vector>::ShyLUBasker(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::ShyLUBasker,Matrix,Vector>(A, X, B)
  , nzvals_()                   // initialize to empty arrays
  , rowind_()
  , colptr_()
  , is_contiguous_(true)
{

  //Nothing

  // Override some default options
  // TODO: use data_ here to init
   
  
#if defined(HAVE_AMESOS2_KOKKOS) && defined(KOKKOS_HAVE_OPENMP)
  /*
  static_assert(std::is_same<kokkos_exe,Kokkos::OpenMP>::value,
  	"Kokkos node type not supported by experimental ShyLUBasker Amesos2");
  */
  typedef Kokkos::OpenMP Exe_Space;

  ShyLUbasker = new ::BaskerNS::BaskerTrilinosInterface<local_ordinal_type, slu_type, Exe_Space>();
  ShyLUbasker->Options.no_pivot      = BASKER_TRUE;
  ShyLUbasker->Options.symmetric     = BASKER_FALSE;
  ShyLUbasker->Options.realloc       = BASKER_FALSE;
  ShyLUbasker->Options.verbose       = BASKER_FALSE;
  ShyLUbasker->Options.matching      = BASKER_TRUE;
  ShyLUbasker->Options.matching_type = 0;
  ShyLUbasker->Options.btf           = BASKER_TRUE;
  ShyLUbasker->Options.amd_btf       = BASKER_TRUE;
  ShyLUbasker->Options.amd_dom       = BASKER_TRUE;
  ShyLUbasker->Options.transpose     = BASKER_FALSE;
  ShyLUbasker->Options.verbose_matrix_out = BASKER_FALSE;
  num_threads = Kokkos::OpenMP::max_hardware_threads();
#else
 TEUCHOS_TEST_FOR_EXCEPTION(1 != 0,
		     std::runtime_error,
	   "Amesos2_ShyLUBasker Exception: Do not have supported Kokkos node type (OpenMP) enabled for ShyLUBasker");
#endif
}


template <class Matrix, class Vector>
ShyLUBasker<Matrix,Vector>::~ShyLUBasker( )
{  
#if defined(HAVE_AMESOS2_KOKKOS) && defined(KOKKOS_HAVE_OPENMP)
  delete ShyLUbasker;
#endif
 
  /* ShyLUBasker will cleanup its own internal memory*/
}

template <class Matrix, class Vector>
bool
ShyLUBasker<Matrix,Vector>::single_proc_optimization() const {
  return (this->root_ && (this->matrixA_->getComm()->getSize() == 1) && is_contiguous_);
}

template<class Matrix, class Vector>
int
ShyLUBasker<Matrix,Vector>::preOrdering_impl()
{
  /* TODO: Define what it means for ShyLUBasker
   */
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif

  return(0);
}


template <class Matrix, class Vector>
int
ShyLUBasker<Matrix,Vector>::symbolicFactorization_impl()
{

  if(this->root_)
  {
    int info = 0;

    ShyLUbasker->SetThreads(num_threads); 

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
    std::cout << "ShyLUBasker:: Before symbolic factorization" << std::endl;
    std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
    std::cout << "rowind_ : " << rowind_.toString() << std::endl;
    std::cout << "colptr_ : " << colptr_.toString() << std::endl;
#endif

    // NDE: Special case 
    // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
    // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
    //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic

    if ( single_proc_optimization() ) {

      // this needs to be checked during loadA_impl...
      auto sp_rowptr = this->matrixA_->returnRowPtr();
      TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null ");
      auto sp_colind = this->matrixA_->returnColInd();
      TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null ");
#ifndef HAVE_TEUCHOS_COMPLEX
      auto sp_values = this->matrixA_->returnValues();
#else
      // NDE: 09/25/2017
      // Cannot convert Kokkos::complex<T>* to std::complex<T>*; in this case, use reinterpret_cast
      using complex_type = typename Util::getStdCplxType< magnitude_type, typename matrix_adapter_type::spmtx_vals_t >::type;
      complex_type * sp_values = nullptr;
      sp_values = reinterpret_cast< complex_type * > ( this->matrixA_->returnValues() );
#endif
      TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_values returned null ");

      // In this case, colptr_, rowind_, nzvals_ are invalid
      info = ShyLUbasker->Symbolic(this->globalNumRows_,
          this->globalNumCols_,
          this->globalNumNonZeros_,
          sp_rowptr,
          sp_colind,
          sp_values,
          true);
    }
    else 
    {   //follow original code path if conditions not met
      // In this case, loadA_impl updates colptr_, rowind_, nzvals_
      info = ShyLUbasker->Symbolic(this->globalNumRows_,
          this->globalNumCols_,
          this->globalNumNonZeros_,
          colptr_.getRawPtr(),
          rowind_.getRawPtr(),
          nzvals_.getRawPtr());
    }
    TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
        std::runtime_error, "Error in ShyLUBasker Symbolic");

  } // end if (this->root_)
 
  /*No symbolic factoriztion*/
  return(0);
}


template <class Matrix, class Vector>
int
ShyLUBasker<Matrix,Vector>::numericFactorization_impl()
{
  using Teuchos::as;

  int info = 0;
  if ( this->root_ ){
    { // Do factorization
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
      std::cout << "ShyLUBasker:: Before numeric factorization" << std::endl;
      std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
      std::cout << "rowind_ : " << rowind_.toString() << std::endl;
      std::cout << "colptr_ : " << colptr_.toString() << std::endl;
#endif

      // NDE: Special case 
      // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
      // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
      //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic

      if ( single_proc_optimization() ) {

        auto sp_rowptr = this->matrixA_->returnRowPtr();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null ");
        auto sp_colind = this->matrixA_->returnColInd();
        TEUCHOS_TEST_FOR_EXCEPTION(sp_colind == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null ");
#ifndef HAVE_TEUCHOS_COMPLEX
        auto sp_values = this->matrixA_->returnValues();
#else
        // NDE: 09/25/2017
        // Cannot convert Kokkos::complex<T>* to std::complex<T>*; in this case, use reinterpret_cast
        using complex_type = typename Util::getStdCplxType< magnitude_type, typename matrix_adapter_type::spmtx_vals_t >::type;
        complex_type * sp_values = nullptr;
        sp_values = reinterpret_cast< complex_type * > ( this->matrixA_->returnValues() );
#endif
        TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_values returned null ");

        // In this case, colptr_, rowind_, nzvals_ are invalid
        info = ShyLUbasker->Factor( this->globalNumRows_,
            this->globalNumCols_,
            this->globalNumNonZeros_,
            sp_rowptr,
            sp_colind,
            sp_values);
      }
      else 
      {
        // In this case, loadA_impl updates colptr_, rowind_, nzvals_
        info = ShyLUbasker->Factor(this->globalNumRows_,
            this->globalNumCols_,
            this->globalNumNonZeros_,
            colptr_.getRawPtr(),
            rowind_.getRawPtr(),
            nzvals_.getRawPtr());
        //We need to handle the realloc options
      }

      //ShyLUbasker->DEBUG_PRINT();

      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, 
          std::runtime_error, "Error ShyLUBasker Factor");

      local_ordinal_type blnnz = local_ordinal_type(0); 
      local_ordinal_type bunnz = local_ordinal_type(0); 
      ShyLUbasker->GetLnnz(blnnz); // Add exception handling?
      ShyLUbasker->GetUnnz(bunnz);

      // This is set after numeric factorization complete as pivoting can be used;
      // In this case, a discrepancy between symbolic and numeric nnz total can occur.
      this->setNnzLU( as<size_t>( blnnz + bunnz ) );

    } // end scope for timer
  } // end if (this->root_)

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);

  //global_size_type info_st = as<global_size_type>(info);
  /* TODO : Proper error messages*/
  TEUCHOS_TEST_FOR_EXCEPTION( (info == -1) ,
    std::runtime_error,
    "ShyLUBasker: Could not alloc space for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION( (info ==  -2),
    std::runtime_error,
    "ShyLUBasker: Could not alloc needed work space");
  TEUCHOS_TEST_FOR_EXCEPTION( (info == -3) ,
    std::runtime_error,
    "ShyLUBasker: Could not alloc additional memory needed for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION( (info > 0) ,
    std::runtime_error,
    "ShyLUBasker: Zero pivot found at: " << info );

  return(info);
}


template <class Matrix, class Vector>
int
ShyLUBasker<Matrix,Vector>::solve_impl(
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
        ierr = ShyLUbasker->Solve(nrhs, b_vector, x_vector);
    } // end scope for timer

      /* All processes should have the same error code */
      Teuchos::broadcast(*(this->getComm()), 0, &ierr);

      TEUCHOS_TEST_FOR_EXCEPTION( ierr  > 0,
          std::runtime_error,
          "Encountered zero diag element at: " << ierr);
      TEUCHOS_TEST_FOR_EXCEPTION( ierr == -1,
          std::runtime_error,
          "Could not alloc needed working memory for solve" );
    } //end if (this->root_)
  } // end if ( single_proc_optimization() && nrhs == 1 )
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
        ierr = ShyLUbasker->Solve(nrhs, bvals_.getRawPtr(), xvals_.getRawPtr());
      } // end scope for timer
    } // end if (this->root_)

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
    } // end scope for timer
  } // end else

  return(ierr);
}


template <class Matrix, class Vector>
bool
ShyLUBasker<Matrix,Vector>::matrixShapeOK_impl() const
{
  // The ShyLUBasker can only handle square for right now
  return( this->globalNumRows_ == this->globalNumCols_ );
}


template <class Matrix, class Vector>
void
ShyLUBasker<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  if(parameterList->isParameter("IsContiguous"))
    {
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }

  if(parameterList->isParameter("num_threads"))
    {
      num_threads = parameterList->get<int>("num_threads");
    }
  if(parameterList->isParameter("pivot"))
    {
      ShyLUbasker->Options.no_pivot = (!parameterList->get<bool>("pivot"));
    }
  if(parameterList->isParameter("pivot_tol"))
    {
      ShyLUbasker->Options.pivot_tol = parameterList->get<double>("pivot_tol");
    }
  if(parameterList->isParameter("symmetric"))
    {
      ShyLUbasker->Options.symmetric = parameterList->get<bool>("symmetric");
    }
  if(parameterList->isParameter("realloc"))
    {
      ShyLUbasker->Options.realloc = parameterList->get<bool>("realloc");
    }
  if(parameterList->isParameter("verbose"))
    {
      ShyLUbasker->Options.verbose = parameterList->get<bool>("verbose");
    }
  if(parameterList->isParameter("verbose_matrix"))
    {
      ShyLUbasker->Options.verbose_matrix_out = parameterList->get<bool>("verbose_matrix");
    }
  if(parameterList->isParameter("matching"))
    {
      ShyLUbasker->Options.matching = parameterList->get<bool>("matching");
    }
  if(parameterList->isParameter("matching_type"))
    {
      ShyLUbasker->Options.matching_type =
        (local_ordinal_type) parameterList->get<int>("matching_type");
    }
  if(parameterList->isParameter("btf"))
    {
      ShyLUbasker->Options.btf = parameterList->get<bool>("btf");
    }
  if(parameterList->isParameter("amd_btf"))
    {
      ShyLUbasker->Options.amd_btf = parameterList->get<bool>("amd_btf");
    }
  if(parameterList->isParameter("amd_dom"))
    {
      ShyLUbasker->Options.amd_dom = parameterList->get<bool>("amd_dom");
    }
  if(parameterList->isParameter("transpose"))
    {
      ShyLUbasker->Options.transpose = parameterList->get<bool>("transpose");
    }

}

template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
ShyLUBasker<Matrix,Vector>::getValidParameters_impl() const
{
  using Teuchos::ParameterList;

  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) )
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
      pl->set("IsContiguous", true, 
	      "Are GIDs contiguous");
      pl->set("num_threads", 1, 
	      "Number of threads");
      pl->set("pivot", false,
	      "Should  not pivot");
      pl->set("pivot_tol", .0001,
	      "Tolerance before pivot, currently not used");
      pl->set("symmetric", false,
	      "Should Symbolic assume symmetric nonzero pattern");
      pl->set("realloc" , false, 
	      "Should realloc space if not enough");
      pl->set("verbose", false,
	      "Information about factoring");
      pl->set("verbose_matrix", false,
	      "Give Permuted Matrices");
      pl->set("matching", true,
	      "Use WC matching (Not Supported)");
      pl->set("matching_type", 0, 
	      "Type of WC matching (Not Supported)");
      pl->set("btf", true, 
	      "Use BTF ordering");
      pl->set("amd_btf", true,
	      "Use AMD on BTF blocks (Not Supported)");
      pl->set("amd_dom", true,
	      "Use CAMD on ND blocks (Not Supported)");
      pl->set("transpose", false,
	      "Solve the transpose A");
      valid_params = pl;
    }
  return valid_params;
}


template <class Matrix, class Vector>
bool
ShyLUBasker<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  using Teuchos::as;
  if(current_phase == SOLVE) return (false);

  #ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
  #endif


  // NDE: Can clean up duplicated code with the #ifdef guards
  if ( single_proc_optimization() ) {
  // NDE: Nothing is done in this special case - CRS raw pointers are passed to SHYLUBASKER and transpose of copies handled there
  // In this case, colptr_, rowind_, nzvals_ are invalid
  }
  else 
  {

    // Only the root image needs storage allocated
    if( this->root_ ){
      nzvals_.resize(this->globalNumNonZeros_);
      rowind_.resize(this->globalNumNonZeros_);
      colptr_.resize(this->globalNumCols_ + 1); //this will be wrong for case of gapped col ids, e.g. 0,2,4,9; num_cols = 10 ([0,10)) but num GIDs = 4...
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
              nnz_ret, ROOTED, ARBITRARY, this->rowIndexBase_); // copies from matrixA_ to ShyLUBasker ConcreteSolver cp, ri, nzval members
      }
      else {
        Util::get_ccs_helper<
          MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
          ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
              nnz_ret, CONTIGUOUS_AND_ROOTED, ARBITRARY, this->rowIndexBase_); // copies from matrixA_ to ShyLUBasker ConcreteSolver cp, ri, nzval members
      }
    }

    if( this->root_ ){
      TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
          std::runtime_error,
          "Amesos2_ShyLUBasker loadA_impl: Did not get the expected number of non-zero vals");
    }

  } //end alternative path 
  return true;
}


template<class Matrix, class Vector>
const char* ShyLUBasker<Matrix,Vector>::name = "ShyLUBasker";


} // end namespace Amesos2

#endif  // AMESOS2_SHYLUBASKER_DEF_HPP
