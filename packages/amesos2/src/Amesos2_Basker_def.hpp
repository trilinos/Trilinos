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
{

  //Nothing

  // Override some default options
  // TODO: use data_ here to init
   
  
#ifdef SHYLUBASKER
#ifdef HAVE_AMESOS2_KOKKOS
#ifdef KOKKOS_HAVE_OPENMP
  /*
  static_assert(std::is_same<kokkos_exe,Kokkos::OpenMP>::value,
  	"Kokkos node type not supported by experimental Basker Amesos2");
  */
  typedef Kokkos::OpenMP Exe_Space;
#elif defined(KOKKOS_HAVE_SERIAL)
  typedef Kokkos::Serial Exe_Space;
#else
 TEUCHOS_TEST_FOR_EXCEPTION(1 != 0,
		     std::runtime_error,
	   "Amesos2_Basker Exception: Do not have supported Kokkos node type for Basker");
#endif
  basker = new ::BaskerNS::Basker<local_ordinal_type, slu_type, Exe_Space>(); 
  basker->Options.no_pivot      = BASKER_TRUE;
  basker->Options.symmetric     = BASKER_FALSE;
  basker->Options.realloc       = BASKER_FALSE;
  basker->Options.verbose       = BASKER_FALSE;
  basker->Options.matching      = BASKER_TRUE;
  basker->Options.matching_type = 0;
  basker->Options.btf           = BASKER_TRUE;
  basker->Options.amd_btf       = BASKER_TRUE;
  basker->Options.amd_dom       = BASKER_TRUE;
  basker->Options.transpose     = BASKER_FALSE;
  basker->Options.verbose_matrix_out = BASKER_FALSE;
  num_threads = Kokkos::OpenMP::max_hardware_threads();
#endif  
#endif
  
}


template <class Matrix, class Vector>
Basker<Matrix,Vector>::~Basker( )
{  
#ifdef SHYLUBASKER
#ifdef HAVE_AMESOS2_KOKKOS
  delete basker;
#endif
#endif
 
  /* Basker will cleanup its own internal memory*/
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
  
#ifdef SHYLUBASKER

  if(this->root_)
    {     
      int info = 0;

      //std::cout << "setting number of threads " 
      //	<< num_threads << std::endl;
      basker->SetThreads(num_threads); 
      //std::cout << "Set Threads Done" << std::endl;

    #ifdef HAVE_AMESOS2_VERBOSE_DEBUG
      std::cout << "Basker:: Before symbolic factorization" << std::endl;
      std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
      std::cout << "rowind_ : " << rowind_.toString() << std::endl;
      std::cout << "colptr_ : " << colptr_.toString() << std::endl;
    #endif

      // NDE: Special case 
      // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
      // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
      //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic
      if ( (this->matrixA_->getComm()->getRank() == 0) && (this->matrixA_->getComm()->getSize() == 1) ) {

        auto sp_rowptr = this->matrixA_->returnRowPtr();
        auto sp_colind = this->matrixA_->returnColInd();
        auto sp_values = this->matrixA_->returnValues();

        // This will require mods and extra impl of Symbolic
        info = basker->Symbolic(this->globalNumRows_, 
                               this->globalNumCols_, 
                               this->globalNumNonZeros_, 
                               sp_rowptr,
                               sp_colind,
                               sp_values,
                               true);
      }
      else {   //follow original code path if conditions not met

        info =basker->Symbolic(this->globalNumRows_, 
                               this->globalNumCols_, 
                               this->globalNumNonZeros_, 
                               colptr_.getRawPtr(), 
                               rowind_.getRawPtr(), 
                               nzvals_.getRawPtr());
      }
      //std::cout << "Symbolic Factorization Done" << std::endl; 
      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
				 std::runtime_error,
				 "Error in Basker Symbolic");
 
    }
#endif
 
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

     
#ifdef SHYLUBASKER
      // NDE: Special case 
      // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
      // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
      //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic
      if ( (this->matrixA_->getComm()->getRank() == 0) && (this->matrixA_->getComm()->getSize() == 1) ) { // outer if check if this is root_

        auto sp_rowptr = this->matrixA_->returnRowPtr();
        auto sp_colind = this->matrixA_->returnColInd();
        auto sp_values = this->matrixA_->returnValues();

        info = basker->Factor( this->globalNumRows_,
            this->globalNumCols_, 
            this->globalNumNonZeros_, 
            sp_rowptr,
            sp_colind,
            sp_values);
      }
      else {

      info = basker->Factor(this->globalNumRows_,
                            this->globalNumCols_, 
                            this->globalNumNonZeros_, 
                            colptr_.getRawPtr(), 
                            rowind_.getRawPtr(), 
                            nzvals_.getRawPtr());
      //We need to handle the realloc options
      }

      //basker->DEBUG_PRINT();

      TEUCHOS_TEST_FOR_EXCEPTION(info != 0, 
				 std::runtime_error,
				 "Error Basker Factor");

#else
      info =basker.factor(this->globalNumRows_, this->globalNumCols_, this->globalNumNonZeros_, colptr_.getRawPtr(), rowind_.getRawPtr(), nzvals_.getRawPtr());
#endif

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
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  const size_t val_store_size = as<size_t>(ld_rhs * nrhs);

  xvals_.resize(val_store_size);
  bvals_.resize(val_store_size);

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif
    Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
      slu_type>::do_get(B, bvals_(),as<size_t>(ld_rhs),
                                               ROOTED);
  }

  int ierr = 0; // returned error code

  if ( this->root_ ) {
    {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

#ifdef SHYLUBASKER
      ierr = basker->Solve(nrhs, bvals_.getRawPtr(), 
			   xvals_.getRawPtr());
#else
    ierr = basker.solveMultiple(nrhs, bvals_.getRawPtr(),xvals_.getRawPtr());
#endif
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

    Util::put_1d_data_helper<
    MultiVecAdapter<Vector>,slu_type>::do_put(X, xvals_(),
                                         as<size_t>(ld_rhs),
                                         ROOTED);
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

#ifdef SHYLUBASKER
  if(parameterList->isParameter("num_threads"))
    {
      num_threads = parameterList->get<int>("num_threads");
    }
  if(parameterList->isParameter("pivot"))
    {
      basker->Options.no_pivot = (!parameterList->get<bool>("pivot"));
    }
  if(parameterList->isParameter("pivot_tol"))
    {
      basker->Options.pivot_tol = parameterList->get<double>("pivot_tol");
    }
  if(parameterList->isParameter("symmetric"))
    {
      basker->Options.symmetric = parameterList->get<bool>("symmetric");
    }
  if(parameterList->isParameter("realloc"))
    {
      basker->Options.realloc = parameterList->get<bool>("realloc");
    }
  if(parameterList->isParameter("verbose"))
    {
      basker->Options.verbose = parameterList->get<bool>("verbose");
    }
  if(parameterList->isParameter("verbose_matrix"))
    {
      basker->Options.verbose_matrix_out = parameterList->get<bool>("verbose_matrix");
    }
  if(parameterList->isParameter("matching"))
    {
      basker->Options.matching = parameterList->get<bool>("matching");
    }
  if(parameterList->isParameter("matching_type"))
    {
      basker->Options.matching_type =
	(local_ordinal_type) parameterList->get<int>("matching_type");
    }
  if(parameterList->isParameter("btf"))
    {
      basker->Options.btf = parameterList->get<bool>("btf");
    }
  if(parameterList->isParameter("amd_btf"))
    {
      basker->Options.amd_btf = parameterList->get<bool>("amd_btf");
    }
  if(parameterList->isParameter("amd_dom"))
    {
      basker->Options.amd_dom = parameterList->get<bool>("amd_dom");
    }
  if(parameterList->isParameter("transpose"))
    {
      basker->Options.transpose = parameterList->get<bool>("transpose");
    }
#endif

}

template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Basker<Matrix,Vector>::getValidParameters_impl() const
{
  using Teuchos::ParameterList;

  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;


#ifdef SHYLUBASKER
  if( is_null(valid_params) )
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
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
#else
  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    pl->set("alnnz",  2, "Approx number of nonzeros in L, default is 2*nnz(A)");
    pl->set("aunnx",  2, "Approx number of nonzeros in I, default is 2*nnz(U)");
    valid_params = pl;
  }
#endif
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

#ifdef SHYLUBASKER
  // NDE: Can clean up duplicated code with the #ifdef guards
  if ( (this->root_) && (this->matrixA_->getComm()->getRank() == 0) && (this->matrixA_->getComm()->getSize() == 1) ) {
  // NDE: Nothing is done in this special case - CRS raw pointers are passed to SHYLUBASKER and transpose of copies handled there
  }
  else {

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
      Util::get_ccs_helper<
        MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
        ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
            nnz_ret, ROOTED, ARBITRARY); // copies from matrixA_ to Basker ConcreteSolver cp, ri, nzval members
    }

    if( this->root_ ){
      TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
          std::runtime_error,
          "Did not get the expected number of non-zero vals");
    }

  } //end alternative path 
#else // Not ShyLUBasker

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
    Util::get_ccs_helper<
    MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
    ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
             nnz_ret, ROOTED, ARBITRARY);
  }

  if( this->root_ ){
    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");
  }
#endif //SHYLUBASKER
  return true;
}


template<class Matrix, class Vector>
const char* Basker<Matrix,Vector>::name = "Basker";


} // end namespace Amesos2

#endif  // AMESOS2_Basker_DEF_HPP
