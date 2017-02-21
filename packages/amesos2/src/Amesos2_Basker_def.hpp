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
	   "Do not have supported Kokkos node type for Basker");
#endif
  //std::cout << "MAKE BASKER" << std::endl;
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
  //std::cout << "DELETE BASKER" << std::endl;
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
  //TODO: SWITCH THIS ONCE CONVERSION COMPLETE!!!!
#if 0
        info = basker->Symbolic(this->globalNumRows_, 
                                this->globalNumCols_, 
                                this->globalNumNonZeros_, 
                                colptr_.getRawPtr(), 
                                rowind_.getRawPtr(), 
                                nzvals_.getRawPtr(),true);
#else
        info = basker->Symbolic(this->globalNumRows_, 
                               this->globalNumCols_, 
                               this->globalNumNonZeros_, 
                               sp_rowptr,
                               sp_colind,
                               sp_values,
                               true);
#endif
      }
      else { //follow original code path if conditions not met

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

#if 1
        auto sp_rowptr = this->matrixA_->returnRowPtr();
        auto sp_colind = this->matrixA_->returnColInd();
        auto sp_values = this->matrixA_->returnValues();
#else
        // NDE: Issues trying cast away const_ness... 
        //using const_rowptr_type = typename super_type::matrix_adapter_type::spmtx_ptr_t;
        //using nonconst_rowptr_type = typename std::remove_const<const_rowptr_type>::type;
        //nonconst_rowptr_type sp_rowptr= const_cast<nonconst_rowptr_type>( this->matrixA_->returnRowPtr );
        //nonconst_rowptr_type ncsp_rowptr = const_cast<nonconst_rowptr_type>( this->matrixA_->returnRowPtr );
        //int* sp_rowptr = static_cast<int*>(ncsp_rowptr);
//        int* sp_rowptr = static_cast<int*>(const_cast<nonconst_rowptr_type>( this->matrixA_->returnRowPtr()) ); // NDE: this forces match with Basker's expected input, but not right way to do this
        //int* sp_rowptr = static_cast<int*>(const_cast<typename super_type::matrix_adapter_type::spmtx_ptr_t>( this->matrixA_->returnRowPtr()) ); 


//        int* sp_rowptr = this->matrixA_->returnRowPtr(); 
//        const long* sp_rowptr = static_cast<const long*>(this->matrixA_->returnRowPtr()); 

        using sp_ptr_type = typename MatrixTraits<matrix_type>::sparse_ptr_type ;
        using nonconst_rowptr_type = typename std::remove_const<sp_ptr_type>::type;

        nonconst_rowptr_type sp_rowptr = const_cast<nonconst_rowptr_type>(this->matrixA_->returnRowPtr() ); //NDE - this works, but no match in Basker
//        const int* sp_rowptr = static_cast<const int*>(this->matrixA_->returnRowPtr() ); //NDE - does not work.. 
        //const long* sp_rowptr = static_cast<const long*>(this->matrixA_->returnRowPtr() ); //NDE - many type conversion problems when trying to pass to Basker
        // Can I wrap this in an ArrayView????

        int* sp_colind = static_cast<int*>( this->matrixA_->returnColInd() ); // NDE: Create Basker overload that receives const pointers
        double* sp_values = static_cast<double*>( this->matrixA_->returnValues() );
#endif
        // At this point, dealing with CRS, will need transpose
        //std::cout << "Calling Basker modded Factor...."<<std::endl;

  //TODO: SWITCH THIS ONCE CONVERSION COMPLETE!!!!
#if 1
        info = basker->Factor( this->globalNumRows_,
            this->globalNumCols_, 
            this->globalNumNonZeros_, 
            sp_rowptr,  // type does not match what Basker is expecting - tpetra returns const long unsigned int*, not int*; what is best thing to do? 
            sp_colind,
            sp_values);
#else
        info = basker->Factor(this->globalNumRows_,
                            this->globalNumCols_, 
                            this->globalNumNonZeros_, 
                            colptr_.getRawPtr(), 
                            rowind_.getRawPtr(), 
                            nzvals_.getRawPtr()); // NDE: Add defaulted parameter to Factor() to indicate transpose needed

#endif
      }
      else { // prep for calling the original way if these conditions are not met

      info = basker->Factor(this->globalNumRows_,
                            this->globalNumCols_, 
                            this->globalNumNonZeros_, 
                            colptr_.getRawPtr(), 
                            rowind_.getRawPtr(), 
                            nzvals_.getRawPtr()); // NDE: Add defaulted parameter to Factor() to indicate transpose needed
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
  // NDE: Once ready to pass CRS pointers to Basker, skip this routine and simply return 'true'

    // NDE: Purpose of this function is to copy SolverCopy matrixA_ to
    // CCS members in the Basker ConcreteSolver (in CCS format)
    // In the existing implementation, this is required for case of reusing
    // symbolic structure of a matrix but with different values in the matrix
    // i.e. after first run, skip the symbolicFactorization stage and
    // call numericFactorization followed by the solve.
    // 
    // This is required because the internal ArrayViews, after being initialized from
    // matrixA_, are dereferenced and then passed to ShyLUBasker or Basker 
    // (via raw pointer) which make a copy of the data prior to performing their own 
    // operations on/with the data.
    //
    // Amesos2 Basker, to stay consistent with design, should it get a copy of the 
    // matrixA_ data? Or can it just access/pass the data from matrixA_ directly
    // (as dereferenced raw pointers or Views) to ShyLUBasker, which then makes its
    // own copy? 
    // It does not appear that the copy to the Basker ConcreteSolver is necessary, as
    // the CRS -> CCS conversion can happen within Basker itself via matrix transpose
    // (and additional sort, if sorted row indices is desired).
    //
    // During the call to symbolic or numeric, a new matrix will be received/copied 
    // to matrixA_ (as Epetra or Tpetra CRS matrix). 
    //

    // NDE: Goal - 'short-circuit' this for SHYLUBASKER
    // pass the Tpetra CRS views (p,i,val) to Basker then transpose there
    // At this point, matrixA_ is already defined and part of SolverCore
    // nzvals_, rowind_, colptr_ are members of this Basker ConcreteSolver
    //
    // Would like to pass the Tpetra Views to ShyLUBasker directly 
    // - what are consequences of not storing/recopying from matrixA_ to cp,ri,val???
    //
    // Important to ensure these aspects are properly updated; however, this should be done in SolverCore, not here, 
    // thus the responsibility should already be properly handled...
    // 1. After preordering:
    // ++status_.numPreOrder_;
    // status_.last_phase_ = PREORDERING;
    //
    // 2. After Symbolic
    // ++status_.numSymbolicFact_;
    // status_.last_phase_ = SYMBFACT;
    //   Skipping loadA depends on status.preOrderingDone() and matrix_loaded_ values
    //   May call loadA(SYMFBACT)
    //
    // 3. After numericFactorization
    // ++status_.numNumericFact_;
    // status_.last_phase_ = NUMFACT;
    //   Skipping loadA depends on status.symbolicFactorizationDone() and matrix_loaded_ values
    //   May call loadA(NUMFACT)
    //
    // Concerns: 
    // 1. What if row index begins at 1, not 0? In do_getCrs, the col indices are adjusted - should 
    // this info be passed along to Basker to deal with this before/after transpose? Or if the starting index
    // is 1, should it be adjusted here???

// NDE: New code attempt
    // This special case will work if rank == 1, numproc = 1, and this->root_ is true...


  if ( (this->root_) && (this->matrixA_->getComm()->getRank() == 0) && (this->matrixA_->getComm()->getSize() == 1) ) {
  // NDE: Should full code get run during symbolicFactorization???
  //TODO: REMOVE ALL THIS ONCE CONVERSION COMPLETE!!!!
#if 1
#else
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
      std::cout << "  Amesos2 Basker: Cal get_ccs_helper for loadA_impl" << std::endl;
      Util::get_ccs_helper<
        MatrixAdapter<Matrix>,slu_type,local_ordinal_type,local_ordinal_type>
        ::do_get(this->matrixA_.ptr(), nzvals_(), rowind_(), colptr_(),
            nnz_ret, ROOTED, ARBITRARY); // copies from matrixA_ to Basker ConcreteSolver cp, ri, nzval members
    }

    // NDE: If skipping do_get call above, this check should be skipped as well as nnz_ret will not be updated
    if( this->root_ ){
      TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<local_ordinal_type>(this->globalNumNonZeros_),
          std::runtime_error,
          "Did not get the expected number of non-zero vals");
    }
#endif
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
#else // Not ShyLU Basker

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
             nnz_ret, ROOTED, ARBITRARY); // copies from matrixA_ to Basker ConcreteSolver cp, ri, nzval members
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
