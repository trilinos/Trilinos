// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2010 Sandia Corporation
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
   \file   Amesos2_Superlumt_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Tue May 24 08:37:17 MDT 2011

   \brief  Definitions for the Amesos2 SuperLU_MT solver interface
*/

#ifndef AMESOS2_SUPERLUMT_DEF_HPP
#define AMESOS2_SUPERLUMT_DEF_HPP

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Util.hpp"


namespace SLUMT {
  /*
   * We have to declare this extern function because of a bug in the
   * SuperLU_MT header files.  In each header is declared a function
   * "extern int xsp_colorder(...)" where x is in {'s','d','c','z'},
   * but this function is never defined anywhere, so if you use the
   * function as defined, it will compile fine, but will choke during
   * linking.  No other code in SuperLU_MT actually calls these
   * functions. Instead, all other SuperLU_MT function just call
   * "sp_colorder", which is defined within the SuperLU_MT library,
   * but not declared.
   */
  extern "C" {
    int sp_colorder(SuperMatrix*,int*,superlumt_options_t*,SuperMatrix*);
  }
} // end namespace SLUMT


namespace Amesos2 {

  /*
   * Note: Although many of the type definitions for SuperLU_MT are
   * identical to those of SuperLU, we do not mix the definitions so
   * that we do not introduce messy coupling between the two
   * interfaces.  Also, there exist enough differences between the two
   * packages to merit dedicated utility functions.
   *
   * We have also taken a different approach to interfacing with
   * SuperLU_MT than with SuperLU which I believe leads to a better
   * seperation of the 4 solution steps.  We may in the future adopt a
   * similar strategy for SuperLU.
   */

template <class Matrix, class Vector>
Superlumt<Matrix,Vector>::Superlumt(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Superlumt,Matrix,Vector>(A, X, B)
  , nzvals_()			// initialize to empty arrays
  , rowind_()
  , colptr_()
{
  // Set some default parameters
  Teuchos::RCP<Teuchos::ParameterList> default_params
    = Teuchos::parameterList( *(this->getValidParameters()) );
  this->setParameters(default_params);
  data_.options.PrintStat = SLUMT::NO;
  data_.options.panel_size = SLUMT::D::sp_ienv(1); // Query optimal panel size
  data_.options.relax = SLUMT::D::sp_ienv(2); // Query optimal relax param from superlumt
  data_.options.lwork = 0;	// Use system memory for factorization

  data_.perm_r.resize(this->globalNumRows_);
  data_.perm_c.resize(this->globalNumCols_);
  data_.options.perm_r = data_.perm_r.getRawPtr();
  data_.options.perm_c = data_.perm_c.getRawPtr();

  data_.R.resize(this->globalNumRows_);
  data_.C.resize(this->globalNumCols_);

  data_.options.refact = SLUMT::NO; // initially we are not refactoring
  data_.options.fact = SLUMT::DOFACT; // Just do fact (without equil) unless user asks
  data_.equed = SLUMT::NOEQUIL;	// No equilibration has yet been performed
  data_.rowequ = false;
  data_.colequ = false;

  data_.A.Store  = NULL;
  data_.AC.Store = NULL;
  data_.BX.Store = NULL;
}


template <class Matrix, class Vector>
Superlumt<Matrix,Vector>::~Superlumt( )
{
  /* Free SuperLU_MT data_types
   * - Matrices
   * - Vectors
   * - Stat object
   */

  // Storage is initialized in numericFactorization_impl()
  if ( data_.A.Store != NULL ){
    // Our Teuchos::Array's will destroy rowind, colptr, and nzval for us
    SLUMT::D::Destroy_SuperMatrix_Store( &(data_.A) );
  }

  // Cleanup memory allocated during last call to sp_colorder if needed
  if( data_.AC.Store != NULL ){
    SLUMT::D::Destroy_CompCol_Permuted( &(data_.AC) ); // free's colbeg, colend, and Store
  }

  if ( this->status_.getNumNumericFact() > 0 && this->root_ ){ // only root allocated these Objects.
    SLUMT::D::Destroy_SuperNode_SCP( &(data_.L) );
    SLUMT::D::Destroy_CompCol_NCP( &(data_.U) );

    // memory alloc'd in sp_colorder
    free( data_.options.etree );
    free( data_.options.colcnt_h );
    free( data_.options.part_super_h );

    SLUMT::D::StatFree( &(data_.stat) ) ;
  }


  // Storage is initialized in solve_impl()
  if ( data_.BX.Store != NULL ){
    /* Cannot use SLU::Destroy_Dense_Matrix routine here, since it attempts to
     * free the array of non-zero values, but that array has already been
     * deallocated by the MultiVector object.  So we release just the Store
     * instead.
     */
    SLUMT::D::Destroy_SuperMatrix_Store( &(data_.BX) );
  }
}

template<class Matrix, class Vector>
int
Superlumt<Matrix,Vector>::preOrdering_impl()
{
  // Use either the column-ordering found in the users perm_c or the requested computed ordering
  int perm_spec = data_.options.ColPerm;
  if( perm_spec != SLUMT::MY_PERMC && this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif
    
    SLUMT::S::get_perm_c(perm_spec, &(data_.A), data_.perm_c.getRawPtr());
  }
  // Ordering will be applied directly before numeric factorization

  return(0);
}



template <class Matrix, class Vector>
int
Superlumt<Matrix,Vector>::symbolicFactorization_impl()
{
  // We assume that a call to symbolicFactorization is indicating that
  // the structure of the matrix has changed.  SuperLU doesn't allow
  // us to do a symbolic factorization directly, but we can leave a
  // flag that tells it it needs to symbolically factor later.
  data_.options.refact = SLUMT::NO;

  if( this->status_.numericFactorizationDone() && this->root_ ){
    // If we've done a numeric factorization already, then we need to
    // cleanup the old L and U. Stores and other data will be
    // allocated during numeric factorization.  Only rank 0 has valid
    // pointers
    SLUMT::D::Destroy_SuperNode_Matrix( &(data_.L) );
    SLUMT::D::Destroy_CompCol_NCP( &(data_.U) );
   }

  return(0);
}


template <class Matrix, class Vector>
int
Superlumt<Matrix,Vector>::numericFactorization_impl()
{
  using Teuchos::as;
  
#ifdef HAVE_AMESOS2_DEBUG
  const int nprocs = data_.options.nprocs;
  TEST_FOR_EXCEPTION( nprocs <= 0,
		      std::invalid_argument,
		      "The number of threads to spawn should be greater than 0." );
#endif

  int info = 0;

  if ( this->root_ ) {

    if( data_.options.fact == SLUMT::EQUILIBRATE ){
      magnitude_type rowcnd, colcnd, amax;
      int info;

      function_map::gsequ(&(data_.A), data_.R.getRawPtr(),
			  data_.C.getRawPtr(), &rowcnd, &colcnd,
			  &amax, &info);
      TEST_FOR_EXCEPTION( info != 0,
			  std::runtime_error,
			  "SuperLU_MT gsequ returned with status " << info );

      function_map::laqgs(&(data_.A), data_.R.getRawPtr(),
			  data_.C.getRawPtr(), rowcnd, colcnd,
			  amax, &(data_.equed));

      data_.rowequ = (data_.equed == SLUMT::ROW) || (data_.equed == SLUMT::BOTH);
      data_.colequ = (data_.equed == SLUMT::COL) || (data_.equed == SLUMT::BOTH);

      data_.options.fact = SLUMT::DOFACT;
    }

    // Cleanup memory allocated during last call to sp_colorder if needed
    if( data_.AC.Store != NULL ){
      SLUMT::D::Destroy_CompCol_Permuted( &(data_.AC) ); // free's colbeg, colend, and Store
      if( data_.options.refact == SLUMT::NO ){		 // then we recompute etree; free the old one
	free( data_.options.etree );
	free( data_.options.colcnt_h );
	free( data_.options.part_super_h );
      }
      data_.AC.Store = NULL;
    }

    // Apply the column ordering, so that AC is the column-permuted A, and compute etree
    SLUMT::sp_colorder(&(data_.A), data_.perm_c.getRawPtr(),
			  &(data_.options), &(data_.AC));


    // Allocate and initialize status variable
    if ( this->status_.getNumNumericFact() > 0 ) SLUMT::D::StatFree( &(data_.stat) );
    const int n = as<int>(this->globalNumCols_); // n is the number of columns in A
    SLUMT::D::StatAlloc(n, data_.options.nprocs, data_.options.panel_size,
			data_.options.relax, &(data_.stat));
    SLUMT::D::StatInit(n, data_.options.nprocs, &(data_.stat));


    { // Do factorization
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer( this->timers_.numFactTime_ );
#endif

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
      std::cout << "SuperLU_MT:: Before numeric factorization" << std::endl;
      std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
      std::cout << "rowind_ : " << rowind_.toString() << std::endl;
      std::cout << "colptr_ : " << colptr_.toString() << std::endl;
#endif

      function_map::gstrf(&(data_.options), &(data_.AC),
			  data_.perm_r.getRawPtr(), &(data_.L), &(data_.U),
			  &(data_.stat), &info);
    }

    // Set the number of non-zero values in the L and U factors
    this->setLNNZ(as<size_t>(((SLUMT::SCformat*)data_.L.Store)->nnz));
    this->setUNNZ(as<size_t>(((SLUMT::NCformat*)data_.U.Store)->nnz));
  }

  // Check output
  const global_size_type info_st = as<global_size_type>(info);
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st <= this->globalNumCols_),
		      std::runtime_error,
		      "Factorization complete, but matrix is singular. Division by zero eminent");
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st > this->globalNumCols_),
		      std::runtime_error,
		      "Memory allocation failure in SuperLU_MT factorization");
  // The other option, that info_st < 0 denotes invalid parameters to
  // the function, but we'll assume for now that that won't happen.

  data_.options.fact = SLUMT::FACTORED;
  data_.options.refact = SLUMT::YES;

  /* All processes should return the same error code */
  Teuchos::broadcast(*(this->getComm()),0,&info);
  return(info);
}


template <class Matrix, class Vector>
int
Superlumt<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
				     const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  Teuchos::Array<slu_type> bxvals_(ld_rhs * nrhs);
  size_t ldbx_ = as<size_t>(ld_rhs);

  {				// Get values from B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer( this->timers_.vecConvTime_ );
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

    Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
                             slu_type>::do_get(B, bxvals_(),
					       ldbx_,
					       ROOTED);
  }

  int info = 0; // returned error code (0 = success)
  if ( this->root_ ) {
    // Clean up old B stores if it has already been created
    if( data_.BX.Store != NULL ){
      SLUMT::D::Destroy_SuperMatrix_Store( &(data_.BX) );
      data_.BX.Store = NULL;
    }

    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer( this->timers_.vecConvTime_ );
#endif
      SLUMT::Dtype_t dtype = type_map::dtype;
      function_map::create_Dense_Matrix(&(data_.BX), as<int>(ld_rhs), as<int>(nrhs),
					bxvals_.getRawPtr(), as<int>(ldbx_),
					SLUMT::SLU_DN, dtype, SLUMT::SLU_GE);
    }

    if( data_.options.trans == SLUMT::NOTRANS ){
      if( data_.rowequ ){		// row equilibration has been done on AC
	// scale bxvals_ by diag(R)
	Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.R(),
		    SLUMT::slu_mt_mult<slu_type,magnitude_type>());
      }
    } else if( data_.colequ ){	// column equilibration has been done on AC
      // scale bxvals_ by diag(C)
      Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.C(),
		  SLUMT::slu_mt_mult<slu_type,magnitude_type>());
    }


    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer( this->timers_.solveTime_ );
#endif

      function_map::gstrs(data_.options.trans, &(data_.L),
			  &(data_.U), data_.perm_r.getRawPtr(),
			  data_.perm_c.getRawPtr(), &(data_.BX),
			  &(data_.stat), &info);
    }
  } // end block for solve time

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()),0,&info);
  
  TEST_FOR_EXCEPTION( info < 0,
		      std::runtime_error,
		      "Argument " << -info << " to gstrs had an illegal value" );

  // "Un-scale" the solution so that it is a solution of the original system
  if( data_.options.trans == SLUMT::NOTRANS ){
    if( data_.colequ ){	// column equilibration has been done on AC
      // scale bxvals_ by diag(C)
      Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.C(),
		  SLUMT::slu_mt_mult<slu_type,magnitude_type>());
    }
  } else if( data_.rowequ ){		// row equilibration has been done on AC
    // scale bxvals_ by diag(R)
    Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.R(),
		SLUMT::slu_mt_mult<slu_type,magnitude_type>());
  }

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

    Util::put_1d_data_helper<
      MultiVecAdapter<Vector>, slu_type>::do_put(X, bxvals_(), ldbx_, ROOTED);
  }

  return(info);
}


template <class Matrix, class Vector>
bool
Superlumt<Matrix,Vector>::matrixShapeOK_impl() const
{
  // The SuperLU_MT factorization routines can handle square as well as
  // rectangular matrices, but SuperLU_MT can only apply the solve routines to
  // square matrices, so we check the matrix for squareness.
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
Superlumt<Matrix,Vector>::setParameters_impl(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::as;

  if( parameterList->isParameter("nprocs") ){
    const int nprocs = parameterList->template get<int>("nprocs");
    data_.options.nprocs = nprocs;
  }

  if( this->control_.useTranspose_ ){
    data_.options.trans = SLUMT::TRANS;
  }
  // The Superlu user guide uses "trans" as the option name for
  // SuperLU_MT, so we will honor that parameter as well.  Since the
  // Control class doesn't recognize this parameter, we check for it
  // ourselves.
  else if ( parameterList->isParameter("trans") ){
    const std::string fact = parameterList->template get<std::string>("trans");
    if( fact == "TRANS" ){
      data_.options.trans = SLUMT::TRANS;
    } else if ( fact == "NOTRANS" ){
      data_.options.trans = SLUMT::NOTRANS;
    } else if ( fact == "CONJ" ) {
      data_.options.trans = SLUMT::CONJ;

      // TODO: Fix this!
      TEST_FOR_EXCEPTION( fact == "CONJ" && Teuchos::ScalarTraits<scalar_type>::isComplex,
			  std::invalid_argument,
			  "Amesos2::Superlumt does not currently support solution of complex "
			  "systems with conjugate transpose");
    }
  } else {                      // default to no transpose if no parameter given
    data_.options.trans = SLUMT::NOTRANS;
  }

  if( parameterList->isParameter("panel_size") ){
    const int panel_size = parameterList->template get<int>("panel_size");
    data_.options.panel_size = panel_size;
  }

  if( parameterList->isParameter("relax") ){
    const int relax = parameterList->template get<int>("relax");
    data_.options.relax = relax;
  }

  if( parameterList->isParameter("Equil") ){
    if ( parameterList->template isType<bool>("Equil") ){
      const bool equil = parameterList->template get<bool>("Equil");
      if( equil ){
	data_.options.fact = SLUMT::EQUILIBRATE;
      } else {
	data_.options.fact = SLUMT::DOFACT;
      }
    } else if ( parameterList->template isType<std::string>("Equil") ) {
      const std::string equil = parameterList->template get<std::string>("Equil");
      if ( equil == "YES" || equil == "yes" ){
	data_.options.fact = SLUMT::EQUILIBRATE;
      } else if ( equil == "NO" || equil == "no" ) {
	data_.options.fact = SLUMT::DOFACT;
      }
    }
  }

  if( parameterList->isParameter("SymmetricMode") ){
    if ( parameterList->template isType<bool>("SymmetricMode") ){
      const bool sym = parameterList->template get<bool>("SymmetricMode");
      if( sym ){
	data_.options.SymmetricMode = SLUMT::YES;
      } else {
	data_.options.SymmetricMode = SLUMT::NO;
      }
    } else if ( parameterList->template isType<std::string>("SymmetricMode") ) {
      const std::string sym = parameterList->template get<std::string>("SymmetricMode");
      if ( sym == "YES" || sym == "yes" ){
	data_.options.SymmetricMode = SLUMT::YES;
      } else if ( sym == "NO" || sym == "no" ) {
	data_.options.SymmetricMode = SLUMT::NO;
      }
    }
  }

  if( parameterList->isParameter("PrintStat") ){
    if ( parameterList->template isType<bool>("PrintStat") ){
      const bool ps = parameterList->template get<bool>("PrintStat");
      if( ps ){
	data_.options.PrintStat = SLUMT::YES;
      } else {
	data_.options.PrintStat = SLUMT::NO;
      }
    } else if ( parameterList->template isType<std::string>("PrintStat") ) {
      std::string ps = parameterList->template get<std::string>("PrintStat");
      if ( ps == "YES" || ps == "yes" ){
	data_.options.PrintStat = SLUMT::YES;
      } else if ( ps == "NO" || ps == "no" ) {
	data_.options.PrintStat = SLUMT::NO;
      }
    }
  }

  if( parameterList->isParameter("diag_pivot_thresh") ){
    data_.options.diag_pivot_thresh = parameterList->template get<double>("diag_pivot_thresh");
  }

  if( parameterList->isParameter("ColPerm") ){
    const std::string method = parameterList->template get<std::string>("ColPerm");
    if( method == "NATURAL" ){
      data_.options.ColPerm = SLUMT::NATURAL;
    } else if ( method == "MMD_AT_PLUS_A" ) {
      data_.options.ColPerm = SLUMT::MMD_AT_PLUS_A;
    } else if ( method == "MMD_ATA" ) {
      data_.options.ColPerm = SLUMT::MMD_ATA;
    } else if ( method == "COLAMD" ) {
      data_.options.ColPerm = SLUMT::COLAMD;
    } else if ( method == "METIS_AT_PLUS_A" ) {
      data_.options.ColPerm = SLUMT::METIS_AT_PLUS_A;
    } else if ( method == "PARMETIS" ) {
      data_.options.ColPerm = SLUMT::PARMETIS;
    } else if ( method == "MY_PERMC" ) {
      data_.options.ColPerm = SLUMT::MY_PERMC;

      // Now we also expect to find a parameter in parameterList called
      // "perm_c"
      TEST_FOR_EXCEPTION(
	!parameterList->isParameter("perm_c"),
	std::invalid_argument,
	"MY_PERMC option specified without accompanying 'perm_c' parameter.");

      data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

      TEST_FOR_EXCEPTION(
	as<global_size_type>(data_.perm_c.size()) != this->globalNumCols_,
	std::length_error,
	"'perm_c' parameter not of correct length.");
    } else {
      TEST_FOR_EXCEPTION(
	true,
	std::invalid_argument,
	"Unrecognized value for 'ColPerm' key.");
    }
  } else {
    data_.options.ColPerm = SLUMT::COLAMD;
  }

  // We also recognize a lone 'perm_c' parameter, assuming that ColPerm = MY_PERMC
  if( parameterList->isParameter("perm_c") ){
    data_.options.ColPerm = SLUMT::MY_PERMC;
    data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

    TEST_FOR_EXCEPTION(
      as<global_size_type>(data_.perm_c.size()) == this->globalNumCols_,
      std::length_error,
      "'perm_c' parameter not of correct length.");
  }

  if( parameterList->isParameter("usepr") ){
    if ( parameterList->template isType<bool>("usepr") ){
      bool upr = parameterList->template get<bool>("usepr");
      if( upr ){
	data_.options.usepr = SLUMT::YES;
      } else {
	data_.options.usepr = SLUMT::NO;
      }
    } else if ( parameterList->template isType<std::string>("SymmetricMode") ) {
      std::string upr = parameterList->template get<std::string>("SymmetricMode");
      if ( upr == "YES" || upr == "yes" ){
	data_.options.SymmetricMode = SLUMT::YES;
      } else if ( upr == "NO" || upr == "no" ) {
	data_.options.SymmetricMode = SLUMT::NO;
      }
    }
    // Now if usepr == YES, then there must also be a perm_r parameter
    if( data_.options.usepr == SLUMT::YES ){
      TEST_FOR_EXCEPTION( !parameterList->isParameter("perm_r"),
			  std::invalid_argument,
			  "Must provide a 'perm_r' parameter if 'usepr' is true");

      data_.perm_r = parameterList->template get<Teuchos::Array<int> >("perm_r");

      TEST_FOR_EXCEPTION(
	as<global_size_type>(data_.perm_r.size()) != this->globalNumRows_,
	std::length_error,
	"'perm_r' parameter array not of the correct length."
	"  Should be same as the number of global rows.");
    }
  }
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Superlumt<Matrix,Vector>::getValidParameters_impl() const
{
  using Teuchos::ParameterList;

  Teuchos::RCP<ParameterList> valid_params = Teuchos::rcp(new ParameterList());

  valid_params->set("nprocs", 1);
  valid_params->set("Trans", "NOTRANS");
  valid_params->set("Equil", true);
  valid_params->set("DiagPivotThresh", 1.0); // partial pivoting
  valid_params->set("ColPerm", "COLAMD");
  valid_params->set("SymmetricMode", false);
  valid_params->set("usepr", false);
  valid_params->set("PrintStat", false);

  return valid_params;
}


template <class Matrix, class Vector>
bool
Superlumt<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  using Teuchos::as;

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer( this->timers_.mtxConvTime_ );
#endif

  if( current_phase == SYMBFACT ) return false;

  // Store is allocated on create_CompCol_Matrix
  if( data_.A.Store != NULL ){
    SLUMT::D::Destroy_SuperMatrix_Store( &(data_.A) );
    data_.A.Store = NULL;
  }

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

    // Use compressed-column store for SuperLU_MT
    Util::get_ccs_helper<
      MatrixAdapter<Matrix>,slu_type,int,int>::do_get(this->matrixA_.ptr(),
						      nzvals_, rowind_, colptr_,
						      nnz_ret, ROOTED, ARBITRARY);
  }

  if( this->root_ ){
    TEST_FOR_EXCEPTION( nnz_ret != as<int>(this->globalNumNonZeros_),
			std::runtime_error,
			"rank=0 failed to get all nonzero values in getCcs()");

    SLUMT::Dtype_t dtype = type_map::dtype;
    function_map::create_CompCol_Matrix(&(data_.A),
					as<int>(this->globalNumRows_),
					as<int>(this->globalNumCols_),
					nnz_ret,
					nzvals_.getRawPtr(),
					rowind_.getRawPtr(),
					colptr_.getRawPtr(),
					SLUMT::SLU_NC,
					dtype, SLUMT::SLU_GE);

    TEST_FOR_EXCEPTION( data_.A.Store == NULL,
			std::runtime_error,
			"Failed to allocate matrix store" );
  }

  return true;
}


template<class Matrix, class Vector>
const char* Superlumt<Matrix,Vector>::name = "SuperLU_MT";


} // end namespace Amesos2

#endif	// AMESOS2_SUPERLUMT_DEF_HPP
