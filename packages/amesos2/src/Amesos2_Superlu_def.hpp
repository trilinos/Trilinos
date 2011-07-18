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
   \file   Amesos2_Superlu_def.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Thu Jul  8 22:43:51 2010

   \brief  Definitions for the Amesos2 Superlu solver interface
*/


#ifndef AMESOS2_SUPERLU_DEF_HPP
#define AMESOS2_SUPERLU_DEF_HPP


namespace Amesos2 {


template <class Matrix, class Vector>
Superlu<Matrix,Vector>::Superlu(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Superlu,Matrix,Vector>(A, X, B)
  , nzvals_(this->globalNumNonZeros_)
  , rowind_(this->globalNumNonZeros_)
  , colptr_(this->globalNumRows_ + 1)
{
  SLU::set_default_options(&(data_.options));
  // Override some default options
  data_.options.PrintStat = SLU::NO;

  SLU::StatInit(&(data_.stat));

  data_.perm_r.resize(this->globalNumRows_);
  data_.perm_c.resize(this->globalNumCols_);
  data_.etree.resize(this->globalNumCols_);
  data_.R.resize(this->globalNumRows_);
  data_.C.resize(this->globalNumCols_);

  data_.relax = SLU::sp_ienv(2); // Query optimal relax param from superlu
  data_.panel_size = SLU::sp_ienv(1); // Query optimal panel size

  data_.equed = 'N';            // No equilibration
  data_.A.Store = NULL;
}


template <class Matrix, class Vector>
Superlu<Matrix,Vector>::~Superlu( )
{
  /* Free Superlu data_types
   * - Matrices
   * - Vectors
   * - Stat object
   */
  SLU::StatFree( &(data_.stat) ) ;

  // Storage is initialized in numericFactorization_impl()
  if ( this->getNumNumericFact() > 0 ){
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );

    if ( this->root_ ){	   // only root allocated these SuperMatrices.
      SLU::Destroy_SuperNode_Matrix( &(data_.L) );
      SLU::Destroy_CompCol_Matrix( &(data_.U) );
    }
  }

  // Storage is initialized in solve_impl()
  if ( this->getNumSolve() > 0 ){
    /* Cannot use SLU::Destroy_Dense_Matrix routine here, since it attempts to
     * free the array of non-zero values, but that array has already been
     * deallocated by the MultiVector object.  So we release just the Store
     * instead.
     */
    SLU::Destroy_SuperMatrix_Store( &(data_.X) );
    SLU::Destroy_SuperMatrix_Store( &(data_.B) );
  }
}

template<class Matrix, class Vector>
int
Superlu<Matrix,Vector>::preOrdering_impl()
{
  // We need a matrix before we can pre-order it
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif
    
    matrix_helper::createCCSMatrix(
      Teuchos::ptrInArg(*this->matrixA_),
      nzvals_(), rowind_(), colptr_(),
      Teuchos::outArg(data_.A),
      this->timers_.mtxRedistTime_);
  } // end matrix conversion block

  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = NATURAL:  natural ordering
   *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
   *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
   *   permc_spec = COLAMD:   approximate minimum degree column ordering
   *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
   */
  int permc_spec = data_.options.ColPerm;
  if ( permc_spec != SLU::MY_PERMC ){
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif
    
    SLU::get_perm_c(permc_spec, &(data_.A), data_.perm_c.getRawPtr());
  }

  // Cleanup SuperMatrix A's Store, will be allocated again when new
  // values are retrieved in numericFactorization.
  SLU::Destroy_SuperMatrix_Store( &(data_.A) );
  data_.A.Store = NULL;

  return(0);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::symbolicFactorization_impl()
{
  /*
   * SuperLU performs symbolic factorization and numeric factorization
   * together, but does leave some options for reusing symbolic
   * structures that have been created on previous factorizations.  If
   * our Amesos2 user calls this function, that is an indication that
   * the symbolic structure of the matrix is no longer valid, and
   * SuperLU should do the factorization from scratch.
   *
   * This can be accomplished by setting the options.Fact flag to
   * DOFACT, as well as setting our own internal flag to false.
   */
  same_symbolic_ = false;
  data_.options.Fact = SLU::DOFACT;

  return(0);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::numericFactorization_impl()
{
  using Teuchos::as;
  
  int info = 0;
  if( data_.A.Store != NULL ){
    // Cleanup old SuperMatrix A's Store, will be allocated again when
    // new values are retrieved.
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );
    data_.A.Store = NULL;

    // Cleanup old L and U matrices if we are not reusing a symbolic
    // factorization.  Stores and other data will be allocated in
    // gstrf.  Only rank 0 has valid pointers
    if ( !same_symbolic_ && this->root_ ){
      SLU::Destroy_SuperNode_Matrix( &(data_.L) );
      SLU::Destroy_CompCol_Matrix( &(data_.U) );
    }
  }

  if( same_symbolic_ ) data_.options.Fact = SLU::SamePattern_SameRowPerm;
  
  {                           // start matrix conversion block
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif
    
    matrix_helper::createCCSMatrix(
      Teuchos::ptrInArg(*this->matrixA_),
      nzvals_(), rowind_(), colptr_(),
      Teuchos::outArg(data_.A),
      this->timers_.mtxRedistTime_);
  } // end matrix conversion block

  if ( this->root_ ){
    
#ifdef HAVE_AMESOS2_DEBUG
    TEST_FOR_EXCEPTION( data_.A.ncol != as<int>(this->globalNumCols_),
			std::runtime_error,
			"Error in converting to SuperLU SuperMatrix: wrong number of global columns." );
    TEST_FOR_EXCEPTION( data_.A.nrow != as<int>(this->globalNumRows_),
			std::runtime_error,
			"Error in converting to SuperLU SuperMatrix: wrong number of global rows." );
#endif

    if( data_.options.Equil == SLU::YES ){
      magnitude_type rowcnd, colcnd, amax;
      int info2 = 0;
      
      // calculate row and column scalings
      function_map::gsequ(&(data_.A), data_.R.getRawPtr(),
			  data_.C.getRawPtr(), &rowcnd, &colcnd,
			  &amax, &info2);
      TEST_FOR_EXCEPTION( info2 != 0,
			  std::runtime_error,
			  "SuperLU gsequ returned with status " << info2 );
      
      // apply row and column scalings if necessary
      function_map::laqgs(&(data_.A), data_.R.getRawPtr(),
			  data_.C.getRawPtr(), rowcnd, colcnd,
			  amax, &(data_.equed));
      
      // // check what types of equilibration was actually done
      // data_.rowequ = (data_.equed == 'R') || (data_.equed == 'B');
      // data_.colequ = (data_.equed == 'C') || (data_.equed == 'B');
    }
    
    // Apply the column permutation computed in preOrdering.  Place the
    // column-permuted matrix in AC
    SLU::sp_preorder(&(data_.options), &(data_.A), data_.perm_c.getRawPtr(),
		     data_.etree.getRawPtr(), &(data_.AC));
    
    { // Do factorization
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);
#endif
      
#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
      std::cout << "Superlu:: Before numeric factorization" << std::endl;
      std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
      std::cout << "rowind_ : " << rowind_.toString() << std::endl;
      std::cout << "colptr_ : " << colptr_.toString() << std::endl;
#endif
      
      function_map::gstrf(&(data_.options), &(data_.AC),
			  data_.relax, data_.panel_size, data_.etree.getRawPtr(),
			  NULL, 0, data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(),
			  &(data_.L), &(data_.U), &(data_.stat), &info);
    }
    // Cleanup. AC data will be alloc'd again for next factorization (if at all)
    SLU::Destroy_CompCol_Permuted( &(data_.AC) );

    // Set the number of non-zero values in the L and U factors
    this->setLNNZ(as<size_t>(((SLU::SCformat*)data_.L.Store)->nnz));
    this->setUNNZ(as<size_t>(((SLU::NCformat*)data_.U.Store)->nnz));
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);

  global_size_type info_st = as<global_size_type>(info);
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st <= this->globalNumCols_),
    std::runtime_error,
    "Factorization complete, but matrix is singular. Division by zero eminent");
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st > this->globalNumCols_),
    std::runtime_error,
    "Memory allocation failure in Superlu factorization");

  data_.options.Fact = SLU::FACTORED;
  same_symbolic_ = true;

  return(info);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> > X,
				   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;
  // root sets up for Solve and calls SuperLU

  const global_size_type len_rhs = X->getGlobalLength();
  const size_t nrhs = X->getGlobalNumVectors();

  data_.ferr.resize(nrhs);
  data_.berr.resize(nrhs);

  const size_t val_store_size = as<size_t>(len_rhs * nrhs);
  Teuchos::Array<slu_type> xValues(val_store_size);
  Teuchos::Array<slu_type> bValues(val_store_size);
  size_t ldx, ldb;

  // We assume the global length of the two vector has already been
  // checked for compatibility

  // Clean up old X and B stores if they have already been created
  if( this->getNumSolve() > 0 ){
    SLU::Destroy_SuperMatrix_Store( &(data_.X) );
    SLU::Destroy_SuperMatrix_Store( &(data_.B) );
  }

  {                 // Convert: Get a SuperMatrix for the B and X multi-vectors
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecConvTime_);
#endif

    matrix_helper::createMVDenseMatrix(
      X,
      xValues(),                // pass as ArrayView
      ldx,
      Teuchos::outArg(data_.X),
      this->timers_.vecRedistTime_);

    matrix_helper::createMVDenseMatrix(
      B,
      bValues(),                // pass as ArrayView
      ldb,
      Teuchos::outArg(data_.B),
      this->timers_.vecRedistTime_);

    // Note: the values of B and X (after solution) are adjusted
    // appropriately within gssvx for row and column scalings.
    
  }         // end block for conversion time

  int ierr = 0; // returned error code

  magnitude_type rpg, rcond;
  if ( this->root_ ) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
    std::cout << "Superlu:: Before solve" << std::endl;
    std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
    std::cout << "rowind_ : " << rowind_.toString() << std::endl;
    std::cout << "colptr_ : " << colptr_.toString() << std::endl;
    std::cout << "B : " << bValues().toString() << std::endl;
    std::cout << "X : " << xValues().toString() << std::endl;
#endif

    function_map::gssvx(&(data_.options), &(data_.A),
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), data_.etree.getRawPtr(),
      &(data_.equed), data_.R.getRawPtr(), data_.C.getRawPtr(), &(data_.L),
      &(data_.U), NULL, 0, &(data_.B), &(data_.X), &rpg, &rcond,
      data_.ferr.getRawPtr(), data_.berr.getRawPtr(), &(data_.mem_usage),
      &(data_.stat), &ierr);

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
    std::cout << "Superlu:: After solve" << std::endl;
    std::cout << "B : " << bValues().toString() << std::endl;
    std::cout << "X : " << xValues().toString() << std::endl;
#endif
  } // end block for solve time

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()),0,&ierr);

  global_size_type ierr_st = as<global_size_type>(ierr);
  TEST_FOR_EXCEPTION( ierr < 0,
		      std::invalid_argument,
		      "Argument " << -ierr << " to SuperLU xgssvx had illegal value" );
  TEST_FOR_EXCEPTION( ierr > 0 && ierr_st <= this->globalNumCols_,
		      std::runtime_error,
		      "Factorization complete, but U is exactly singular" );
  TEST_FOR_EXCEPTION( ierr > 0 && ierr_st > this->globalNumCols_ + 1,
		      std::runtime_error,
		      "SuperLU allocated " << ierr - this->globalNumCols_ << " bytes of "
		      "memory before allocation failure occured." );

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    Util::put_1d_data_helper<
      MultiVecAdapter<Vector> ,slu_type>::do_put(X, xValues(), ldx, ROOTED);
  }

  return(ierr);
}


template <class Matrix, class Vector>
bool
Superlu<Matrix,Vector>::matrixShapeOK_impl() const
{
  // The Superlu factorization routines can handle square as well as
  // rectangular matrices, but Superlu can only apply the solve routines to
  // square matrices, so we check the matrix for squareness.
  return( this->matrixA_->getGlobalNumRows() == this->matrixA_->getGlobalNumCols() );
}


template <class Matrix, class Vector>
void
Superlu<Matrix,Vector>::setParameters_impl(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  if( this->control_.useTranspose_ ){
    data_.options.Trans = SLU::TRANS;
  }
  // The Superlu user guide uses "Trans" as the option name, so we will honor
  // that parameter as well.  Since the Control class doesn't recognize this
  // parameter, we check for it ourselves.
  else if ( parameterList->isParameter("Trans") ){
    std::string fact = parameterList->template get<std::string>("Trans");
    if( fact == "TRANS" ){
      data_.options.Trans = SLU::TRANS;
    } else if ( fact == "NOTRANS" ){
      data_.options.Trans = SLU::NOTRANS;
    } else if ( fact == "CONJ" ) {
      data_.options.Trans = SLU::CONJ;

      // TODO: Fix this!
      TEST_FOR_EXCEPTION( fact == "CONJ" && Teuchos::ScalarTraits<scalar_type>::isComplex,
			  std::invalid_argument,
			  "Amesos2::Superlu does not currently support solution of complex "
			  "systems with conjugate transpose" );
    }
  } else {                      // default to no transpose if no parameter given
    data_.options.Trans = SLU::NOTRANS;
  }

  if( parameterList->isParameter("Equil") ){
    if ( parameterList->template isType<bool>("Equil") ){
      bool equil = parameterList->template get<bool>("Equil");
      if( equil ){
	data_.options.Equil = SLU::YES;
      } else {
	data_.options.Equil = SLU::NO;
      }
    } else if ( parameterList->template isType<std::string>("Equil") ) {
      std::string equil = parameterList->template get<std::string>("Equil");
      if ( equil == "YES" || equil == "yes" ){
	data_.options.Equil = SLU::YES;
      } else if ( equil == "NO" || equil == "no" ) {
	data_.options.Equil = SLU::NO;
      }
    }
  }

  if( parameterList->isParameter("IterRefine") ){
    std::string refine = parameterList->template get<std::string>("IterRefine");
    if( refine == "NO" ){
      data_.options.IterRefine = SLU::NOREFINE;
    } else if ( refine == "SINGLE" ) {
      data_.options.IterRefine = SLU::SINGLE;
    } else if ( refine == "DOUBLE" ) {
      data_.options.IterRefine = SLU::DOUBLE;
    } else if ( refine == "EXTRA" ) {
      data_.options.IterRefine = SLU::EXTRA;
    } else {
      TEST_FOR_EXCEPTION(
	true,
	std::invalid_argument,
	"Unrecognized value for 'IterRefine' key.");
    }
  }

  if( parameterList->isParameter("SymmetricMode") ){
    if ( parameterList->template isType<bool>("SymmetricMode") ){
      bool sym = parameterList->template get<bool>("SymmetricMode");
      if( sym ){
	data_.options.SymmetricMode = SLU::YES;
      } else {
	data_.options.SymmetricMode = SLU::NO;
      }
    } else if ( parameterList->template isType<std::string>("SymmetricMode") ) {
      std::string sym = parameterList->template get<std::string>("SymmetricMode");
      if ( sym == "YES" || sym == "yes" ){
	data_.options.SymmetricMode = SLU::YES;
      } else if ( sym == "NO" || sym == "no" ) {
	data_.options.SymmetricMode = SLU::NO;
      }
    }
  }

  if( parameterList->isParameter("DiagPivotThresh") ){
    double diag_pivot_thresh = parameterList->template get<double>("DiagPivotThresh");
    data_.options.DiagPivotThresh = diag_pivot_thresh;
    TEST_FOR_EXCEPTION( diag_pivot_thresh < 0 || diag_pivot_thresh > 1,
			std::invalid_argument,
			"Invalid value given for 'DiagPivotThresh' parameter" );
  }

  if( parameterList->isParameter("ColPerm") ){
    std::string method = parameterList->template get<std::string>("ColPerm");
    if( method == "NATURAL" ){
      data_.options.ColPerm = SLU::NATURAL;
    } else if ( method == "MMD_AT_PLUS_A" ) {
      data_.options.ColPerm = SLU::MMD_AT_PLUS_A;
    } else if ( method == "MMD_ATA" ) {
      data_.options.ColPerm = SLU::MMD_ATA;
    } else if ( method == "COLAMD" ) {
      data_.options.ColPerm = SLU::COLAMD;
    } else if ( method == "MY_PERMC" ) {
      data_.options.ColPerm = SLU::MY_PERMC;

      // Now we also expect to find a parameter in parameterList called
      // "perm_c"
      TEST_FOR_EXCEPTION(
	!parameterList->isParameter("perm_c"),
	std::invalid_argument,
	"MY_PERMC option specified without accompanying 'perm_c' parameter.");

      data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

      TEST_FOR_EXCEPTION(
	Teuchos::as<global_size_type>(data_.perm_c.size()) == this->globalNumCols_,
	std::length_error,
	"'perm_c' parameter not of correct length.");
    } else {
      TEST_FOR_EXCEPTION(
	true,
	std::invalid_argument,
	"Unrecognized value for 'ColPerm' key.");
    }
  }

  // We also recognize a lone 'perm_c' parameter, assuming that ColPerm = MY_PERMC
  if( parameterList->isParameter("perm_c") ){
    data_.options.ColPerm = SLU::MY_PERMC;
    data_.perm_c = parameterList->template get<Teuchos::Array<int> >("perm_c");

    TEST_FOR_EXCEPTION(
      Teuchos::as<global_size_type>(data_.perm_c.size()) == this->globalNumCols_,
      std::length_error,
      "'perm_c' parameter not of correct length.");
  }
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Superlu<Matrix,Vector>::getValidParameters_impl() const
{
  using Teuchos::ParameterList;

  ParameterList valid_params;

  valid_params.set("Trans","NOTRANS");
  valid_params.set("Equil",true);
  valid_params.set("IterRefine","NO");
  valid_params.set("DiagPivotThresh",1.0); // partial pivoting
  valid_params.set("ColPerm","COLAMD");
  valid_params.set("SymmetricMode",false);

  return Teuchos::rcpFromRef( valid_params );
}


template<class Matrix, class Vector>
const char* Superlu<Matrix,Vector>::name = "SuperLU";


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_DEF_HPP
