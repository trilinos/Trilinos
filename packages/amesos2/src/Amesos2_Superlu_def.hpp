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


namespace Amesos {


template <class Matrix, class Vector>
Superlu<Matrix,Vector>::Superlu(
  Teuchos::RCP<Matrix> A,
  Teuchos::RCP<Vector> X,
  Teuchos::RCP<Vector> B)
  : Solver<Amesos::Superlu,Matrix,Vector>(A, X, B)
  , nzvals_(this->globalNumNonZeros_)
  , rowind_(this->globalNumNonZeros_)
  , colptr_(this->globalNumRows_ + 1)
  , factorizationDone_(false)
{
  SLU::set_default_options(&(data_.options));
  // Override some default options
  data_.options.PrintStat = SLU::NO;

  SLU::StatInit(&(data_.stat));

  data_.perm_r.resize(this->globalNumRows_);
  data_.perm_c.resize(this->globalNumRows_);
  data_.etree.resize(this->globalNumRows_);
  data_.R.resize(this->globalNumRows_);
  data_.C.resize(this->globalNumRows_);

  data_.relax = SLU::sp_ienv(2); // Query optimal relax param from superlu
  data_.panel_size = SLU::sp_ienv(1); // Query optimal panel size

  data_.equed = 'N';            // No equilibration
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

    if ( this->status_.root_ ){       // only root allocated these SuperMatrices.
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
  return(0);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::symbolicFactorization_impl()
{
  // Nothing happens here, all is done in
  // numericFactorization_impl()
  factorizationDone_ = false;
  return(0);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::numericFactorization_impl(){
  int info = 0;
  if( !factorizationDone_ ){
    // First time factorization, factor from scratch
    data_.options.Fact = SLU::DOFACT;

    typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;
    {                           // start matrix conversion block
      Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);

      MatrixHelper<Amesos::Superlu>::createCRSMatrix(
        this->matrixA_.ptr(),
        (Teuchos::ArrayView<slu_type>)nzvals_,
        (Teuchos::ArrayView<int>)rowind_,
        (Teuchos::ArrayView<int>)colptr_,
        Teuchos::ptrFromRef(data_.A),
        this->timers_.mtxRedistTime_);
    } // end matrix conversion block
  } else {                      // Factorization has been performed already
    // Cleanup old SuperMatrix A's Store, will be allocated again when
    // new values are retrieved.
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );

    // We have factored once before, check new values for same pattern

    // Get values from matrix in temporary storage
    typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;

    /*
     * Keep a copy of the old values for temporary safekeeping.  We will check
     * later whether the new values have the same non-zero structure as the
     * old.
     */
    // Temporary store for nonzero values
    Teuchos::Array<slu_type> nzvals_tmp(nzvals_);
    // Temporary store for row indices
    Teuchos::Array<int> rowind_tmp(rowind_);
    // Temporary store for column pointers
    Teuchos::Array<int> colptr_tmp(colptr_);

    {
      Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);

      MatrixHelper<Amesos::Superlu>::createCRSMatrix(
        this->matrixA_.ptr(),
        (Teuchos::ArrayView<slu_type>)nzvals_,
        (Teuchos::ArrayView<int>)rowind_,
        (Teuchos::ArrayView<int>)colptr_,
        Teuchos::ptrFromRef(data_.A),
        this->timers_.mtxRedistTime_);
    }
    // Test whether the structure of the matrix has changed
    bool samePattern = (rowind_tmp == rowind_) && (colptr_tmp == colptr_);

    if( samePattern ){
      // If symbolic structure is the same, we can reuse the old pattern.  The
      // existing L and U SuperMatrices will be accessed.
      // SamePattern_SameRowPerm assumes that the non-zero structure is the
      // same and that the numeric values themselves are similar.  If a
      // pivoting threshold is exceeded, then gstrf may actually overwrite
      // perm_c and perm_r in favour of a better permutation.
      data_.options.Fact = SLU::SamePattern_SameRowPerm;
    } else {
      // If symbolic structure has been changed, then we must factor from
      // scratch
      data_.options.Fact = SLU::DOFACT;
      // Cleanup old L and U, Stores and other data will be allocated in gstrf.
      // Only rank 0 has valid pointers
      if ( this->status_.root_ ){
        SLU::Destroy_SuperNode_Matrix( &(data_.L) );
        SLU::Destroy_CompCol_Matrix( &(data_.U) );
      }
    }
  }

  if ( this->status_.root_ ) { // Do factorization
    Teuchos::TimeMonitor numFactTimer(this->timers_.numFactTime_);

    // std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
    // std::cout << "rowind_ : " << rowind_.toString() << std::endl;
    // std::cout << "colptr_ : " << colptr_.toString() << std::endl;

    FunctionMap<Amesos::Superlu,scalar_type>::gstrf(&(data_.options), &(data_.A),
      data_.relax, data_.panel_size, data_.etree.getRawPtr(), NULL, 0,
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), &(data_.L), &(data_.U),
      &(data_.stat), &info);
  }

  // Check output
  global_size_type info_st = Teuchos::as<global_size_type>(info);
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st <= this->globalNumCols_),
    std::runtime_error,
    "Factorization complete, but matrix is singular. Division by zero eminent");
  TEST_FOR_EXCEPTION( (info_st > 0) && (info_st > this->globalNumCols_),
    std::runtime_error,
    "Memory allocation failure in Superlu factorization");

  factorizationDone_ = true;
  data_.options.Fact = SLU::FACTORED;

  /* All processes should return the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()),0,&info);
  return(info);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::solve_impl()
{
  if( !factorizationDone_ ){
    TEST_FOR_EXCEPTION( !numericFactorization_impl(),
      std::runtime_error,
      "Numeric Factorization failed!");
  }
  // root sets up for Solve and calls SuperLU

  typedef typename MatrixAdapter<Matrix>::scalar_type scalar_type;
  typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;

  size_t nrhs = this->multiVecX_->getGlobalNumVectors();

  data_.ferr.resize(nrhs);
  data_.berr.resize(nrhs);

  Teuchos::ArrayRCP<slu_type> xValues;
  Teuchos::ArrayRCP<slu_type> bValues;
  // We assume the global length of the two vector has already been
  // checked for compatibility

  // Clean up old X and B stores if they have already been created
  if( this->getNumSolve() > 0 ){
    SLU::Destroy_SuperMatrix_Store( &(data_.X) );
    SLU::Destroy_SuperMatrix_Store( &(data_.B) );
  }

  {                 // Convert: Get a SuperMatrix for the B and X multi-vectors
    Teuchos::TimeMonitor redistTimer(this->timers_.vecConvTime_);

    xValues = MatrixHelper<Amesos::Superlu>::createMVDenseMatrix(
      this->multiVecX_.ptr(), Teuchos::ptrFromRef(data_.X),
      this->timers_.vecRedistTime_);

    bValues = MatrixHelper<Amesos::Superlu>::createMVDenseMatrix(
      this->multiVecB_.ptr(), Teuchos::ptrFromRef(data_.B),
      this->timers_.vecRedistTime_);
  }         // end block for conversion time
  int ierr = 0; // returned error code

  typedef typename TypeMap<Amesos::Superlu,scalar_type>::magnitude_type magnitude_type;
  magnitude_type rpg, rcond;
  if ( this->status_.root_ ) {
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);

    // std::cout << "nzvals_ : " << nzvals_.toString() << std::endl;
    // std::cout << "rowind_ : " << rowind_.toString() << std::endl;
    // std::cout << "colptr_ : " << colptr_.toString() << std::endl;
    // std::cout << "B : " << bValues().toString() << std::endl;
    // std::cout << "X : " << xValues().toString() << std::endl;

    FunctionMap<Amesos::Superlu,scalar_type>::gssvx(&(data_.options), &(data_.A),
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), data_.etree.getRawPtr(),
      &(data_.equed), data_.R.getRawPtr(), data_.C.getRawPtr(), &(data_.L),
      &(data_.U), NULL, 0, &(data_.B), &(data_.X), &rpg, &rcond,
      data_.ferr.getRawPtr(), data_.berr.getRawPtr(), &(data_.mem_usage),
      &(data_.stat), &ierr);

    // std::cout << "B : " << bValues().toString() << std::endl;
    // std::cout << "X : " << xValues().toString() << std::endl;
  } // end block for solve time

  /* Update X's global values */
  {
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
    // broadcast solution to everyone
    Teuchos::broadcast(*(this->matrixA_->getComm()),0,xValues());

    this->multiVecX_->globalize(xValues()); // operator() does conversion from ArrayRCP to ArrayView
  }

  /* All processes should return the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()),0,&ierr);
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
    data_.options.DiagPivotThresh = parameterList->template get<double>("DiagPivotThresh");
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
const char* Superlu<Matrix,Vector>::name = "Amesos2::SuperLU";


} // end namespace Amesos

#endif	// AMESOS2_SUPERLU_DEF_HPP
