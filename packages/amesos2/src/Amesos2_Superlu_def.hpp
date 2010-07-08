// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/*
 * symbolicFactorization_impl() performs no action other than making sure that
 * Factorization is performed.  numericFactorization_impl() performs the
 * factorization but no solve (right hand side is set to 0 vectors) reuses
 * factors only if ReuseFactorization_ is set.
 *  
 *  If FactorizationOK_() && ReuseSymbolic_ 
 *     reFactor()
 *  else
 *     factor() 
 *
 * solve() 
 *
 * factor() does everything from scratch:
 *  - Redistributes the data if necessary
 *  - Deletes any data structures left over from the previous call to
 *    Factor()
 *  - Copies the data into the format that SuperLU wants it
 *  - Calls dgssvx to factor the matrix with factor set to true
 * reFactor()
 *  - Redistributes the data if necessary
 *    + Attempting to check to make sure that the non-zero structure is
 *      unchanged
 *  - Copies the data into the format that SuperLU already has it
 *    FIRST PASS - assert( false ) 
 */

#ifndef AMESOS2_SUPERLU_DEF_HPP
#define AMESOS2_SUPERLU_DEF_HPP

#include "Amesos2_Solver.hpp"
#include "Amesos2_Superlu_MatrixHelper.hpp"
#include "Amesos2_Superlu_TypeMap.hpp"
#include "Amesos2_Superlu_FunctionMap.hpp"


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
   *
   * TODO:: Make sure the SLU matrix frees do not actually try to free raw
   * storage
   */
  
  SLU::StatFree( &(data_.stat) ) ;

  // Storage is initialized in numericFactorization_impl()
  if ( this->getNumNumericFact() > 0 ){
//    SLU::Destroy_CompCol_Matrix( &(data_.A) );
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );

    SLU::Destroy_SuperNode_Matrix( &(data_.L) );
    SLU::Destroy_CompCol_Matrix( &(data_.U) );
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


int preOrdering_impl()
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
  int info;
  if (! factorizationDone_){
    // First time factorization, factor from scratch
    data_.options.Fact = SLU::DOFACT;

    typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;
    MatrixHelper<Amesos::Superlu>::createCRSMatrix(
      this->matrixA_.ptr(),
      (Teuchos::ArrayView<slu_type>)nzvals_,
      (Teuchos::ArrayView<int>)rowind_,
      (Teuchos::ArrayView<int>)colptr_,
      Teuchos::ptrFromRef(data_.A));

    // Do factorization
    FunctionMap<Amesos::Superlu,scalar_type>::gstrf(&(data_.options), &(data_.A),
      data_.relax, data_.panel_size, data_.etree.getRawPtr(), NULL, 0,
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), &(data_.L), &(data_.U),
      &(data_.stat), &info);

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

  } else {                      // We have factored once before
    // Get values from matrix in temporary storage

    // Temporary store for nonzero values
    Teuchos::Array<typename TypeMap<Amesos::Superlu,scalar_type>::type> nzvals_tmp;
    // Temporary store for row indices
    Teuchos::Array<int> rowind_tmp;
    // Temporary store for column pointers
    Teuchos::Array<int> colptr_tmp;

    typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;
    MatrixHelper<Amesos::Superlu>::createCRSMatrix(
      this->matrixA_.ptr(),
      (Teuchos::ArrayView<slu_type>)nzvals_tmp,
      (Teuchos::ArrayView<int>)rowind_tmp,
      (Teuchos::ArrayView<int>)colptr_tmp,
      Teuchos::ptrFromRef(data_.A));

    // Test whether the structure of the matrix has changed
    bool samePattern = (rowind_tmp == rowind_) && (colptr_tmp == colptr_);

    rowind_ = rowind_tmp;
    colptr_ = colptr_tmp;

    if( samePattern ){
      data_.options.Fact = SLU::SamePattern;
    } else {
      // If symbolic structure has been changed, then we must factor from
      // scratch
      data_.options.Fact = SLU::DOFACT;
    }
    
    FunctionMap<Amesos::Superlu,scalar_type>::gstrf(&(data_.options), &(data_.A),
      data_.relax, data_.panel_size, data_.etree.getRawPtr(), NULL, 0,
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), &(data_.L), &(data_.U),
      &(data_.stat), &info);

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
  }

  /* All precesses should return the same error code */
  if( this->status_.numProcs_ != 1 ) Teuchos::broadcast(*(this->matrixA_->getComm()),0,&info);
  return(info);               
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::solve_impl()
{
  if( !factorizationDone_ ){
    factorizationDone_ = false;
    TEST_FOR_EXCEPTION( !numericFactorization_impl(),
      std::runtime_error,
      "Numeric Factorization failed!");
  }
  // Processor 0 sets up for Solve and calls SuperLU

  // TODO: Amesos2::Solver should check that the matrices are
  // properly initialized (and that X and B have equal numbers
  // of vectors) upon construction, or before calling solve().
  // This check should not be necessary within the concrete
  // solver code.

  typedef typename MatrixAdapter<Matrix>::scalar_type scalar_type;
  typedef typename TypeMap<Amesos::Superlu,scalar_type>::type slu_type;

  size_t nrhs = this->multiVecX_->getGlobalNumVectors();

  data_.ferr.resize(nrhs);
  data_.berr.resize(nrhs);

  Teuchos::ArrayRCP<slu_type> xValues;
  Teuchos::ArrayRCP<slu_type> bValues;
  // We assume the global length of the two vector has already been
  // checked for compatibility
  //size_t ldx = Teuchos::as<size_t>(this->multiVecX_->getStride());

  {
    // Get a SuperMatrix for the B and X multiVectors
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);

    xValues = MatrixHelper<Amesos::Superlu>::createMVDenseMatrix(
      this->multiVecX_.ptr(), Teuchos::ptrFromRef(data_.X));

    bValues = MatrixHelper<Amesos::Superlu>::createMVDenseMatrix(
      this->multiVecB_.ptr(), Teuchos::ptrFromRef(data_.B));
  }
  int ierr;		// returned error code
  {
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);

    typedef typename TypeMap<Amesos::Superlu,scalar_type>::magnitude_type magnitude_type;
    magnitude_type rpg, rcond;

    /* Superlu does by default a Transpose solve when the matrix A is stored
     * in compressed row format.  So, if we actually want to use no
     * transpose, we have to tell it the oposite.
     */
    if( this->control_.useTranspose_ ){
      data_.options.Trans = SLU::TRANS;
    } else {
      data_.options.Trans = SLU::NOTRANS; 
    }
    
    FunctionMap<Amesos::Superlu,scalar_type>::gssvx(&(data_.options), &(data_.A),
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), data_.etree.getRawPtr(),
      &(data_.equed), data_.R.getRawPtr(), data_.C.getRawPtr(), &(data_.L),
      &(data_.U), NULL, 0, &(data_.B), &(data_.X), &rpg, &rcond,
      data_.ferr.getRawPtr(), data_.berr.getRawPtr(), &(data_.mem_usage),
      &(data_.stat), &ierr);

    /* Update X's global values */
    {
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
      this->multiVecX_->globalize(xValues()); // operator() does conversion from ArrayRCP to ArrayView
    }
  } // end block for solve time

  /* All precesses should return the same error code */
  if( this->status_.numProcs_ != 1 ) Teuchos::broadcast(*(this->matrixA_->getComm()),0,&ierr);
  return(ierr);               
}


template <class Matrix, class Vector>
bool
Superlu<Matrix,Vector>::matrixShapeOK_impl() const
{
  // Superlu can handle square as well as rectangular matrices
  return(true);
}


template <class Matrix, class Vector>
void
Superlu<Matrix,Vector>::setParameters_impl(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{}


/* Updates the CRS form from matrix_.  It also checks whether the nonzero
 * structure has changed since the last time.
 *
 * \return \c true if the nonzero structure of matrixA_ is the same as before
 */
template <class Matrix, class Vector>
bool
Superlu<Matrix,Vector>::recreateSLUMatrixAndCheckSameStructure()
{
  using Teuchos::ArrayView;
  using Teuchos::Array;
  
  TEST_FOR_EXCEPTION( !factorizationDone_,
    std::logic_error,
    "Refactor without factorization being done!");
  Teuchos::TimeMonitor localtimer(this->timers_.mtxconvtime_);
  if (this->status_->myPID_ == 0){
    if (this->globalNumRows_ != this->matrixA_->getGlobalNumRows() ||
      this->globalNumRows_ != this->matrixA_->getGlobalNumCols())
    {
      // something fishy here.  The dimensions of the matrix have changed
      //
      // TODO: Is the above checking that the matrix is square?
      // Because SuperLU does support general rectangular
      // matrices.
      TEST_FOR_EXCEPTION(-1, std::runtime_error,
        "Unsupported matrix dimension in SuperLU factor");
    }

    Array<scalar_type> nzval(this->globalNumNonZeros_);
    Array<global_ordinal_type> colind(this->globalNumNonZeros_);
    Array<global_ordinal_type> rowptr(this->globalNumRows_ + 1);
    global_size_type nnz;

    this->matrixA_->getCrs( (ArrayView<scalar_type>)nzval,
      (ArrayView<global_ordinal_type>)colind,
      (ArrayView<global_ordinal_type>)rowptr,
      nnz);

    // Check for equality with existing, then return bool
    
    // Timer stop here for Amesos2 overhead time

    Destroy_SuperMatrix_Store(&data_.A);
    Destroy_SuperNode_Matrix(&data_.L);
    Destroy_CompCol_Matrix(&data_.U);
  }
  return(0);
}
    

template<class Matrix, class Vector>
const char* Superlu<Matrix,Vector>::name = "Amesos2::SuperLU";


} // end namespace Amesos

#endif	// AMESOS2_SUPERLU_DEF_HPP
