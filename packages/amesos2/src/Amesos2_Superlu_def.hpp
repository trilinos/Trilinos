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
   \file   Amesos2_Superlu_def.hpp
   \author Eric Bavier <etbavier@sandia.gov>
   \date   Thu Jul  8 22:43:51 2010

   \brief  Definitions for the Amesos2 Superlu solver interface
*/


#ifndef AMESOS2_SUPERLU_DEF_HPP
#define AMESOS2_SUPERLU_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Superlu_decl.hpp"

namespace Amesos2 {


template <class Matrix, class Vector>
Superlu<Matrix,Vector>::Superlu(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Superlu,Matrix,Vector>(A, X, B)
  , nzvals_()                   // initialize to empty arrays
  , rowind_()
  , colptr_()
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
  data_.L.Store = NULL;
  data_.U.Store = NULL;
  data_.X.Store = NULL;
  data_.B.Store = NULL;
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
  if ( data_.A.Store != NULL ){
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );
  }

  // only root allocated these SuperMatrices.
  if ( data_.L.Store != NULL ){	// will only be true for this->root_
    SLU::Destroy_SuperNode_Matrix( &(data_.L) );
    SLU::Destroy_CompCol_Matrix( &(data_.U) );
  }
}

template<class Matrix, class Vector>
int
Superlu<Matrix,Vector>::preOrdering_impl()
{
  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = NATURAL:  natural ordering
   *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
   *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
   *   permc_spec = COLAMD:   approximate minimum degree column ordering
   *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
   */
  int permc_spec = data_.options.ColPerm;
  if ( permc_spec != SLU::MY_PERMC && this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif

    SLU::get_perm_c(permc_spec, &(data_.A), data_.perm_c.getRawPtr());
  }

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

  // Cleanup old L and U matrices if we are not reusing a symbolic
  // factorization.  Stores and other data will be allocated in gstrf.
  // Only rank 0 has valid pointers
  if ( !same_symbolic_ && data_.L.Store != NULL ){
    SLU::Destroy_SuperNode_Matrix( &(data_.L) );
    SLU::Destroy_CompCol_Matrix( &(data_.U) );
    data_.L.Store = NULL;
    data_.U.Store = NULL;
  }

  if( same_symbolic_ ) data_.options.Fact = SLU::SamePattern_SameRowPerm;

  int info = 0;
  if ( this->root_ ){

#ifdef HAVE_AMESOS2_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( data_.A.ncol != as<int>(this->globalNumCols_),
                        std::runtime_error,
                        "Error in converting to SuperLU SuperMatrix: wrong number of global columns." );
    TEUCHOS_TEST_FOR_EXCEPTION( data_.A.nrow != as<int>(this->globalNumRows_),
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
      TEUCHOS_TEST_FOR_EXCEPTION( info2 != 0,
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
    this->setNnzLU( as<size_t>(((SLU::SCformat*)data_.L.Store)->nnz +
                               ((SLU::NCformat*)data_.U.Store)->nnz) );
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);

  global_size_type info_st = as<global_size_type>(info);
  TEUCHOS_TEST_FOR_EXCEPTION( (info_st > 0) && (info_st <= this->globalNumCols_),
    std::runtime_error,
    "Factorization complete, but matrix is singular. Division by zero eminent");
  TEUCHOS_TEST_FOR_EXCEPTION( (info_st > 0) && (info_st > this->globalNumCols_),
    std::runtime_error,
    "Memory allocation failure in Superlu factorization");

  data_.options.Fact = SLU::FACTORED;
  same_symbolic_ = true;

  return(info);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  const size_t val_store_size = as<size_t>(ld_rhs * nrhs);
  Teuchos::Array<slu_type> xValues(val_store_size);
  Teuchos::Array<slu_type> bValues(val_store_size);

  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif
    Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
                             slu_type>::do_get(B, bValues(),
                                               as<size_t>(ld_rhs),
                                               ROOTED, this->rowIndexBase_);
  }

  int ierr = 0; // returned error code

  magnitude_type rpg, rcond;
  if ( this->root_ ) {
    data_.ferr.resize(nrhs);
    data_.berr.resize(nrhs);

    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
#endif
      SLU::Dtype_t dtype = type_map::dtype;
      int i_ld_rhs = as<int>(ld_rhs);
      function_map::create_Dense_Matrix(&(data_.B), i_ld_rhs, as<int>(nrhs),
                                        bValues.getRawPtr(), i_ld_rhs,
                                        SLU::SLU_DN, dtype, SLU::SLU_GE);
      function_map::create_Dense_Matrix(&(data_.X), i_ld_rhs, as<int>(nrhs),
                                        xValues.getRawPtr(), i_ld_rhs,
                                        SLU::SLU_DN, dtype, SLU::SLU_GE);
    }

    // Note: the values of B and X (after solution) are adjusted
    // appropriately within gssvx for row and column scalings.

    {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

    function_map::gssvx(&(data_.options), &(data_.A),
      data_.perm_c.getRawPtr(), data_.perm_r.getRawPtr(), data_.etree.getRawPtr(),
      &(data_.equed), data_.R.getRawPtr(), data_.C.getRawPtr(), &(data_.L),
      &(data_.U), NULL, 0, &(data_.B), &(data_.X), &rpg, &rcond,
      data_.ferr.getRawPtr(), data_.berr.getRawPtr(), &(data_.mem_usage),
      &(data_.stat), &ierr);
    }

    // Cleanup X and B stores
    SLU::Destroy_SuperMatrix_Store( &(data_.X) );
    SLU::Destroy_SuperMatrix_Store( &(data_.B) );
    data_.X.Store = NULL;
    data_.B.Store = NULL;
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  global_size_type ierr_st = as<global_size_type>(ierr);
  TEUCHOS_TEST_FOR_EXCEPTION( ierr < 0,
                      std::invalid_argument,
                      "Argument " << -ierr << " to SuperLU xgssvx had illegal value" );
  TEUCHOS_TEST_FOR_EXCEPTION( ierr > 0 && ierr_st <= this->globalNumCols_,
                      std::runtime_error,
                      "Factorization complete, but U is exactly singular" );
  TEUCHOS_TEST_FOR_EXCEPTION( ierr > 0 && ierr_st > this->globalNumCols_ + 1,
                      std::runtime_error,
                      "SuperLU allocated " << ierr - this->globalNumCols_ << " bytes of "
                      "memory before allocation failure occured." );

  /* Update X's global values */
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,slu_type>::do_put(X, xValues(),
                                         as<size_t>(ld_rhs),
                                         ROOTED, this->rowIndexBase_);
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
Superlu<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  using Teuchos::RCP;
  using Teuchos::getIntegralValue;
  using Teuchos::ParameterEntryValidator;

  RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

  data_.options.Trans = this->control_.useTranspose_ ? SLU::TRANS : SLU::NOTRANS;
  // The SuperLU transpose option can override the Amesos2 option
  if( parameterList->isParameter("Trans") ){
    RCP<const ParameterEntryValidator> trans_validator = valid_params->getEntry("Trans").validator();
    parameterList->getEntry("Trans").setValidator(trans_validator);

    data_.options.Trans = getIntegralValue<SLU::trans_t>(*parameterList, "Trans");
  }

  if( parameterList->isParameter("IterRefine") ){
    RCP<const ParameterEntryValidator> refine_validator = valid_params->getEntry("IterRefine").validator();
    parameterList->getEntry("IterRefine").setValidator(refine_validator);

    data_.options.IterRefine = getIntegralValue<SLU::IterRefine_t>(*parameterList, "IterRefine");
  }

  if( parameterList->isParameter("ColPerm") ){
    RCP<const ParameterEntryValidator> colperm_validator = valid_params->getEntry("ColPerm").validator();
    parameterList->getEntry("ColPerm").setValidator(colperm_validator);

    data_.options.ColPerm = getIntegralValue<SLU::colperm_t>(*parameterList, "ColPerm");
  }

  data_.options.DiagPivotThresh = parameterList->get<double>("DiagPivotThresh", 1.0);

  bool equil = parameterList->get<bool>("Equil", true);
  data_.options.Equil = equil ? SLU::YES : SLU::NO;

  bool symmetric_mode = parameterList->get<bool>("SymmetricMode", false);
  data_.options.SymmetricMode = symmetric_mode ? SLU::YES : SLU::NO;
}


template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
Superlu<Matrix,Vector>::getValidParameters_impl() const
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

    setStringToIntegralParameter<SLU::trans_t>("Trans", "NOTRANS",
                                               "Solve for the transpose system or not",
                                               tuple<string>("TRANS","NOTRANS","CONJ"),
                                               tuple<string>("Solve with transpose",
                                                             "Do not solve with transpose",
                                                             "Solve with the conjugate transpose"),
                                               tuple<SLU::trans_t>(SLU::TRANS,
                                                                   SLU::NOTRANS,
                                                                   SLU::CONJ),
                                               pl.getRawPtr());

    setStringToIntegralParameter<SLU::IterRefine_t>("IterRefine", "NOREFINE",
                                                    "Type of iterative refinement to use",
                                                    tuple<string>("NOREFINE", "SLU_SINGLE", "SLU_DOUBLE"),
                                                    tuple<string>("Do not use iterative refinement",
                                                                  "Do single iterative refinement",
                                                                  "Do double iterative refinement"),
                                                    tuple<SLU::IterRefine_t>(SLU::NOREFINE,
                                                                             SLU::SLU_SINGLE,
                                                                             SLU::SLU_DOUBLE),
                                                    pl.getRawPtr());

    // Note: MY_PERMC not yet supported
    setStringToIntegralParameter<SLU::colperm_t>("ColPerm", "COLAMD",
                                                 "Specifies how to permute the columns of the "
                                                 "matrix for sparsity preservation",
                                                 tuple<string>("NATURAL", "MMD_AT_PLUS_A",
                                                               "MMD_ATA", "COLAMD"),
                                                 tuple<string>("Natural ordering",
                                                               "Minimum degree ordering on A^T + A",
                                                               "Minimum degree ordering on A^T A",
                                                               "Approximate minimum degree column ordering"),
                                                 tuple<SLU::colperm_t>(SLU::NATURAL,
                                                                       SLU::MMD_AT_PLUS_A,
                                                                       SLU::MMD_ATA,
                                                                       SLU::COLAMD),
                                                 pl.getRawPtr());

    Teuchos::RCP<EnhancedNumberValidator<double> > diag_pivot_thresh_validator
      = Teuchos::rcp( new EnhancedNumberValidator<double>(0.0, 1.0) );
    pl->set("DiagPivotThresh", 1.0,
            "Specifies the threshold used for a diagonal entry to be an acceptable pivot",
            diag_pivot_thresh_validator); // partial pivoting

    pl->set("Equil", true, "Whether to equilibrate the system before solve");

    pl->set("SymmetricMode", false,
            "Specifies whether to use the symmetric mode. "
            "Gives preference to diagonal pivots and uses "
            "an (A^T + A)-based column permutation.");

    valid_params = pl;
  }

  return valid_params;
}


template <class Matrix, class Vector>
bool
Superlu<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
  using Teuchos::as;

#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

  // SuperLU does not need the matrix at symbolicFactorization()
  if( current_phase == SYMBFACT ) return false;

  // Cleanup old store memory if it's non-NULL (should only ever be non-NULL at root_)
  if( data_.A.Store != NULL ){
    SLU::Destroy_SuperMatrix_Store( &(data_.A) );
    data_.A.Store = NULL;
  }

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
    Util::get_ccs_helper<
    MatrixAdapter<Matrix>,slu_type,int,int>::do_get(this->matrixA_.ptr(),
                                                    nzvals_(), rowind_(),
                                                    colptr_(), nnz_ret, ROOTED,
                                                    ARBITRARY,
                                                    this->rowIndexBase_);
  }

  // Get the SLU data type for this type of matrix
  SLU::Dtype_t dtype = type_map::dtype;

  if( this->root_ ){
    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<int>(this->globalNumNonZeros_),
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");

    function_map::create_CompCol_Matrix( &(data_.A),
                                         this->globalNumRows_, this->globalNumCols_,
                                         nnz_ret,
                                         nzvals_.getRawPtr(),
                                         rowind_.getRawPtr(),
                                         colptr_.getRawPtr(),
                                         SLU::SLU_NC, dtype, SLU::SLU_GE);
  }

  return true;
}


template<class Matrix, class Vector>
const char* Superlu<Matrix,Vector>::name = "SuperLU";
  

} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_DEF_HPP
