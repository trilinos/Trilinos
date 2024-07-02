// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_Superlumt_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Tue May 24 08:37:17 MDT 2011

   \brief  Definitions for the Amesos2 SuperLU_MT solver interface
*/

#ifndef AMESOS2_SUPERLUMT_DEF_HPP
#define AMESOS2_SUPERLUMT_DEF_HPP

#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

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
  Superlumt<Matrix,Vector>::Superlumt(Teuchos::RCP<const Matrix> A,
                                      Teuchos::RCP<Vector>       X,
                                      Teuchos::RCP<const Vector> B )
    : SolverCore<Amesos2::Superlumt,Matrix,Vector>(A, X, B)
    , nzvals_()                 // initialize to empty arrays
    , rowind_()
    , colptr_()
    , is_contiguous_(true)
  {
    Teuchos::RCP<Teuchos::ParameterList> default_params
      = Teuchos::parameterList( *(this->getValidParameters()) );
    this->setParameters(default_params);

    data_.options.lwork = 0;    // Use system memory for factorization

    data_.perm_r.resize(this->globalNumRows_);
    data_.perm_c.resize(this->globalNumCols_);
    data_.options.perm_r = data_.perm_r.getRawPtr();
    data_.options.perm_c = data_.perm_c.getRawPtr();

    data_.R.resize(this->globalNumRows_);
    data_.C.resize(this->globalNumCols_);

    data_.options.refact = SLUMT::NO; // initially we are not refactoring
    data_.equed = SLUMT::NOEQUIL;       // No equilibration has yet been performed
    data_.rowequ = false;
    data_.colequ = false;

    data_.A.Store  = NULL;
    data_.AC.Store = NULL;
    data_.BX.Store = NULL;
    data_.L.Store  = NULL;
    data_.U.Store  = NULL;

    data_.stat.ops = NULL;
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

    if ( data_.L.Store != NULL ){ // will only ever be true for this->root_
      SLUMT::D::Destroy_SuperNode_SCP( &(data_.L) );
      SLUMT::D::Destroy_CompCol_NCP( &(data_.U) );

      // memory alloc'd in sp_colorder
      free( data_.options.etree );
      free( data_.options.colcnt_h );
      free( data_.options.part_super_h );
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
    
    if ( data_.stat.ops != NULL )
      SLUMT::D::StatFree( &(data_.stat) );
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
      data_.L.Store = NULL;
      data_.U.Store = NULL;
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
    TEUCHOS_TEST_FOR_EXCEPTION( nprocs <= 0,
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
        TEUCHOS_TEST_FOR_EXCEPTION( info != 0,
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
        if( data_.options.refact == SLUMT::NO ){                 // then we recompute etree; free the old one
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
      const int n = as<int>(this->globalNumCols_); // n is the number of columns in A
      if( data_.stat.ops != NULL ){ SLUMT::D::StatFree( &(data_.stat) ); data_.stat.ops = NULL; }
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
      this->setNnzLU( as<size_t>(((SLUMT::SCformat*)data_.L.Store)->nnz +
                                 ((SLUMT::NCformat*)data_.U.Store)->nnz) );
    }

    // Check output
    const global_size_type info_st = as<global_size_type>(info);
    TEUCHOS_TEST_FOR_EXCEPTION( (info_st > 0) && (info_st <= this->globalNumCols_),
                        std::runtime_error,
                        "Factorization complete, but matrix is singular. Division by zero eminent");
    TEUCHOS_TEST_FOR_EXCEPTION( (info_st > 0) && (info_st > this->globalNumCols_),
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

    {                           // Get values from B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor convTimer( this->timers_.vecConvTime_ );
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      if ( is_contiguous_ == true ) {
        Util::get_1d_copy_helper<MultiVecAdapter<Vector>,
          slu_type>::do_get(B, bxvals_(),
              ldbx_,
              ROOTED, this->rowIndexBase_);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "NON_CONTIGUOUS_AND_ROOTED not supported.");
      }
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
        if( data_.rowequ ){             // row equilibration has been done on AC
          // scale bxvals_ by diag(R)
          Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.R(),
                      SLUMT::slu_mt_mult<slu_type,magnitude_type>());
        }
      } else if( data_.colequ ){        // column equilibration has been done on AC
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

    TEUCHOS_TEST_FOR_EXCEPTION( info < 0,
                        std::runtime_error,
                        "Argument " << -info << " to gstrs had an illegal value" );

    // "Un-scale" the solution so that it is a solution of the original system
    if( data_.options.trans == SLUMT::NOTRANS ){
      if( data_.colequ ){       // column equilibration has been done on AC
        // scale bxvals_ by diag(C)
        Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.C(),
                    SLUMT::slu_mt_mult<slu_type,magnitude_type>());
      }
    } else if( data_.rowequ ){          // row equilibration has been done on AC
      // scale bxvals_ by diag(R)
      Util::scale(bxvals_(), as<size_t>(ld_rhs), ldbx_, data_.R(),
                  SLUMT::slu_mt_mult<slu_type,magnitude_type>());
    }

    /* Update X's global values */
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      if ( is_contiguous_ == true ) {
        Util::put_1d_data_helper<
          MultiVecAdapter<Vector>, slu_type>::do_put(X, bxvals_(), ldbx_, ROOTED);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "NON_CONTIGUOUS_AND_ROOTED not supported.");
      }
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
  Superlumt<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterEntryValidator;

    RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

    int default_nprocs = Kokkos::DefaultHostExecutionSpace().concurrency();
    data_.options.nprocs = parameterList->get<int>("nprocs", default_nprocs);

    data_.options.trans = this->control_.useTranspose_ ? SLUMT::TRANS : SLUMT::NOTRANS;
    // SuperLU_MT "trans" option can override the Amesos2 option
    if( parameterList->isParameter("trans") ){
      RCP<const ParameterEntryValidator> trans_validator = valid_params->getEntry("trans").validator();
      parameterList->getEntry("trans").setValidator(trans_validator);

      data_.options.trans = getIntegralValue<SLUMT::trans_t>(*parameterList, "trans");
    }

    data_.options.panel_size = parameterList->get<int>("panel_size", SLUMT::D::sp_ienv(1));

    data_.options.relax = parameterList->get<int>("relax", SLUMT::D::sp_ienv(2));

    const bool equil = parameterList->get<bool>("Equil", true);
    data_.options.fact = equil ? SLUMT::EQUILIBRATE : SLUMT::DOFACT;

    const bool symmetric_mode = parameterList->get<bool>("SymmetricMode", false);
    data_.options.SymmetricMode = symmetric_mode ? SLUMT::YES : SLUMT::NO;

    const bool print_stat = parameterList->get<bool>("PrintStat", false);
    data_.options.PrintStat = print_stat ? SLUMT::YES : SLUMT::NO;

    data_.options.diag_pivot_thresh = parameterList->get<double>("diag_pivot_thresh", 1.0);

    if( parameterList->isParameter("ColPerm") ){
      RCP<const ParameterEntryValidator> colperm_validator = valid_params->getEntry("ColPerm").validator();
      parameterList->getEntry("ColPerm").setValidator(colperm_validator);

      data_.options.ColPerm = getIntegralValue<SLUMT::colperm_t>(*parameterList, "ColPerm");
    }

    // TODO: until we support retrieving col/row permutations from the user
    data_.options.usepr = SLUMT::NO;

    if( parameterList->isParameter("IsContiguous") ){
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }
  }


  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  Superlumt<Matrix,Vector>::getValidParameters_impl() const
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

      int default_nprocs = Kokkos::DefaultHostExecutionSpace().concurrency();
      Teuchos::RCP<EnhancedNumberValidator<int> > nprocs_validator
        = Teuchos::rcp( new EnhancedNumberValidator<int>() );
      nprocs_validator->setMin(default_nprocs);
      pl->set("nprocs", default_nprocs, "The number of processors to factorize with", nprocs_validator);

      setStringToIntegralParameter<SLUMT::trans_t>("trans", "NOTRANS",
                                                   "Solve for the transpose system or not",
                                                   tuple<string>("TRANS","NOTRANS","CONJ"),
                                                   tuple<string>("Solve with transpose",
                                                                 "Do not solve with transpose",
                                                                 "Solve with the conjugate transpose"),
                                                   tuple<SLUMT::trans_t>(SLUMT::TRANS,
                                                                         SLUMT::NOTRANS,
                                                                         SLUMT::CONJ),
                                                   pl.getRawPtr());

      Teuchos::RCP<EnhancedNumberValidator<int> > panel_size_validator
        = Teuchos::rcp( new EnhancedNumberValidator<int>() );
      panel_size_validator->setMin(1);
      pl->set("panel_size", SLUMT::D::sp_ienv(1),
              "Specifies the number of consecutive columns to be treated as a unit of task.",
              panel_size_validator);

      Teuchos::RCP<EnhancedNumberValidator<int> > relax_validator
        = Teuchos::rcp( new EnhancedNumberValidator<int>() );
      relax_validator->setMin(1);
      pl->set("relax", SLUMT::D::sp_ienv(2),
              "Specifies the number of columns to be grouped as a relaxed supernode",
              relax_validator);

      pl->set("Equil", true, "Whether to equilibrate the system before solve");

      Teuchos::RCP<EnhancedNumberValidator<double> > diag_pivot_thresh_validator
        = Teuchos::rcp( new EnhancedNumberValidator<double>(0.0, 1.0) );
      pl->set("diag_pivot_thresh", 1.0,
              "Specifies the threshold used for a diagonal entry to be an acceptable pivot",
              diag_pivot_thresh_validator); // partial pivoting

      // Note: MY_PERMC not yet supported
      setStringToIntegralParameter<SLUMT::colperm_t>("ColPerm", "COLAMD",
                                                     "Specifies how to permute the columns of the "
                                                     "matrix for sparsity preservation",
                                                     tuple<string>("NATURAL",
                                                                   "MMD_AT_PLUS_A",
                                                                   "MMD_ATA",
                                                                   "COLAMD"),
                                                     tuple<string>("Natural ordering",
                                                                       "Minimum degree ordering on A^T + A",
                                                                   "Minimum degree ordering on A^T A",
                                                                   "Approximate minimum degree column ordering"),
                                                     tuple<SLUMT::colperm_t>(SLUMT::NATURAL,
                                                                             SLUMT::MMD_AT_PLUS_A,
                                                                             SLUMT::MMD_ATA,
                                                                             SLUMT::COLAMD),
                                                     pl.getRawPtr());

      pl->set("SymmetricMode", false, "Specifies whether to use the symmetric mode");

      // TODO: until we support getting row/col permutations from user
      //      pl->set("usepr", false);

      pl->set("PrintStat", false, "Specifies whether to print the solver's statistics");

      pl->set("IsContiguous", true, "Whether GIDs contiguous");

      valid_params = pl;
    }

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
      if ( is_contiguous_ == true ) {
        Util::get_ccs_helper<
          MatrixAdapter<Matrix>,slu_type,int,int>::do_get(this->matrixA_.ptr(),
              nzvals_, rowind_, colptr_,
              nnz_ret, ROOTED, ARBITRARY, this->rowIndexBase_);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
          "NON_CONTIGUOUS_AND_ROOTED not supported.");
      }
    }

    if( this->root_ ){
      TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<int>(this->globalNumNonZeros_),
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

      TEUCHOS_TEST_FOR_EXCEPTION( data_.A.Store == NULL,
                          std::runtime_error,
                          "Failed to allocate matrix store" );
    }

    return true;
  }


  template<class Matrix, class Vector>
  const char* Superlumt<Matrix,Vector>::name = "SuperLU_MT";


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLUMT_DEF_HPP
