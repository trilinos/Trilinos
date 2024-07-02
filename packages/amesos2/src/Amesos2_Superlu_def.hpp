// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
// TODO: This 'using namespace SLU' is not a good thing.
// It was added because kernels does not namespace SuperLU but Amesos2 does.
// Declaring the namespace SLU allows us to reuse that file without duplication.
// We need to determine a uniform system to avoid this this but for now, at least
// this issue is only present when TRSV is enabled.
using namespace SLU;
#include "KokkosSparse_sptrsv_superlu.hpp"
#endif

namespace Amesos2 {


template <class Matrix, class Vector>
Superlu<Matrix,Vector>::Superlu(
  Teuchos::RCP<const Matrix> A,
  Teuchos::RCP<Vector>       X,
  Teuchos::RCP<const Vector> B )
  : SolverCore<Amesos2::Superlu,Matrix,Vector>(A, X, B)
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  , sptrsv_invert_diag_(true)
  , sptrsv_invert_offdiag_ (false)
  , sptrsv_u_in_csr_ (true)
  , sptrsv_merge_supernodes_ (false)
  , sptrsv_use_spmv_ (false)
#endif
  , is_contiguous_(true) // default is set by params
  , use_triangular_solves_(false) // default is set by params
  , use_metis_(false)
  , symmetrize_metis_(true)
{
  // ilu_set_default_options is called later in set parameter list if required.
  // This is not the ideal way, but the other option to always call
  // ilu_set_default_options here and assuming it won't have any side effect
  // in the TPL is more dangerous. It is not a good idea to rely on external
  // libraries' internal "features".
  SLU::set_default_options(&(data_.options));
  // Override some default options
  data_.options.PrintStat = SLU::NO;

  SLU::StatInit(&(data_.stat));

  Kokkos::resize(data_.perm_r, this->globalNumRows_);
  Kokkos::resize(data_.perm_c, this->globalNumCols_);
  Kokkos::resize(data_.etree, this->globalNumCols_);
  Kokkos::resize(data_.R, this->globalNumRows_);
  Kokkos::resize(data_.C, this->globalNumCols_);

  data_.relax = SLU::sp_ienv(2); // Query optimal relax param from superlu
  data_.panel_size = SLU::sp_ienv(1); // Query optimal panel size

  data_.equed = 'N';            // No equilibration
  data_.rowequ = false;
  data_.colequ = false;
  data_.A.Store = NULL;
  data_.L.Store = NULL;
  data_.U.Store = NULL;
  data_.X.Store = NULL;
  data_.B.Store = NULL;

  ILU_Flag_=false; // default: turn off ILU
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
  if ( data_.L.Store != NULL ){ // will only be true for this->root_
    SLU::Destroy_SuperNode_Matrix( &(data_.L) );
    SLU::Destroy_CompCol_Matrix( &(data_.U) );
  }
}

template <class Matrix, class Vector>
std::string
Superlu<Matrix,Vector>::description() const
{
  std::ostringstream oss;
  oss << "SuperLU solver interface";
  if (ILU_Flag_) {
    oss << ", \"ILUTP\" : {";
    oss << "drop tol = " << data_.options.ILU_DropTol;
    oss << ", fill factor = " << data_.options.ILU_FillFactor;
    oss << ", fill tol = " << data_.options.ILU_FillTol;
    switch(data_.options.ILU_MILU) {
      case SLU::SMILU_1 :
         oss << ", MILU 1";
         break;
      case SLU::SMILU_2  :
         oss << ", MILU 2";
         break;
      case SLU::SMILU_3  :
         oss << ", MILU 3";
         break;
      case SLU::SILU     :
      default:
         oss << ", regular ILU";
    }
    switch(data_.options.ILU_Norm) {
      case SLU::ONE_NORM :
         oss << ", 1-norm";
         break;
      case SLU::TWO_NORM  :
         oss << ", 2-norm";
         break;
      case SLU::INF_NORM  :
      default:
         oss << ", infinity-norm";
    }
    oss << "}";
  } else {
    oss << ", direct solve";
  }
  return oss.str();
  /*

  // ILU parameters
  if( parameterList->isParameter("RowPerm") ){
    RCP<const ParameterEntryValidator> rowperm_validator = valid_params->getEntry("RowPerm").validator();
    parameterList->getEntry("RowPerm").setValidator(rowperm_validator);
    data_.options.RowPerm = getIntegralValue<SLU::rowperm_t>(*parameterList, "RowPerm");
  }


  */
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
  if (!use_metis_ && permc_spec != SLU::MY_PERMC && this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor preOrderTimer(this->timers_.preOrderTime_);
#endif

    SLU::get_perm_c(permc_spec, &(data_.A), data_.perm_c.data());
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
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor symFactTime( this->timers_.symFactTime_ );
#endif
  same_symbolic_ = false;
  data_.options.Fact = SLU::DOFACT;

  if (use_metis_) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::RCP< Teuchos::Time > MetisTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for METIS");
    Teuchos::TimeMonitor numFactTimer(*MetisTimer_);
#endif
    size_t n  = this->globalNumRows_;
    size_t nz = this->globalNumNonZeros_;
    host_size_type_array    host_rowind_view_ = host_rows_view_;
    host_ordinal_type_array host_colptr_view_ = host_col_ptr_view_;

    /* symmetrize the matrix (A + A') if needed */
    if (symmetrize_metis_) {
      int new_nz = 0;
      int *new_col_ptr;
      int *new_row_ind;

      // NOTE: both size_type and ordinal_type are defined as int 
      //       Consider using "symmetrize_graph" in KokkosKernels, if this becomes significant in time.
      SLU::at_plus_a(n, nz, host_col_ptr_view_.data(), host_rows_view_.data(),
                     &new_nz, &new_col_ptr, &new_row_ind);

      // reallocate (so that we won't overwrite the original)
      // and copy for output
      nz = new_nz;
      Kokkos::resize (host_rowind_view_, nz);
      Kokkos::realloc(host_colptr_view_, n+1);
      for (size_t i = 0; i <= n; i++) {
        host_colptr_view_(i) = new_col_ptr[i];
      }
      for (size_t i = 0; i < nz; i++) {
        host_rowind_view_(i) = new_row_ind[i];
      }

      // free
      SLU::SUPERLU_FREE(new_col_ptr);
      if (new_nz > 0) {
        SLU::SUPERLU_FREE(new_row_ind);
      }
    }

    // reorder will convert both graph and perm/iperm to the internal METIS integer type
    using metis_int_array_type = host_ordinal_type_array;
    metis_int_array_type metis_perm  = metis_int_array_type(Kokkos::ViewAllocateWithoutInitializing("metis_perm"),  n);
    metis_int_array_type metis_iperm = metis_int_array_type(Kokkos::ViewAllocateWithoutInitializing("metis_iperm"), n);

    // call METIS
    size_t new_nnz = 0;
    Amesos2::Util::reorder(
      host_colptr_view_, host_rowind_view_, 
      metis_perm, metis_iperm, new_nnz, false);

    for (size_t i = 0; i < n; i++) {
      data_.perm_r(i) = metis_iperm(i);
      data_.perm_c(i) = metis_iperm(i);
    }

    // SLU will use user-specified METIS ordering
    data_.options.ColPerm = SLU::MY_PERMC;
    data_.options.RowPerm = SLU::MY_PERMR;
    // Turn off Equil not to mess up METIS ordering?
    //data_.options.Equil = SLU::NO;
  }
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
      function_map::gsequ(&(data_.A), data_.R.data(),
                          data_.C.data(), &rowcnd, &colcnd,
                          &amax, &info2);
      TEUCHOS_TEST_FOR_EXCEPTION
        (info2 < 0, std::runtime_error,
         "SuperLU's gsequ function returned with status " << info2 << " < 0."
         "  This means that argument " << (-info2) << " given to the function"
         " had an illegal value.");
      if (info2 > 0) {
        if (info2 <= data_.A.nrow) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, "SuperLU's gsequ function returned with "
             "info = " << info2 << " > 0, and info <= A.nrow = " << data_.A.nrow
             << ".  This means that row " << info2 << " of A is exactly zero.");
        }
        else if (info2 > data_.A.ncol) {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, "SuperLU's gsequ function returned with "
             "info = " << info2 << " > 0, and info > A.ncol = " << data_.A.ncol
             << ".  This means that column " << (info2 - data_.A.nrow) << " of "
             "A is exactly zero.");
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION
            (true, std::runtime_error, "SuperLU's gsequ function returned "
             "with info = " << info2 << " > 0, but its value is not in the "
             "range permitted by the documentation.  This should never happen "
             "(it appears to be SuperLU's fault).");
        }
      }

      // apply row and column scalings if necessary
      function_map::laqgs(&(data_.A), data_.R.data(),
                          data_.C.data(), rowcnd, colcnd,
                          amax, &(data_.equed));

      // check what types of equilibration was actually done
      data_.rowequ = (data_.equed == 'R') || (data_.equed == 'B');
      data_.colequ = (data_.equed == 'C') || (data_.equed == 'B');
    }

    // Apply the column permutation computed in preOrdering.  Place the
    // column-permuted matrix in AC
    SLU::sp_preorder(&(data_.options), &(data_.A), data_.perm_c.data(),
                     data_.etree.data(), &(data_.AC));

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

      if(ILU_Flag_==false) {
        function_map::gstrf(&(data_.options), &(data_.AC),
            data_.relax, data_.panel_size, data_.etree.data(),
            NULL, 0, data_.perm_c.data(), data_.perm_r.data(),
            &(data_.L), &(data_.U), 
#ifdef HAVE_AMESOS2_SUPERLU5_API
            &(data_.lu), 
#endif
            &(data_.stat), &info);
      }
      else {
        function_map::gsitrf(&(data_.options), &(data_.AC),
            data_.relax, data_.panel_size, data_.etree.data(),
            NULL, 0, data_.perm_c.data(), data_.perm_r.data(),
            &(data_.L), &(data_.U), 
#ifdef HAVE_AMESOS2_SUPERLU5_API
            &(data_.lu), 
#endif
            &(data_.stat), &info);
      }

      if (data_.options.ConditionNumber == SLU::YES) {
        char norm[1];
        if (data_.options.Trans == SLU::NOTRANS) {
            *(unsigned char *)norm = '1';
        } else {
            *(unsigned char *)norm = 'I';
        }

        data_.anorm = function_map::langs(norm, &(data_.A));
        function_map::gscon(norm, &(data_.L), &(data_.U),
                            data_.anorm, &(data_.rcond),
                            &(data_.stat), &info);
      }
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

  if(use_triangular_solves_) {
#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
    #if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
    if (this->getComm()->getRank()) {
      std::cout << " > Metis           : " << (use_metis_ ? "YES" : "NO") << std::endl;
      std::cout << " > Equil           : " << (data_.options.Equil == SLU::YES ? "YES" : "NO") << std::endl;
      std::cout << " > Cond Number     : " << (data_.options.ConditionNumber == SLU::YES ? "YES" : "NO") << std::endl;
      std::cout << " > Invert diag     : " << sptrsv_invert_diag_ << std::endl;
      std::cout << " > Invert off-diag : " << sptrsv_invert_offdiag_ << std::endl;
      std::cout << " > U in CSR        : " << sptrsv_u_in_csr_ << std::endl;
      std::cout << " > Merge           : " << sptrsv_merge_supernodes_ << std::endl;
      std::cout << " > Use SpMV        : " << sptrsv_use_spmv_ << std::endl;
    }
    //std::cout << myRank << " : siize(A) " << (data_.A.nrow) << " x " << (data_.A.ncol) << std::endl;
    //std::cout << myRank << " : nnz(A)   " << ((SLU::NCformat*)data_.A.Store)->nnz << std::endl;
    //std::cout << myRank << " : nnz(L)   " << ((SLU::SCformat*)data_.L.Store)->nnz << std::endl;
    //std::cout << myRank << " : nnz(U)   " << ((SLU::SCformat*)data_.U.Store)->nnz << std::endl;
    #endif
#endif
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::RCP< Teuchos::Time > SpTrsvTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for SpTrsv setup");
    Teuchos::TimeMonitor numFactTimer(*SpTrsvTimer_);
#endif
    triangular_solve_factor();
  }

  return(info);
}


template <class Matrix, class Vector>
int
Superlu<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                   const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
{
  using Teuchos::as;
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::RCP< Teuchos::Time > Amesos2SolveTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for Amesos2");
    Teuchos::TimeMonitor solveTimer(*Amesos2SolveTimer_);
#endif

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  bool bDidAssignX = false; // will be set below
  bool bDidAssignB = false; // will be set below
  {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
    Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

    // In general we may want to write directly to the x space without a copy.
    // So we 'get' x which may be a direct view assignment to the MV.
    const bool initialize_data = true;
    const bool do_not_initialize_data = false;
    if(use_triangular_solves_) { // to device
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
      bDidAssignX = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        device_solve_array_t>::do_get(do_not_initialize_data, X, device_xValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
      bDidAssignB = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        device_solve_array_t>::do_get(initialize_data, B, device_bValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
#endif
    }
    else { // to host
      bDidAssignX = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        host_solve_array_t>::do_get(do_not_initialize_data, X, host_xValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
      bDidAssignB = Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
        host_solve_array_t>::do_get(initialize_data, B, host_bValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
    }
  }

  // If equilibration was applied at numeric, then gssvx and gsisx are going to
  // modify B, so we can't use the optimized assignment to B since we'll change
  // the source test vector and then fail the subsequent cycle. We need a copy.
  // If bDidAssignB is false, then we skip the copy since it was copied already.
  if(bDidAssignB && data_.equed != 'N') {
    if(use_triangular_solves_) { // to device
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
      device_solve_array_t copyB(Kokkos::ViewAllocateWithoutInitializing("copyB"),
        device_bValues_.extent(0), device_bValues_.extent(1));
      Kokkos::deep_copy(copyB, device_bValues_);
      device_bValues_ = copyB;
#endif
    }
    else {
      host_solve_array_t copyB(Kokkos::ViewAllocateWithoutInitializing("copyB"),
        host_bValues_.extent(0), host_bValues_.extent(1));
      Kokkos::deep_copy(copyB, host_bValues_);
      host_bValues_ = copyB;
    }
  }

  int ierr = 0; // returned error code

  magnitude_type rpg, rcond;
  if ( this->root_ ) {
    // Note: the values of B and X (after solution) are adjusted
    // appropriately within gssvx for row and column scalings.
    {                           // Do solve!
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

    if(use_triangular_solves_) {
      triangular_solve();
    }
    else {
      Kokkos::resize(data_.ferr, nrhs);
      Kokkos::resize(data_.berr, nrhs);

      {
  #ifdef HAVE_AMESOS2_TIMERS
        Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
  #endif
        SLU::Dtype_t dtype = type_map::dtype;
        int i_ld_rhs = as<int>(ld_rhs);
        function_map::create_Dense_Matrix(&(data_.B), i_ld_rhs, as<int>(nrhs),
                                          convert_bValues_, host_bValues_, i_ld_rhs,
                                          SLU::SLU_DN, dtype, SLU::SLU_GE);
        function_map::create_Dense_Matrix(&(data_.X), i_ld_rhs, as<int>(nrhs),
                                          convert_xValues_, host_xValues_, i_ld_rhs,
                                          SLU::SLU_DN, dtype, SLU::SLU_GE);
      }

      if(ILU_Flag_==false) {
        function_map::gssvx(&(data_.options), &(data_.A),
            data_.perm_c.data(), data_.perm_r.data(),
            data_.etree.data(), &(data_.equed), data_.R.data(),
            data_.C.data(), &(data_.L), &(data_.U), NULL, 0, &(data_.B),
            &(data_.X), &rpg, &rcond, data_.ferr.data(),
            data_.berr.data(),
  #ifdef HAVE_AMESOS2_SUPERLU5_API
            &(data_.lu),
  #endif
            &(data_.mem_usage), &(data_.stat), &ierr);
      }
      else {
        function_map::gsisx(&(data_.options), &(data_.A),
            data_.perm_c.data(), data_.perm_r.data(),
            data_.etree.data(), &(data_.equed), data_.R.data(),
            data_.C.data(), &(data_.L), &(data_.U), NULL, 0, &(data_.B),
            &(data_.X), &rpg, &rcond,
  #ifdef HAVE_AMESOS2_SUPERLU5_API
            &(data_.lu),
  #endif
            &(data_.mem_usage), &(data_.stat), &ierr);
      }

      // Cleanup X and B stores
      SLU::Destroy_SuperMatrix_Store( &(data_.X) );
      SLU::Destroy_SuperMatrix_Store( &(data_.B) );
      data_.X.Store = NULL;
      data_.B.Store = NULL;
    }

    } // do solve

    // convert back to Kokkos views
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer(this->timers_.vecConvTime_);
#endif
      function_map::convert_back_Dense_Matrix(convert_bValues_, host_bValues_);
      function_map::convert_back_Dense_Matrix(convert_xValues_, host_xValues_);
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

  // if bDidAssignX, then we solved straight to the adapter's X memory space without
  // requiring additional memory allocation, so the x data is already in place.
  if(!bDidAssignX) {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

    if(use_triangular_solves_) { // to device
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
      Util::put_1d_data_helper_kokkos_view<
        MultiVecAdapter<Vector>,device_solve_array_t>::do_put(X, device_xValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
#endif
    }
    else {
      Util::put_1d_data_helper_kokkos_view<
        MultiVecAdapter<Vector>,host_solve_array_t>::do_put(X, host_xValues_,
            as<size_t>(ld_rhs),
            (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
            this->rowIndexBase_);
    }
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

  ILU_Flag_ = parameterList->get<bool>("ILU_Flag",false);
  if (ILU_Flag_) {
      SLU::ilu_set_default_options(&(data_.options));
      // Override some default options
      data_.options.PrintStat = SLU::NO;
  }

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

  bool condNum = parameterList->get<bool>("ConditionNumber", false);
  data_.options.ConditionNumber = condNum ? SLU::YES : SLU::NO;

  bool symmetric_mode = parameterList->get<bool>("SymmetricMode", false);
  data_.options.SymmetricMode = symmetric_mode ? SLU::YES : SLU::NO;

  // ILU parameters
  if( parameterList->isParameter("RowPerm") ){
    RCP<const ParameterEntryValidator> rowperm_validator = valid_params->getEntry("RowPerm").validator();
    parameterList->getEntry("RowPerm").setValidator(rowperm_validator);
    data_.options.RowPerm = getIntegralValue<SLU::rowperm_t>(*parameterList, "RowPerm");
  }

  /*if( parameterList->isParameter("ILU_DropRule") ) {
    RCP<const ParameterEntryValidator> droprule_validator = valid_params->getEntry("ILU_DropRule").validator();
    parameterList->getEntry("ILU_DropRule").setValidator(droprule_validator);
    data_.options.ILU_DropRule = getIntegralValue<SLU::rule_t>(*parameterList, "ILU_DropRule");
  }*/

  data_.options.ILU_DropTol = parameterList->get<double>("ILU_DropTol", 0.0001);

  data_.options.ILU_FillFactor = parameterList->get<double>("ILU_FillFactor", 10.0);

  if( parameterList->isParameter("ILU_Norm") ) {
    RCP<const ParameterEntryValidator> norm_validator = valid_params->getEntry("ILU_Norm").validator();
    parameterList->getEntry("ILU_Norm").setValidator(norm_validator);
    data_.options.ILU_Norm = getIntegralValue<SLU::norm_t>(*parameterList, "ILU_Norm");
  }

  if( parameterList->isParameter("ILU_MILU") ) {
    RCP<const ParameterEntryValidator> milu_validator = valid_params->getEntry("ILU_MILU").validator();
    parameterList->getEntry("ILU_MILU").setValidator(milu_validator);
    data_.options.ILU_MILU = getIntegralValue<SLU::milu_t>(*parameterList, "ILU_MILU");
  }

  data_.options.ILU_FillTol = parameterList->get<double>("ILU_FillTol", 0.01);

  is_contiguous_ = parameterList->get<bool>("IsContiguous", true);

  use_metis_ = parameterList->get<bool>("UseMetis", false);
  symmetrize_metis_ = parameterList->get<bool>("SymmetrizeMetis", true);

  use_triangular_solves_ = parameterList->get<bool>("Enable_KokkosKernels_TriangularSolves", false);
  if(use_triangular_solves_) {
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
    // specify whether to invert diagonal blocks
    sptrsv_invert_diag_ = parameterList->get<bool>("SpTRSV_Invert_Diag", true);
    // specify whether to apply diagonal-inversion to off-diagonal blocks (optional, default is false)
    sptrsv_invert_offdiag_ = parameterList->get<bool>("SpTRSV_Invert_OffDiag", false);
    // specify whether to store U in CSR format
    sptrsv_u_in_csr_ = parameterList->get<bool>("SpTRSV_U_CSR", true);
    // specify whether to merge supernodes with similar sparsity structures
    sptrsv_merge_supernodes_ = parameterList->get<bool>("SpTRSV_Merge_Supernodes", false);
    // specify whether to use spmv for sptrsv
    sptrsv_use_spmv_ = parameterList->get<bool>("SpTRSV_Use_SpMV", false);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
      "Calling for triangular solves but KokkosKernels_ENABLE_SUPERNODAL_SPTRSV and KokkosKernels_ENABLE_TPL_SUPERLU were not configured." );
#endif
  }
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
    pl->set("ConditionNumber", false, "Whether to approximate condition number");

    pl->set("SymmetricMode", false,
            "Specifies whether to use the symmetric mode. "
            "Gives preference to diagonal pivots and uses "
            "an (A^T + A)-based column permutation.");

    // ILU parameters

    setStringToIntegralParameter<SLU::rowperm_t>("RowPerm", "LargeDiag",
            "Type of row permutation strategy to use",
            tuple<string>("NOROWPERM","LargeDiag","MY_PERMR"),
            tuple<string>("Use natural ordering",
            "Use weighted bipartite matching algorithm",
            "Use the ordering given in perm_r input"),
            tuple<SLU::rowperm_t>(SLU::NOROWPERM,
#if SUPERLU_MAJOR_VERSION > 5 || (SUPERLU_MAJOR_VERSION == 5 && SUPERLU_MINOR_VERSION > 2) || (SUPERLU_MAJOR_VERSION == 5 && SUPERLU_MINOR_VERSION == 2 && SUPERLU_PATCH_VERSION >= 2)
            SLU::LargeDiag_MC64,
#else
            SLU::LargeDiag,
#endif
            SLU::MY_PERMR),
            pl.getRawPtr());

    /*setStringToIntegralParameter<SLU::rule_t>("ILU_DropRule", "DROP_BASIC",
            "Type of dropping strategy to use",
            tuple<string>("DROP_BASIC","DROP_PROWS",
            "DROP_COLUMN","DROP_AREA",
            "DROP_DYNAMIC","DROP_INTERP"),
            tuple<string>("ILUTP(t)","ILUTP(p,t)",
            "Variant of ILUTP(p,t) for j-th column",
            "Variant of ILUTP to control memory",
            "Dynamically adjust threshold",
            "Compute second dropping threshold by interpolation"),
            tuple<SLU::rule_t>(SLU::DROP_BASIC,SLU::DROP_PROWS,SLU::DROP_COLUMN,
            SLU::DROP_AREA,SLU::DROP_DYNAMIC,SLU::DROP_INTERP),
            pl.getRawPtr());*/

    pl->set("ILU_DropTol", 0.0001, "ILUT drop tolerance");

    pl->set("ILU_FillFactor", 10.0, "ILUT fill factor");

    setStringToIntegralParameter<SLU::norm_t>("ILU_Norm", "INF_NORM",
            "Type of norm to use",
            tuple<string>("ONE_NORM","TWO_NORM","INF_NORM"),
            tuple<string>("1-norm","2-norm","inf-norm"),
            tuple<SLU::norm_t>(SLU::ONE_NORM,SLU::TWO_NORM,SLU::INF_NORM),
            pl.getRawPtr());

    setStringToIntegralParameter<SLU::milu_t>("ILU_MILU", "SILU",
            "Type of modified ILU to use",
            tuple<string>("SILU","SMILU_1","SMILU_2","SMILU_3"),
            tuple<string>("Regular ILU","MILU 1","MILU 2","MILU 3"),
            tuple<SLU::milu_t>(SLU::SILU,SLU::SMILU_1,SLU::SMILU_2,
            SLU::SMILU_3),
            pl.getRawPtr());

    pl->set("ILU_FillTol", 0.01, "ILUT fill tolerance");

    pl->set("ILU_Flag", false, "ILU flag: if true, run ILU routines");

    pl->set("IsContiguous", true, "Whether GIDs contiguous");

    // call METIS before SuperLU
    pl->set("UseMetis", false, "Whether to call METIS before SuperLU");
    pl->set("SymmetrizeMetis", true, "Whether to symmetrize matrix before METIS");

    pl->set("Enable_KokkosKernels_TriangularSolves", false, "Whether to use triangular solves.");
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
    pl->set("SpTRSV_Invert_Diag", true, "specify whether to invert diagonal blocks for supernodal sparse-trianguular solve");
    pl->set("SpTRSV_Invert_OffDiag", false, "specify whether to apply diagonal-inversion to off-diagonal blocks for supernodal sparse-trianguular solve");
    pl->set("SpTRSV_U_CSR", true, "specify whether to store U in CSR forma for supernodal sparse-trianguular solve");
    pl->set("SpTRSV_Merge_Supernodes", false, "specify whether to merge supernodes with similar sparsity structures for supernodal sparse-trianguular solve");
    pl->set("SpTRSV_Use_SpMV", false, "specify whether to use spmv for supernodal sparse-trianguular solve");
#endif

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
    Kokkos::resize(host_nzvals_view_, this->globalNumNonZeros_);
    Kokkos::resize(host_rows_view_, this->globalNumNonZeros_);
    Kokkos::resize(host_col_ptr_view_, this->globalNumRows_ + 1);
  }

  int nnz_ret = 0;
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    TEUCHOS_TEST_FOR_EXCEPTION( this->rowIndexBase_ != this->columnIndexBase_,
                        std::runtime_error,
                        "Row and column maps have different indexbase ");

    Util::get_ccs_helper_kokkos_view<
      MatrixAdapter<Matrix>,host_value_type_array,host_ordinal_type_array,
        host_size_type_array>::do_get(this->matrixA_.ptr(),
          host_nzvals_view_, host_rows_view_,
          host_col_ptr_view_, nnz_ret,
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          ARBITRARY,
          this->rowIndexBase_);
  }

  // Get the SLU data type for this type of matrix
  SLU::Dtype_t dtype = type_map::dtype;

  if( this->root_ ){
    TEUCHOS_TEST_FOR_EXCEPTION( nnz_ret != as<int>(this->globalNumNonZeros_),
                        std::runtime_error,
                        "Did not get the expected number of non-zero vals");

    function_map::template create_CompCol_Matrix<host_value_type_array>( &(data_.A),
                                         this->globalNumRows_, this->globalNumCols_,
                                         nnz_ret,
                                         convert_nzvals_, host_nzvals_view_,
                                         host_rows_view_.data(),
                                         host_col_ptr_view_.data(),
                                         SLU::SLU_NC, dtype, SLU::SLU_GE);
  }

  return true;
}

template <class Matrix, class Vector>
void
Superlu<Matrix,Vector>::triangular_solve_factor()
{
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  size_t ld_rhs = this->matrixA_->getGlobalNumRows();

  // convert etree to parents
  SLU::SCformat * Lstore = (SLU::SCformat*)(data_.L.Store);
  int nsuper = 1 + Lstore->nsuper;     // # of supernodal columns
  Kokkos::resize(data_.parents, nsuper);
  for (int s = 0; s < nsuper; s++) {
    int j = Lstore->sup_to_col[s+1]-1; // the last column index of this supernode
    if (data_.etree[j] == static_cast<int>(ld_rhs)) {
        data_.parents(s) = -1;
    } else {
        data_.parents(s) = Lstore->col_to_sup[data_.etree[j]];
    }
  }

  deep_copy_or_assign_view(device_trsv_perm_r_, data_.perm_r); // will use device to permute
  deep_copy_or_assign_view(device_trsv_perm_c_, data_.perm_c); // will use device to permute
  if (data_.options.Equil == SLU::YES) {
    deep_copy_or_assign_view(device_trsv_R_, data_.R); // will use device to scale
    deep_copy_or_assign_view(device_trsv_C_, data_.C); // will use device to scale
  }

  bool condition_flag = false;
  if (data_.options.ConditionNumber == SLU::YES) {
    using STM = Teuchos::ScalarTraits<magnitude_type>;
    const magnitude_type eps = STM::eps ();

    SCformat *Lstore = (SCformat*)(data_.L.Store);
    int nsuper = 1 + Lstore->nsuper;
    int *nb = Lstore->sup_to_col;
    int max_cols = 0;
    for (int i = 0; i < nsuper; i++) {
      if (nb[i+1] - nb[i] > max_cols) {
        max_cols = nb[i+1] - nb[i];
      }
    }

    // when rcond is small, it is ill-conditioned and flag is true
    const magnitude_type multiply_fact (10.0); // larger the value, more likely flag is true, and no invert
    condition_flag = (((double)max_cols * nsuper) * eps * multiply_fact >= data_.rcond);

#ifdef HAVE_AMESOS2_VERBOSE_DEBUG
    int n = data_.perm_r.extent(0);
    std::cout << this->getComm()->getRank()
              << " : anorm = " << data_.anorm << ", rcond = " << data_.rcond << ", n = " << n
              << ", num super cols = " << nsuper << ", max super cols = " << max_cols
              << " -> " << ((double)max_cols * nsuper) * eps / data_.rcond
              << (((double)max_cols * nsuper) * eps * multiply_fact < data_.rcond ? " (okay)" : " (warn)") << std::endl;
#endif
  }

  // Create handles for U and U^T solves
  if (sptrsv_use_spmv_ && !condition_flag) {
    device_khL_.create_sptrsv_handle(
      KokkosSparse::Experimental::SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, ld_rhs, true);
    device_khU_.create_sptrsv_handle(
      KokkosSparse::Experimental::SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG, ld_rhs, false);
  } else {
    device_khL_.create_sptrsv_handle(
      KokkosSparse::Experimental::SPTRSVAlgorithm::SUPERNODAL_DAG, ld_rhs, true);
    device_khU_.create_sptrsv_handle(
      KokkosSparse::Experimental::SPTRSVAlgorithm::SUPERNODAL_DAG, ld_rhs, false);
  }

  // specify whether to store U in CSR format
  device_khU_.set_sptrsv_column_major (!sptrsv_u_in_csr_);

  // specify whether to merge supernodes with similar sparsity structure
  device_khL_.set_sptrsv_merge_supernodes (sptrsv_merge_supernodes_);
  device_khU_.set_sptrsv_merge_supernodes (sptrsv_merge_supernodes_);

  // invert only if flag is not on
  bool sptrsv_invert_diag    = (!condition_flag && sptrsv_invert_diag_);
  bool sptrsv_invert_offdiag = (!condition_flag && sptrsv_invert_offdiag_);

  // specify whether to invert diagonal blocks
  device_khL_.set_sptrsv_invert_diagonal (sptrsv_invert_diag);
  device_khU_.set_sptrsv_invert_diagonal (sptrsv_invert_diag);

  // specify wheather to apply diagonal-inversion to off-diagonal blocks (optional, default is false)
  device_khL_.set_sptrsv_invert_offdiagonal (sptrsv_invert_offdiag);
  device_khU_.set_sptrsv_invert_offdiagonal (sptrsv_invert_offdiag);

  // set etree
  device_khL_.set_sptrsv_etree(data_.parents.data());
  device_khU_.set_sptrsv_etree(data_.parents.data());

  // set permutation
  device_khL_.set_sptrsv_perm(data_.perm_r.data());
  device_khU_.set_sptrsv_perm(data_.perm_c.data());

  int block_size  = -1; // TODO: add needed params
  if (block_size >= 0) {
    std::cout << " Block Size         : " << block_size << std::endl;
    device_khL_.set_sptrsv_diag_supernode_sizes (block_size, block_size);
    device_khU_.set_sptrsv_diag_supernode_sizes (block_size, block_size);
  }

  // Do symbolic analysis
  {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::RCP< Teuchos::Time > SymboTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for SpTrsv symbolic");
      Teuchos::TimeMonitor numFactTimer(*SymboTimer_);
#endif
    KokkosSparse::Experimental::sptrsv_symbolic
      (&device_khL_, &device_khU_, data_.L, data_.U);
  }

  // Do numerical compute
  {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::RCP< Teuchos::Time > NumerTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for SpTrsv numeric");
      Teuchos::TimeMonitor numFactTimer(*NumerTimer_);
#endif
    KokkosSparse::Experimental::sptrsv_compute
      (&device_khL_, &device_khU_, data_.L, data_.U);
  }
#endif // HAVE_AMESOS2_TRIANGULAR_SOLVE
}

template <class Matrix, class Vector>
void
Superlu<Matrix,Vector>::triangular_solve() const
{
#if defined(KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV) && defined(KOKKOSKERNELS_ENABLE_TPL_SUPERLU)
  size_t ld_rhs = device_xValues_.extent(0);
  size_t nrhs = device_xValues_.extent(1);

  Kokkos::resize(device_trsv_rhs_, ld_rhs, nrhs);
  Kokkos::resize(device_trsv_sol_, ld_rhs, nrhs);

  // forward pivot
  auto local_device_bValues = device_bValues_;
  auto local_device_trsv_perm_r = device_trsv_perm_r_;
  auto local_device_trsv_rhs = device_trsv_rhs_;

  if (data_.rowequ) {
    auto local_device_trsv_R = device_trsv_R_;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceExecSpaceType>(0, ld_rhs),
      KOKKOS_LAMBDA(size_t j) {
      for(size_t k = 0; k < nrhs; ++k) {
        local_device_trsv_rhs(local_device_trsv_perm_r(j),k) = local_device_trsv_R(j) * local_device_bValues(j,k);
      }
    });
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceExecSpaceType>(0, ld_rhs),
      KOKKOS_LAMBDA(size_t j) {
      for(size_t k = 0; k < nrhs; ++k) {
        local_device_trsv_rhs(local_device_trsv_perm_r(j),k) = local_device_bValues(j,k);
      }
    });
  }

  for(size_t k = 0; k < nrhs; ++k) { // sptrsv_solve does not batch
    auto sub_sol = Kokkos::subview(device_trsv_sol_, Kokkos::ALL, k);
    auto sub_rhs = Kokkos::subview(device_trsv_rhs_, Kokkos::ALL, k);

    // do L solve= - numeric (only rhs is modified) on the default device/host space
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::RCP< Teuchos::Time > SpLtrsvTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for L solve");
      Teuchos::TimeMonitor solveTimer(*SpLtrsvTimer_);
#endif
      KokkosSparse::Experimental::sptrsv_solve(&device_khL_, sub_sol, sub_rhs);
    }
    // do L^T solve - numeric (only rhs is modified) on the default device/host space
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::RCP< Teuchos::Time > SpUtrsvTimer_ = Teuchos::TimeMonitor::getNewCounter ("Time for U solve");
      Teuchos::TimeMonitor solveTimer(*SpUtrsvTimer_);
#endif
      KokkosSparse::Experimental::sptrsv_solve(&device_khU_, sub_rhs, sub_sol);
    }
  } // end loop over rhs vectors

  // backward pivot
  auto local_device_xValues = device_xValues_;
  auto local_device_trsv_perm_c = device_trsv_perm_c_;
  if (data_.colequ) {
    auto local_device_trsv_C = device_trsv_C_;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceExecSpaceType>(0, ld_rhs),
      KOKKOS_LAMBDA(size_t j) {
      for(size_t k = 0; k < nrhs; ++k) {
        local_device_xValues(j,k) = local_device_trsv_C(j) * local_device_trsv_rhs(local_device_trsv_perm_c(j),k);
      }
    });
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceExecSpaceType>(0, ld_rhs),
      KOKKOS_LAMBDA(size_t j) {
      for(size_t k = 0; k < nrhs; ++k) {
        local_device_xValues(j,k) = local_device_trsv_rhs(local_device_trsv_perm_c(j),k);
      }
    });
  }
#endif // HAVE_AMESOS2_TRIANGULAR_SOLVE
}

template<class Matrix, class Vector>
const char* Superlu<Matrix,Vector>::name = "SuperLU";


} // end namespace Amesos2

#endif  // AMESOS2_SUPERLU_DEF_HPP
