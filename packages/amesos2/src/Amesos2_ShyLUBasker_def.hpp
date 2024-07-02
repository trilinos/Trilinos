// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  , is_contiguous_(true)
{

  //Nothing

  // Override some default options
  // TODO: use data_ here to init
#if defined(HAVE_AMESOS2_KOKKOS) && defined(KOKKOS_ENABLE_OPENMP)
  /*
  static_assert(std::is_same<kokkos_exe,Kokkos::OpenMP>::value,
  "Kokkos node type not supported by experimental ShyLUBasker Amesos2");
  */
  typedef Kokkos::OpenMP Exe_Space;

  ShyLUbasker = new ::BaskerNS::BaskerTrilinosInterface<local_ordinal_type, shylubasker_dtype, Exe_Space>();
  ShyLUbasker->Options.no_pivot      = BASKER_FALSE;
  ShyLUbasker->Options.static_delayed_pivot = 0;
  ShyLUbasker->Options.symmetric      = BASKER_FALSE;
  ShyLUbasker->Options.realloc        = BASKER_TRUE;
  ShyLUbasker->Options.verbose        = BASKER_FALSE;
  ShyLUbasker->Options.prune          = BASKER_TRUE;
  ShyLUbasker->Options.btf_matching   = 2; // use cardinary matching from Trilinos, globally
  ShyLUbasker->Options.blk_matching   = 1; // use max-weight matching from Basker on each diagonal block
  ShyLUbasker->Options.matrix_scaling = 0; // use matrix scaling on a big A block
  ShyLUbasker->Options.min_block_size = 0; // no merging small blocks
  ShyLUbasker->Options.amd_dom           = BASKER_TRUE;  // use block-wise AMD
  ShyLUbasker->Options.use_metis         = BASKER_TRUE;  // use scotch/metis for ND (TODO: should METIS optional?)
  ShyLUbasker->Options.use_nodeNDP       = BASKER_TRUE;  // use nodeNDP to compute ND partition
  ShyLUbasker->Options.run_nd_on_leaves  = BASKER_TRUE;  // run ND on the final leaf-nodes
  ShyLUbasker->Options.run_amd_on_leaves = BASKER_FALSE; // run AMD on the final leaf-nodes
  ShyLUbasker->Options.transpose     = BASKER_FALSE;
  ShyLUbasker->Options.replace_tiny_pivot = BASKER_FALSE;
  ShyLUbasker->Options.verbose_matrix_out = BASKER_FALSE;

  ShyLUbasker->Options.user_fill     = (double)BASKER_FILL_USER;
  ShyLUbasker->Options.use_sequential_diag_facto = BASKER_FALSE;
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  num_threads = Kokkos::OpenMP::max_hardware_threads();
#else
  num_threads = Kokkos::OpenMP::impl_max_hardware_threads();
#endif

#else
 TEUCHOS_TEST_FOR_EXCEPTION(1 != 0,
     std::runtime_error,
     "Amesos2_ShyLUBasker Exception: Do not have supported Kokkos node type (OpenMP) enabled for ShyLUBasker");
#endif
}


template <class Matrix, class Vector>
ShyLUBasker<Matrix,Vector>::~ShyLUBasker( )
{  
  /* ShyLUBasker will cleanup its own internal memory*/
#if defined(HAVE_AMESOS2_KOKKOS) && defined(KOKKOS_ENABLE_OPENMP)
  ShyLUbasker->Finalize();
  delete ShyLUbasker;
#endif
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

  int info = 0;
  if(this->root_)
  {
    ShyLUbasker->SetThreads(num_threads); 


    // NDE: Special case 
    // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
    // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
    //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic

    if ( single_proc_optimization() ) {

      host_ordinal_type_array sp_rowptr;
      host_ordinal_type_array sp_colind;
      // this needs to be checked during loadA_impl...
      this->matrixA_->returnRowPtr_kokkos_view(sp_rowptr);
      TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr.data() == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null ");
      this->matrixA_->returnColInd_kokkos_view(sp_colind);
      TEUCHOS_TEST_FOR_EXCEPTION(sp_colind.data() == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null ");

      host_value_type_array hsp_values;
      this->matrixA_->returnValues_kokkos_view(hsp_values);
      shylubasker_dtype * sp_values = function_map::convert_scalar(hsp_values.data());
      //shylubasker_dtype * sp_values = function_map::convert_scalar(nzvals_view_.data());
      TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
          std::runtime_error, "Amesos2 Runtime Error: sp_values returned null ");

      // In this case, colptr_, rowind_, nzvals_ are invalid
      info = ShyLUbasker->Symbolic(this->globalNumRows_,
          this->globalNumCols_,
          this->globalNumNonZeros_,
          sp_rowptr.data(),
          sp_colind.data(),
          sp_values,
          true);

      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
          std::runtime_error, "Error in ShyLUBasker Symbolic");
    }
    else 
    { //follow original code path if conditions not met
      // In this case, loadA_impl updates colptr_, rowind_, nzvals_
      shylubasker_dtype * sp_values = function_map::convert_scalar(nzvals_view_.data());
      info = ShyLUbasker->Symbolic(this->globalNumRows_,
          this->globalNumCols_,
          this->globalNumNonZeros_,
          colptr_view_.data(),
          rowind_view_.data(),
          sp_values);

      TEUCHOS_TEST_FOR_EXCEPTION(info != 0,
          std::runtime_error, "Error in ShyLUBasker Symbolic");
    }
  } // end if (this->root_)
  /*No symbolic factoriztion*/

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->matrixA_->getComm()), 0, &info);
  return(info);
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


      // NDE: Special case 
      // Rather than going through the Amesos2 machinery to convert the matrixA_ CRS pointer data to CCS and store in Teuchos::Arrays,
      // in this special case we pass the CRS raw pointers directly to ShyLUBasker which copies+transposes+sorts the data for CCS format
      //   loadA_impl is essentially an empty function in this case, as the raw pointers are handled here and similarly in Symbolic

      if ( single_proc_optimization() ) {

        host_ordinal_type_array sp_rowptr;
        host_ordinal_type_array sp_colind;
        this->matrixA_->returnRowPtr_kokkos_view(sp_rowptr);
        TEUCHOS_TEST_FOR_EXCEPTION(sp_rowptr.data() == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_rowptr returned null ");
        this->matrixA_->returnColInd_kokkos_view(sp_colind);
        TEUCHOS_TEST_FOR_EXCEPTION(sp_colind.data() == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_colind returned null ");

        host_value_type_array hsp_values;
        this->matrixA_->returnValues_kokkos_view(hsp_values);
        shylubasker_dtype * sp_values = function_map::convert_scalar(hsp_values.data());
        //shylubasker_dtype * sp_values = function_map::convert_scalar(nzvals_view_.data());

        TEUCHOS_TEST_FOR_EXCEPTION(sp_values == nullptr,
            std::runtime_error, "Amesos2 Runtime Error: sp_values returned null ");

        // In this case, colptr_, rowind_, nzvals_ are invalid
        info = ShyLUbasker->Factor( this->globalNumRows_,
            this->globalNumCols_,
            this->globalNumNonZeros_,
            sp_rowptr.data(),
            sp_colind.data(),
            sp_values);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, 
            std::runtime_error, "Error ShyLUBasker Factor");
      }
      else 
      {
        // In this case, loadA_impl updates colptr_, rowind_, nzvals_
        shylubasker_dtype * sp_values = function_map::convert_scalar(nzvals_view_.data());
        info = ShyLUbasker->Factor(this->globalNumRows_,
            this->globalNumCols_,
            this->globalNumNonZeros_,
            colptr_view_.data(),
            rowind_view_.data(),
            sp_values);

        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, 
            std::runtime_error, "Error ShyLUBasker Factor");
      }

      //ShyLUbasker->DEBUG_PRINT();

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
  TEUCHOS_TEST_FOR_EXCEPTION(info == -1,
    std::runtime_error,
    "ShyLUBasker: Could not alloc space for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION(info == -2,
    std::runtime_error,
    "ShyLUBasker: Could not alloc needed work space");
  TEUCHOS_TEST_FOR_EXCEPTION(info == -3,
    std::runtime_error,
    "ShyLUBasker: Could not alloc additional memory needed for L and U");
  TEUCHOS_TEST_FOR_EXCEPTION(info > 0,
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
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor solveTimer(this->timers_.solveTime_);
#endif

  int ierr = 0; // returned error code

  using Teuchos::as;

  const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
  const size_t nrhs = X->getGlobalNumVectors();

  const bool ShyluBaskerTransposeRequest = this->control_.useTranspose_;
  const bool initialize_data = true;
  const bool do_not_initialize_data = false;

  if ( single_proc_optimization() && nrhs == 1 ) {

    // no msp creation
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(initialize_data, B, bValues_, as<size_t>(ld_rhs));

    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(do_not_initialize_data, X, xValues_, as<size_t>(ld_rhs));

  } // end if ( single_proc_optimization() && nrhs == 1 )
  else {

    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(initialize_data, B, bValues_,
        as<size_t>(ld_rhs),
        (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
        this->rowIndexBase_);

    // See Amesos2_Tacho_def.hpp for notes on why we 'get' x here.
    Util::get_1d_copy_helper_kokkos_view<MultiVecAdapter<Vector>,
      host_solve_array_t>::do_get(do_not_initialize_data, X, xValues_,
        as<size_t>(ld_rhs),
        (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
        this->rowIndexBase_);
  }

  if ( this->root_ ) { // do solve
    shylubasker_dtype * pxValues = function_map::convert_scalar(xValues_.data());
    shylubasker_dtype * pbValues = function_map::convert_scalar(bValues_.data());
    if (!ShyluBaskerTransposeRequest)
      ierr = ShyLUbasker->Solve(nrhs, pbValues, pxValues);
    else
      ierr = ShyLUbasker->Solve(nrhs, pbValues, pxValues, true);
  }

  /* All processes should have the same error code */
  Teuchos::broadcast(*(this->getComm()), 0, &ierr);

  TEUCHOS_TEST_FOR_EXCEPTION( ierr  > 0,
      std::runtime_error,
      "Encountered zero diag element at: " << ierr);
  TEUCHOS_TEST_FOR_EXCEPTION( ierr == -1,
      std::runtime_error,
      "Could not alloc needed working memory for solve" );

  Util::put_1d_data_helper_kokkos_view<
    MultiVecAdapter<Vector>,host_solve_array_t>::do_put(X, xValues_,
      as<size_t>(ld_rhs),
      (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED);

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
  if(parameterList->isParameter("delayed pivot"))
    {
      ShyLUbasker->Options.static_delayed_pivot = (parameterList->get<int>("delayed pivot"));
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
  if(parameterList->isParameter("btf"))
    {
      ShyLUbasker->Options.btf = parameterList->get<bool>("btf");
    }
  if(parameterList->isParameter("use_metis"))
    {
      ShyLUbasker->Options.use_metis = parameterList->get<bool>("use_metis");
    }
  if(parameterList->isParameter("use_nodeNDP"))
    {
      ShyLUbasker->Options.use_nodeNDP = parameterList->get<bool>("use_nodeNDP");
    }
  if(parameterList->isParameter("run_nd_on_leaves"))
    {
      ShyLUbasker->Options.run_nd_on_leaves = parameterList->get<bool>("run_nd_on_leaves");
    }
  if(parameterList->isParameter("run_amd_on_leaves"))
    {
      ShyLUbasker->Options.run_amd_on_leaves = parameterList->get<bool>("run_amd_on_leaves");
    }
  if(parameterList->isParameter("amd_on_blocks"))
    {
      ShyLUbasker->Options.amd_dom = parameterList->get<bool>("amd_on_blocks");
    }
  if(parameterList->isParameter("transpose"))
    {
      // NDE: set transpose vs non-transpose mode as bool; track separate shylubasker objects
      const auto transpose = parameterList->get<bool>("transpose");
      if (transpose == true)
        this->control_.useTranspose_ = true;
    }
  if(parameterList->isParameter("use_sequential_diag_facto"))
    {
      ShyLUbasker->Options.use_sequential_diag_facto = parameterList->get<bool>("use_sequential_diag_facto");
    }
  if(parameterList->isParameter("user_fill"))
    {
      ShyLUbasker->Options.user_fill = parameterList->get<double>("user_fill");
    }
  if(parameterList->isParameter("prune"))
    {
      ShyLUbasker->Options.prune = parameterList->get<bool>("prune");
    }
  if(parameterList->isParameter("replace_tiny_pivot"))
    {
      ShyLUbasker->Options.replace_tiny_pivot = parameterList->get<bool>("replace_tiny_pivot");
    }
  if(parameterList->isParameter("btf_matching"))
    {
      ShyLUbasker->Options.btf_matching = parameterList->get<int>("btf_matching");
      if (ShyLUbasker->Options.btf_matching == 1 || ShyLUbasker->Options.btf_matching == 2) {
        ShyLUbasker->Options.matching = true;
      } else {
        ShyLUbasker->Options.matching = false;
      }
    }
  if(parameterList->isParameter("blk_matching"))
    {
      ShyLUbasker->Options.blk_matching = parameterList->get<int>("blk_matching");
    }
  if(parameterList->isParameter("matrix_scaling"))
    {
      ShyLUbasker->Options.matrix_scaling = parameterList->get<int>("matrix_scaling");
    }
  if(parameterList->isParameter("min_block_size"))
    {
      ShyLUbasker->Options.min_block_size = parameterList->get<int>("min_block_size");
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
              "Should not pivot");
      pl->set("delayed pivot", 0,
              "Apply static delayed pivot on a big block");
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
      pl->set("btf", true, 
              "Use BTF ordering");
      pl->set("prune", false,
              "Use prune on BTF blocks (Not Supported)");
      pl->set("btf_matching",  2, 
              "Matching option for BTF: 0 = none, 1 = Basker, 2 = Trilinos (default), (3 = MC64 if enabled)");
      pl->set("blk_matching", 1, 
              "Matching optioon for block: 0 = none, 1 or anything else = Basker (default), (2 = MC64 if enabled)");
      pl->set("matrix_scaling", 0, 
              "Use matrix scaling to biig A BTF block: 0 = no-scaling, 1 = symmetric diagonal scaling, 2 = row-max, and then col-max scaling");
      pl->set("min_block_size",  0, 
              "Size of the minimum diagonal blocks");
      pl->set("replace_tiny_pivot",  false, 
              "Replace tiny pivots during the numerical factorization");
      pl->set("use_metis", true,
              "Use METIS for ND");
      pl->set("use_nodeNDP", true,
              "Use nodeND to compute ND partition");
      pl->set("run_nd_on_leaves", false,
              "Run ND on the final leaf-nodes for ND factorization");
      pl->set("run_amd_on_leaves", false,
              "Run AMD on the final leaf-nodes for ND factorization");
      pl->set("amd_on_blocks", true,
              "Run AMD on each diagonal blocks");
      pl->set("transpose", false,
              "Solve the transpose A");
      pl->set("use_sequential_diag_facto", false,
              "Use sequential algorithm to factor each diagonal block");
      pl->set("user_fill", (double)BASKER_FILL_USER,
              "User-provided padding for the fill ratio");
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
      Kokkos::resize(nzvals_view_, this->globalNumNonZeros_);
      Kokkos::resize(rowind_view_, this->globalNumNonZeros_);
      Kokkos::resize(colptr_view_, this->globalNumCols_ + 1); //this will be wrong for case of gapped col ids, e.g. 0,2,4,9; num_cols = 10 ([0,10)) but num GIDs = 4...
    }

    local_ordinal_type nnz_ret = 0;
    {
    #ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
    #endif

      Util::get_ccs_helper_kokkos_view<
        MatrixAdapter<Matrix>, host_value_type_array, host_ordinal_type_array, host_ordinal_type_array>
        ::do_get(this->matrixA_.ptr(), nzvals_view_, rowind_view_, colptr_view_, nnz_ret,
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          ARBITRARY,
          this->rowIndexBase_); // copies from matrixA_ to ShyLUBasker ConcreteSolver cp, ri, nzval members
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
