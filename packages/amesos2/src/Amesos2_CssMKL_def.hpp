// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_CssMKL_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:52:00 MDT 2011

   \brief  Definitions for the Amesos2 CssMKL interface.
*/

#ifndef AMESOS2_CSSMKL_DEF_HPP
#define AMESOS2_CSSMKL_DEF_HPP

#include <map>

#include <Teuchos_Tuple.hpp>
#include <Teuchos_toString.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_CssMKL_decl.hpp"


namespace Amesos2 {

  namespace PMKL {
#   include <mkl.h>
#   include <mkl_pardiso.h>
  }

  template <class Matrix, class Vector>
  CssMKL<Matrix,Vector>::CssMKL(Teuchos::RCP<const Matrix> A,
                                        Teuchos::RCP<Vector>       X,
                                        Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::CssMKL,Matrix,Vector>(A, X, B) // instantiate superclass
    , n_(Teuchos::as<int_t>(this->globalNumRows_))
    , perm_(this->globalNumRows_)
    , nrhs_(0)
    , css_initialized_(false)
    , is_contiguous_(true)
  {
    // set the default matrix type
    set_css_mkl_matrix_type();
    set_css_mkl_default_parameters(pt_, iparm_);

    // index base
    const global_ordinal_type indexBase = this->matrixA_->getRowMap ()->getIndexBase ();
    iparm_[34] = (indexBase == 0 ? 1 : 0);  /* Use one or zero-based indexing */
    // 1D block-row distribution
    auto frow  = this->matrixA_->getRowMap()->getMinGlobalIndex();
    auto nrows = this->matrixA_->getLocalNumRows();
    iparm_[39] = 2;  /* Matrix input format. */
    iparm_[40] = frow;          /* > Beginning of input domain. */
    iparm_[41] = frow+nrows-1;  /* > End of input domain. */

    // get MPI Comm
    Teuchos::RCP<const Teuchos::Comm<int> > matComm = this->matrixA_->getComm ();
    TEUCHOS_TEST_FOR_EXCEPTION(
        matComm.is_null (), std::logic_error, "Amesos2::CssMKL "
        "constructor: The matrix's communicator is null!");
    Teuchos::RCP<const Teuchos::MpiComm<int> > matMpiComm =
      Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> > (matComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      matMpiComm.is_null (), std::logic_error, "Amesos2::CssMKL "
      "constructor: The matrix's communicator is not an MpiComm!");
    TEUCHOS_TEST_FOR_EXCEPTION(
      matMpiComm->getRawMpiComm ().is_null (), std::logic_error, "Amesos2::"
      "CssMKL constructor: The matrix's communicator claims to be a "
      "Teuchos::MpiComm<int>, but its getRawPtrComm() method returns "
      "Teuchos::null!  This means that the underlying MPI_Comm doesn't even "
      "exist, which likely implies that the Teuchos::MpiComm was constructed "
      "incorrectly.  It means something different than if the MPI_Comm were "
      "MPI_COMM_NULL.");
    MPI_Comm CssComm = *(matMpiComm->getRawMpiComm ());
    CssComm_ = MPI_Comm_c2f(CssComm);

    // rowmap for loadA (to have locally contiguous)
    css_rowmap_ =
      Teuchos::rcp (new map_type (this->globalNumRows_, nrows, indexBase, matComm));
  }


  template <class Matrix, class Vector>
  CssMKL<Matrix,Vector>::~CssMKL( )
  {
    /*
     * Free any memory allocated by the CssMKL library functions
     */
    int_t error = 0;
    if (css_initialized_)
    {
      int_t phase = -1;         // release all internal solver memory
      void *bdummy, *xdummy;
      const MPI_Fint CssComm = CssComm_;
      function_map::cluster_sparse_solver( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
      css_initialized_ = false;
    }
    check_css_mkl_error(Amesos2::CLEAN, error);
  }


  template<class Matrix, class Vector>
  int
  CssMKL<Matrix,Vector>::preOrdering_impl()
  {
    // preOrdering done during "Analysis" (aka symbolic
    // factorization) phase
    return(0);
  }


  template <class Matrix, class Vector>
  int
  CssMKL<Matrix,Vector>::symbolicFactorization_impl()
  {
    int_t error = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor symbFactTimer( this->timers_.symFactTime_ );
#endif

      int_t phase = 11; // Analysis
      void *bdummy, *xdummy;
      const MPI_Fint CssComm = CssComm_;
      function_map::cluster_sparse_solver( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
    }
    check_css_mkl_error(Amesos2::SYMBFACT, error);

    // Pardiso only lets you retrieve the total number of factor
    // non-zeros, not for each individually.  We should document how
    // such a situation is reported.
    this->setNnzLU(iparm_[17]);
    css_initialized_ = true;

    return(0);
  }


  template <class Matrix, class Vector>
  int
  CssMKL<Matrix,Vector>::numericFactorization_impl()
  {
    int_t error = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer( this->timers_.numFactTime_ );
#endif

      //int_t phase = 12; // Analysis, numerical factorization
      int_t phase = 22; // Numerical factorization
      void *bdummy, *xdummy;
      const MPI_Fint CssComm = CssComm_;
      function_map::cluster_sparse_solver( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &CssComm, &error );
    }
    check_css_mkl_error(Amesos2::NUMFACT, error);

    return( 0 );
  }


  template <class Matrix, class Vector>
  int
  CssMKL<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                    const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    using Teuchos::as;

    // Get B data
    const local_ordinal_type ld_rhs = this->matrixA_->getLocalNumRows();
    nrhs_ = as<int_t>(X->getGlobalNumVectors());

    const size_t val_store_size = as<size_t>(ld_rhs * nrhs_);
    xvals_.resize(val_store_size);
    bvals_.resize(val_store_size);
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer( this->timers_.vecConvTime_ );
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      Util::get_1d_copy_helper<
        MultiVecAdapter<Vector>,
        solver_scalar_type>::do_get(B, bvals_(),
          as<size_t>(ld_rhs),
          DISTRIBUTED_NO_OVERLAP,
          this->rowIndexBase_);
    }

    int_t error = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer( this->timers_.solveTime_ );
#endif

      const int_t phase = 33; // Solve, iterative refinement
      const MPI_Fint CssComm = CssComm_;
      function_map::cluster_sparse_solver( pt_,
                             const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_),
                             const_cast<int_t*>(&mtype_),
                             const_cast<int_t*>(&phase),
                             const_cast<int_t*>(&n_),
                             const_cast<solver_scalar_type*>(nzvals_view_.data()),
                             const_cast<int_t*>(rowptr_view_.data()),
                             const_cast<int_t*>(colind_view_.data()),
                             const_cast<int_t*>(perm_.getRawPtr()),
                             &nrhs_,
                             const_cast<int_t*>(iparm_),
                             const_cast<int_t*>(&msglvl_),
                             as<void*>(bvals_.getRawPtr()),
                             as<void*>(xvals_.getRawPtr()), &CssComm, &error );
    }
    check_css_mkl_error(Amesos2::SOLVE, error);

    /* Get values to X */
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,
        solver_scalar_type>::do_put(X, xvals_(),
          as<size_t>(ld_rhs),
          DISTRIBUTED_NO_OVERLAP);
    }

    return( 0 );
}


  template <class Matrix, class Vector>
  bool
  CssMKL<Matrix,Vector>::matrixShapeOK_impl() const
  {
    // CssMKL supports square matrices
    return( this->globalNumRows_ == this->globalNumCols_ );
  }


  template <class Matrix, class Vector>
  void
  CssMKL<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    using Teuchos::RCP;
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterEntryValidator;

    RCP<const Teuchos::ParameterList> valid_params = getValidParameters_impl();

    // Fill-in reordering: 0 = minimum degree, 2 = METIS 4.0.1 (default), 3 = METIS 5.1, 4 = AMD,
    if( parameterList->isParameter("IPARM(2)") )
    {
      RCP<const ParameterEntryValidator> fillin_validator = valid_params->getEntry("IPARM(2)").validator();
      parameterList->getEntry("IPARM(2)").setValidator(fillin_validator);
      iparm_[1] = getIntegralValue<int>(*parameterList, "IPARM(2)");
    }

    // Max numbers of iterative refinement steps
    if( parameterList->isParameter("IPARM(8)") )
    {
      RCP<const ParameterEntryValidator> refine_validator = valid_params->getEntry("IPARM(8)").validator();
      parameterList->getEntry("IPARM(8)").setValidator(refine_validator);
      iparm_[7] = getIntegralValue<int>(*parameterList, "IPARM(8)");
    }

    // Perturb the pivot elements
    if( parameterList->isParameter("IPARM(10)") )
    {
      RCP<const ParameterEntryValidator> pivot_perturb_validator = valid_params->getEntry("IPARM(10)").validator();
      parameterList->getEntry("IPARM(10)").setValidator(pivot_perturb_validator);
      iparm_[9] = getIntegralValue<int>(*parameterList, "IPARM(10)");
    }

    // First check if the control object requests a transpose solve.
    // Then solver specific options can override this.
    iparm_[11] = this->control_.useTranspose_ ? 2 : 0;
    // Normal solve (0), or a transpose solve (1)
    if( parameterList->isParameter("IPARM(12)") )
    {
      RCP<const ParameterEntryValidator> trans_validator = valid_params->getEntry("IPARM(12)").validator();
      parameterList->getEntry("IPARM(12)").setValidator(trans_validator);
      iparm_[11] = getIntegralValue<int>(*parameterList, "IPARM(12)");
    }

    // (Non-)symmetric matchings : detault 1 for nonsymmetric and 0 for symmetric matrix (default is nonsymmetric)
    if( parameterList->isParameter("IPARM(13)") )
    {
      RCP<const ParameterEntryValidator> trans_validator = valid_params->getEntry("IPARM(13)").validator();
      parameterList->getEntry("IPARM(13)").setValidator(trans_validator);
      iparm_[12] = getIntegralValue<int>(*parameterList, "IPARM(13)");
    }

    // Output: Number of nonzeros in the factor LU
    if( parameterList->isParameter("IPARM(18)") )
    {
      RCP<const ParameterEntryValidator> report_validator = valid_params->getEntry("IPARM(18)").validator();
      parameterList->getEntry("IPARM(18)").setValidator(report_validator);
      iparm_[17] = getIntegralValue<int>(*parameterList, "IPARM(18)");
    }
   
    if( parameterList->isParameter("IsContiguous") ){
      is_contiguous_ = parameterList->get<bool>("IsContiguous");
    }
  }


/*
 * TODO: It would be nice if the parameters could be expressed as
 * either all string or as all integers.  I see no way of doing this
 * at present with the standard validators.  However, we could create
 * our own validators or kindly ask the Teuchos team to add some
 * features for use.
 *
 * The issue is that with the current validators we cannot specify
 * arbitrary sets of numbers that are the only allowed parameters.
 * For example the IPARM(2) parameter can take only the values 0, 2,
 * and 3.  The EnhancedNumberValidator can take a min value, and max
 * value, and a step size, but with those options there is no way to
 * specify the needed set.
 *
 * Another missing feature is the ability to give docstrings for such
 * numbers.  For example IPARM(25) can take on the values 0 and 1.
 * This would be easy enough to accomplish with just a number
 * validator, but then have no way to document the effect of each
 * value.
 */
template <class Matrix, class Vector>
Teuchos::RCP<const Teuchos::ParameterList>
CssMKL<Matrix,Vector>::getValidParameters_impl() const
{
  using std::string;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::tuple;
  using Teuchos::toString;
  using Teuchos::EnhancedNumberValidator;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::anyNumberParameterEntryValidator;

  static Teuchos::RCP<const Teuchos::ParameterList> valid_params;

  if( is_null(valid_params) ){
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

    void* pt_temp[64];
    int_t iparm_temp[64];
    set_css_mkl_default_parameters(pt_temp, iparm_temp);
    setStringToIntegralParameter<int>("IPARM(2)", toString(iparm_temp[1]),
                                      "Fill-in reducing ordering for the input matrix",
                                      tuple<string>("2", "3", "10"),
                                      tuple<string>("Nested dissection algorithm from METIS",
                                      "Parallel version of the nested dissection algorithm",
                                      "MPI version of the nested dissection and symbolic factorization algorithms"),
                                      tuple<int>(2, 3, 10),
                                      pl.getRawPtr());
    
    setStringToIntegralParameter<int>("IPARM(12)", toString(iparm_temp[11]),
                                      "Solve with transposed or conjugate transposed matrix A",
                                      tuple<string>("0", "1", "2"),
                                      tuple<string>("Non-transposed",
                                      "Conjugate-transposed",
                                      "Transposed"),
                                      tuple<int>(0, 1, 2),
                                      pl.getRawPtr());
 
    setStringToIntegralParameter<int>("IPARM(13)", toString(iparm_temp[12]),
                                      "Use weighted matching",
                                      tuple<string>("0", "1"),
                                      tuple<string>("No matching", "Use matching"),
                                      tuple<int>(0, 1),
                                      pl.getRawPtr());

    Teuchos::AnyNumberParameterEntryValidator::EPreferredType preferred_int =
      Teuchos::AnyNumberParameterEntryValidator::PREFER_INT;

    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes accept_int( false );
    accept_int.allowInt( true );

    pl->set("IPARM(8)" , as<int>(iparm_temp[7]) , "Iterative refinement step",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IPARM(10)", as<int>(iparm_temp[9]) , "Pivoting perturbation",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IPARM(18)", as<int>(iparm_temp[17]), "Report the number of non-zero elements in the factors",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IsContiguous", true, "Whether GIDs contiguous");

    valid_params = pl;
  }

  return valid_params;
}



template <class Matrix, class Vector>
bool
CssMKL<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

  // CssMKL does not need matrix data in the pre-ordering phase
  if( current_phase == PREORDERING ) return( false );

  EDistribution dist_option = (iparm_[39] != 0 ? DISTRIBUTED_NO_OVERLAP : ((is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED));
  if (current_phase == SYMBFACT) {
    if (dist_option == DISTRIBUTED_NO_OVERLAP) {
      Kokkos::resize(nzvals_temp_, this->matrixA_->getLocalNNZ());
      Kokkos::resize(nzvals_view_, this->matrixA_->getLocalNNZ());
      Kokkos::resize(colind_view_, this->matrixA_->getLocalNNZ());
      Kokkos::resize(rowptr_view_, this->matrixA_->getLocalNumRows() + 1);
    } else {
      if( this->root_ ) {
        Kokkos::resize(nzvals_temp_, this->matrixA_->getGlobalNNZ());
        Kokkos::resize(nzvals_view_, this->matrixA_->getGlobalNNZ());
        Kokkos::resize(colind_view_, this->matrixA_->getGlobalNNZ());
        Kokkos::resize(rowptr_view_, this->matrixA_->getGlobalNumRows() + 1);
      } else {
        Kokkos::resize(nzvals_temp_, 0);
        Kokkos::resize(nzvals_view_, 0);
        Kokkos::resize(colind_view_, 0);
        Kokkos::resize(rowptr_view_, 0);
      }
    }
  }

  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif
    int_t nnz_ret = 0;
    Util::get_crs_helper_kokkos_view<MatrixAdapter<Matrix>,
      host_value_type_array,host_ordinal_type_array, host_size_type_array >::do_get(
                                       this->matrixA_.ptr(),
                                       nzvals_temp_, colind_view_, rowptr_view_,
                                       nnz_ret,
                                       Teuchos::ptrInArg(*css_rowmap_),
                                       dist_option,
                                       SORTED_INDICES);
    Kokkos::deep_copy(nzvals_view_, nzvals_temp_);
  }
  return( true );
}


template <class Matrix, class Vector>
void
CssMKL<Matrix,Vector>::check_css_mkl_error(EPhase phase,
                                                   int_t error) const
{
  int error_i = error;
  Teuchos::broadcast(*(this->getComm()), 0, &error_i); // We only care about root's value

  if( error == 0 ) return;      // No error

  std::string errmsg = "Other error";
  switch( error ){
  case -1:
    errmsg = "CssMKL reported error: 'Input inconsistent'";
    break;
  case -2:
    errmsg = "CssMKL reported error: 'Not enough memory'";
    break;
  case -3:
    errmsg = "CssMKL reported error: 'Reordering problem'";
    break;
  case -4:
    errmsg =
      "CssMKL reported error: 'Zero pivot, numerical "
      "factorization or iterative refinement problem'";
    break;
  case -5:
    errmsg = "CssMKL reported error: 'Unclassified (internal) error'";
    break;
  case -6:
    errmsg = "CssMKL reported error: 'Reordering failed'";
    break;
  case -7:
    errmsg = "CssMKL reported error: 'Diagonal matrix is singular'";
    break;
  case -8:
    errmsg = "CssMKL reported error: '32-bit integer overflow problem'";
    break;
  case -9:
    errmsg = "CssMKL reported error: 'Not enough memory for OOC'";
    break;
  case -10:
    errmsg = "CssMKL reported error: 'Problems with opening OOC temporary files'";
    break;
  case -11:
    errmsg = "CssMKL reported error: 'Read/write problem with OOC data file'";
    break;
  }
  errmsg += (" at phase = "+std::to_string(phase));

  TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, errmsg );
}


template <class Matrix, class Vector>
void
CssMKL<Matrix,Vector>::set_css_mkl_matrix_type(int_t mtype)
{
  if( mtype == 0 ){
    if( complex_ ){
      mtype_ = 13;              // complex, unsymmetric
    } else {
      mtype_ = 11;              // real, unsymmetric
    }
  } else {
    switch( mtype ){
    case 11:
      TEUCHOS_TEST_FOR_EXCEPTION( complex_,
                          std::invalid_argument,
                          "Cannot set a real Pardiso matrix type with scalar type complex" );
      mtype_ = 11; break;
    case 13:
      TEUCHOS_TEST_FOR_EXCEPTION( !complex_,
                          std::invalid_argument,
                          "Cannot set a complex Pardiso matrix type with non-complex scalars" );
      mtype_ = 13; break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( true,
                          std::invalid_argument,
                          "Symmetric matrices are not yet supported by the Amesos2 interface" );
    }
  }
}

template <class Matrix, class Vector>
void
CssMKL<Matrix,Vector>::set_css_mkl_default_parameters(void* pt[], int_t iparm[]) const
{
  for( int i = 0; i < 64; ++i ){
    pt[i] = nullptr;
    iparm[i] = 0;
  }
  iparm[0] = 1; /* No solver default */
  // Reset some of the default parameters
  iparm[1] = 10;  /* 2: Fill-in reordering from METIS, 3: thread dissection, 10: MPI version of the nested dissection and symbolic factorization*/
  iparm[7] = 0;   /* Max numbers of iterative refinement steps */
  iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
  iparm[10] = 0;  /* Disable nonsymmetric permutation and scaling MPS */
  iparm[11] = 0;  /* Normal solve (0), or a transpose solve (1) */
  iparm[12] = 0;  /* Do not use (non-)symmetric matchings */
  iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  iparm[20] = -1; /* Pivoting for symmetric indefinite matrices */
  iparm[26] = 1;  /* Check input matrix is sorted */

  // set single or double precision
  if constexpr ( std::is_same_v<solver_magnitude_type, PMKL::_REAL_t> ) {
    iparm[27] = 1;           // single-precision
  } else {
    iparm[27] = 0;           // double-precision
  }
  iparm[34] = 1;  /* Use zero-based indexing */
}

template <class Matrix, class Vector>
const char* CssMKL<Matrix,Vector>::name = "CSSMKL";

template <class Matrix, class Vector>
const typename CssMKL<Matrix,Vector>::int_t
CssMKL<Matrix,Vector>::msglvl_ = 0;  // set to be one, for more CSS messages

template <class Matrix, class Vector>
const typename CssMKL<Matrix,Vector>::int_t
CssMKL<Matrix,Vector>::maxfct_ = 1;

template <class Matrix, class Vector>
const typename CssMKL<Matrix,Vector>::int_t
CssMKL<Matrix,Vector>::mnum_ = 1;


} // end namespace Amesos

#endif  // AMESOS2_CSSMKL_DEF_HPP
