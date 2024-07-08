// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/**
   \file   Amesos2_PardisoMKL_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Wed Jul 27 12:52:00 MDT 2011

   \brief  Definitions for the Amesos2 PardisoMKL interface.
*/

#ifndef AMESOS2_PARDISOMKL_DEF_HPP
#define AMESOS2_PARDISOMKL_DEF_HPP

#include <map>

#include <Teuchos_Tuple.hpp>
#include <Teuchos_toString.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_PardisoMKL_decl.hpp"


namespace Amesos2 {

  namespace PMKL {
#   include <mkl.h>
#   include <mkl_pardiso.h>
  }

  template <class Matrix, class Vector>
  PardisoMKL<Matrix,Vector>::PardisoMKL(Teuchos::RCP<const Matrix> A,
                                        Teuchos::RCP<Vector>       X,
                                        Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::PardisoMKL,Matrix,Vector>(A, X, B) // instantiate superclass
    , n_(Teuchos::as<int_t>(this->globalNumRows_))
    , perm_(this->globalNumRows_)
    , nrhs_(0)
    , is_contiguous_(true)
  {
    // set the default matrix type
    set_pardiso_mkl_matrix_type();

    PMKL::_INTEGER_t iparm_temp[64];
    PMKL::_INTEGER_t mtype_temp = mtype_;
    PMKL::pardisoinit(pt_, &mtype_temp, iparm_temp);

    for( int i = 0; i < 64; ++i ){
      iparm_[i] = iparm_temp[i];
    }

    // set single or double precision
    if constexpr ( std::is_same_v<solver_magnitude_type, PMKL::_REAL_t> ) {
      iparm_[27] = 1;           // single-precision
    } else {
      iparm_[27] = 0;           // double-precision
    }

    // Reset some of the default parameters
    iparm_[34] = 1;             // Use zero-based indexing
#ifdef HAVE_AMESOS2_DEBUG
    iparm_[26] = 1;             // turn the Pardiso matrix checker on
#endif
  }


  template <class Matrix, class Vector>
  PardisoMKL<Matrix,Vector>::~PardisoMKL( )
  {
    /*
     * Free any memory allocated by the PardisoMKL library functions
     */
    int_t error = 0;
    void *bdummy, *xdummy;

    if( this->root_ ){
      int_t phase = -1;         // release all internal solver memory
      function_map::pardiso( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &error );
    }

    check_pardiso_mkl_error(Amesos2::CLEAN, error);
  }


  template<class Matrix, class Vector>
  int
  PardisoMKL<Matrix,Vector>::preOrdering_impl()
  {
    // preOrdering done in PardisoMKL during "Analysis" (aka symbolic
    // factorization) phase

    return(0);
  }


  template <class Matrix, class Vector>
  int
  PardisoMKL<Matrix,Vector>::symbolicFactorization_impl()
  {
    int_t error = 0;

    if( this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor symbFactTimer( this->timers_.symFactTime_ );
#endif

      int_t phase = 11;
      void *bdummy, *xdummy;

      function_map::pardiso( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &error );
    }

    check_pardiso_mkl_error(Amesos2::SYMBFACT, error);

    // Pardiso only lets you retrieve the total number of factor
    // non-zeros, not for each individually.  We should document how
    // such a situation is reported.
    this->setNnzLU(iparm_[17]);

    return(0);
  }


  template <class Matrix, class Vector>
  int
  PardisoMKL<Matrix,Vector>::numericFactorization_impl()
  {
    int_t error = 0;

    if( this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor numFactTimer( this->timers_.numFactTime_ );
#endif

      int_t phase = 22;
      void *bdummy, *xdummy;

      function_map::pardiso( pt_, const_cast<int_t*>(&maxfct_),
                             const_cast<int_t*>(&mnum_), &mtype_, &phase, &n_,
                             nzvals_view_.data(), rowptr_view_.data(),
                             colind_view_.data(), perm_.getRawPtr(), &nrhs_, iparm_,
                             const_cast<int_t*>(&msglvl_), &bdummy, &xdummy, &error );
    }

    check_pardiso_mkl_error(Amesos2::NUMFACT, error);

    return( 0 );
  }


  template <class Matrix, class Vector>
  int
  PardisoMKL<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
                                        const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    using Teuchos::as;

    int_t error = 0;

    // Get B data
    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    nrhs_ = as<int_t>(X->getGlobalNumVectors());

    const size_t val_store_size = as<size_t>(ld_rhs * nrhs_);
    xvals_.resize(val_store_size);
    bvals_.resize(val_store_size);

    {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer( this->timers_.vecConvTime_ );
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      Util::get_1d_copy_helper<
        MultiVecAdapter<Vector>,
        solver_scalar_type>::do_get(B, bvals_(),
          as<size_t>(ld_rhs),
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
          this->rowIndexBase_);
    }

    if( this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer( this->timers_.solveTime_ );
#endif

      const int_t phase = 33;

      function_map::pardiso( pt_,
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
                             as<void*>(xvals_.getRawPtr()), &error );
    }

    check_pardiso_mkl_error(Amesos2::SOLVE, error);

    /* Export X from root to the global space */
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer(this->timers_.vecRedistTime_);
#endif

      Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,
        solver_scalar_type>::do_put(X, xvals_(),
          as<size_t>(ld_rhs),
          (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED);
  }

    return( 0 );
}


  template <class Matrix, class Vector>
  bool
  PardisoMKL<Matrix,Vector>::matrixShapeOK_impl() const
  {
    // PardisoMKL supports square matrices
    return( this->globalNumRows_ == this->globalNumCols_ );
  }


  template <class Matrix, class Vector>
  void
  PardisoMKL<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
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

    // Iterative-direct algorithm
    if( parameterList->isParameter("IPARM(4)") )
    {
      RCP<const ParameterEntryValidator> prec_validator = valid_params->getEntry("IPARM(4)").validator();
      parameterList->getEntry("IPARM(4)").setValidator(prec_validator);
      iparm_[3] = getIntegralValue<int>(*parameterList, "IPARM(4)");
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
   
    // Scheduling method for the parallel numerical factorization
    if( parameterList->isParameter("IPARM(24)") )
    {
      RCP<const ParameterEntryValidator> par_fact_validator = valid_params->getEntry("IPARM(24)").validator();
      parameterList->getEntry("IPARM(24)").setValidator(par_fact_validator);
      iparm_[23] = getIntegralValue<int>(*parameterList, "IPARM(24)");
    }
  
    // Parallelization scheme for the forward and backward solv
    if( parameterList->isParameter("IPARM(25)") )
    {
      RCP<const ParameterEntryValidator> par_fbsolve_validator = valid_params->getEntry("IPARM(25)").validator();
      parameterList->getEntry("IPARM(25)").setValidator(par_fbsolve_validator);
      iparm_[24] = getIntegralValue<int>(*parameterList, "IPARM(25)");
    } 

    // Graph compression scheme for METIS.
    if( parameterList->isParameter("IPARM(60)") )
    {
      RCP<const ParameterEntryValidator> ooc_validator = valid_params->getEntry("IPARM(60)").validator();
      parameterList->getEntry("IPARM(60)").setValidator(ooc_validator);
      iparm_[59] = getIntegralValue<int>(*parameterList, "IPARM(60)");
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
PardisoMKL<Matrix,Vector>::getValidParameters_impl() const
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

    // Use pardisoinit to get some default values;
    void *pt_dummy[64];
    PMKL::_INTEGER_t mtype_temp = mtype_;
    PMKL::_INTEGER_t iparm_temp[64];
    PMKL::pardisoinit(pt_dummy,
                      const_cast<PMKL::_INTEGER_t*>(&mtype_temp),
                      const_cast<PMKL::_INTEGER_t*>(iparm_temp));

    setStringToIntegralParameter<int>("IPARM(2)", toString(iparm_temp[1]),
                                      "Fill-in reducing ordering for the input matrix",
                                      tuple<string>("0", "2", "3"),
                                      tuple<string>("The minimum degree algorithm",
                                      "Nested dissection algorithm from METIS",
                                      "OpenMP parallel nested dissection algorithm"),
                                      tuple<int>(0, 2, 3),
                                      pl.getRawPtr());
    
    Teuchos::RCP<EnhancedNumberValidator<int> > iparm_4_validator
      = Teuchos::rcp( new EnhancedNumberValidator<int>() );
    iparm_4_validator->setMin(0);
    pl->set("IPARM(4)" , as<int>(iparm_temp[3]) , "Preconditioned CGS/CG",
            iparm_4_validator);

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

    setStringToIntegralParameter<int>("IPARM(24)", toString(iparm_temp[23]),
                                      "Parallel factorization control",
                                      tuple<string>("0", "1"),
                                      tuple<string>("PARDISO uses the previous algorithm for factorization",
                                      "PARDISO uses the new two-level factorization algorithm"),
                                      tuple<int>(0, 1),
                                      pl.getRawPtr());

    setStringToIntegralParameter<int>("IPARM(25)", toString(iparm_temp[24]),
                                      "Parallel forward/backward solve control",
                                      tuple<string>("0", "1"),
                                      tuple<string>("PARDISO uses the parallel algorithm for the solve step",
                                      "PARDISO uses the sequential forward and backward solve"),
                                      tuple<int>(0, 1),
                                      pl.getRawPtr());

    setStringToIntegralParameter<int>("IPARM(60)", toString(iparm_temp[59]),
                                      "PARDISO mode (OOC mode)",
                                      tuple<string>("0", "2"),
                                      tuple<string>("In-core PARDISO",
                                      "Out-of-core PARDISO.  The OOC PARDISO can solve very "
                                      "large problems by holding the matrix factors in files "
                                      "on the disk. Hence the amount of RAM required by OOC "
                                      "PARDISO is significantly reduced."),
                                      tuple<int>(0, 2),
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
PardisoMKL<Matrix,Vector>::loadA_impl(EPhase current_phase)
{
#ifdef HAVE_AMESOS2_TIMERS
  Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

  // PardisoMKL does not need matrix data in the pre-ordering phase
  if( current_phase == PREORDERING ) return( false );

  if( this->root_ ){
    Kokkos::resize(nzvals_view_, this->globalNumNonZeros_);
    Kokkos::resize(colind_view_, this->globalNumNonZeros_);
    Kokkos::resize(rowptr_view_, this->globalNumRows_ + 1);
  }

  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    int_t nnz_ret = 0;
    Util::get_crs_helper_kokkos_view<
      MatrixAdapter<Matrix>,
      host_value_type_array, host_ordinal_type_array, host_size_type_array>::do_get(
        this->matrixA_.ptr(),
        nzvals_view_, colind_view_, rowptr_view_, nnz_ret, 
        (is_contiguous_ == true) ? ROOTED : CONTIGUOUS_AND_ROOTED,
        SORTED_INDICES,
        this->rowIndexBase_);
  }

  return( true );
}


template <class Matrix, class Vector>
void
PardisoMKL<Matrix,Vector>::check_pardiso_mkl_error(EPhase phase,
                                                   int_t error) const
{
  int error_i = error;
  Teuchos::broadcast(*(this->getComm()), 0, &error_i); // We only care about root's value

  if( error == 0 ) return;      // No error

  std::string errmsg = "Other error";
  switch( error ){
  case -1:
    errmsg = "PardisoMKL reported error: 'Input inconsistent'";
    break;
  case -2:
    errmsg = "PardisoMKL reported error: 'Not enough memory'";
    break;
  case -3:
    errmsg = "PardisoMKL reported error: 'Reordering problem'";
    break;
  case -4:
    errmsg =
      "PardisoMKL reported error: 'Zero pivot, numerical "
      "factorization or iterative refinement problem'";
    break;
  case -5:
    errmsg = "PardisoMKL reported error: 'Unclassified (internal) error'";
    break;
  case -6:
    errmsg = "PardisoMKL reported error: 'Reordering failed'";
    break;
  case -7:
    errmsg = "PardisoMKL reported error: 'Diagonal matrix is singular'";
    break;
  case -8:
    errmsg = "PardisoMKL reported error: '32-bit integer overflow problem'";
    break;
  case -9:
    errmsg = "PardisoMKL reported error: 'Not enough memory for OOC'";
    break;
  case -10:
    errmsg = "PardisoMKL reported error: 'Problems with opening OOC temporary files'";
    break;
  case -11:
    errmsg = "PardisoMKL reported error: 'Read/write problem with OOC data file'";
    break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, errmsg );
}


template <class Matrix, class Vector>
void
PardisoMKL<Matrix,Vector>::set_pardiso_mkl_matrix_type(int_t mtype)
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
const char* PardisoMKL<Matrix,Vector>::name = "PARDISOMKL";

template <class Matrix, class Vector>
const typename PardisoMKL<Matrix,Vector>::int_t
PardisoMKL<Matrix,Vector>::msglvl_ = 0;

template <class Matrix, class Vector>
const typename PardisoMKL<Matrix,Vector>::int_t
PardisoMKL<Matrix,Vector>::maxfct_ = 1;

template <class Matrix, class Vector>
const typename PardisoMKL<Matrix,Vector>::int_t
PardisoMKL<Matrix,Vector>::mnum_ = 1;


} // end namespace Amesos

#endif  // AMESOS2_PARDISOMKL_DEF_HPP
