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
    , nzvals_()
    , colind_()
    , rowptr_()
    , n_(Teuchos::as<int_t>(this->globalNumRows_))
    , perm_(this->globalNumRows_)
    , nrhs_(0)
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
    if( Meta::is_same<solver_magnitude_type, PMKL::_REAL_t>::value ){
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
                             nzvals_.getRawPtr(), rowptr_.getRawPtr(),
                             colind_.getRawPtr(), perm_.getRawPtr(), &nrhs_, iparm_,
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
                             nzvals_.getRawPtr(), rowptr_.getRawPtr(),
                             colind_.getRawPtr(), perm_.getRawPtr(), &nrhs_, iparm_,
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
                             nzvals_.getRawPtr(), rowptr_.getRawPtr(),
                             colind_.getRawPtr(), perm_.getRawPtr(), &nrhs_, iparm_,
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
                                    ROOTED);
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
                             const_cast<solver_scalar_type*>(nzvals_.getRawPtr()),
                             const_cast<int_t*>(rowptr_.getRawPtr()),
                             const_cast<int_t*>(colind_.getRawPtr()),
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
                                    ROOTED);
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
  iparm_[1]  = validators[2]->getIntegralValue(*parameterList, "IPARM(2)",
					       validators[2]->getDefaultParameterName());
  iparm_[3]  = parameterList->get<int>("IPARM(4)" , (int)iparm_[3]);
  iparm_[7]  = parameterList->get<int>("IPARM(8)" , (int)iparm_[7]);
  iparm_[9]  = parameterList->get<int>("IPARM(10)", (int)iparm_[9]);
  iparm_[17] = parameterList->get<int>("IPARM(18)", (int)iparm_[17]);
  iparm_[23] = validators[24]->getIntegralValue(*parameterList, "IPARM(24)",
						validators[24]->getDefaultParameterName());
  iparm_[24] = validators[25]->getIntegralValue(*parameterList, "IPARM(25)",
						validators[25]->getDefaultParameterName());
  iparm_[59] = validators[60]->getIntegralValue(*parameterList, "IPARM(60)",
						validators[60]->getDefaultParameterName());
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
  using Teuchos::stringToIntegralParameterEntryValidator;
  typedef Teuchos::StringToIntegralParameterEntryValidator<int> STIPEV;
  Teuchos::AnyNumberParameterEntryValidator::EPreferredType preferred_int =
    Teuchos::AnyNumberParameterEntryValidator::PREFER_INT;

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

    // Initialize our parameter validators, saving the string to int validators for later
    RCP<STIPEV> iparm_2_validator
      = stringToIntegralParameterEntryValidator<int>(tuple<string>("0", "2", "3"),
						     tuple<string>("The minimum degree algorithm",
								   "Nested dissection algorithm from METIS",
								   "OpenMP parallel nested dissection algorithm"),
						     tuple<int>(0, 2, 3),
						     toString(iparm_temp[1]));
    validators.insert( std::pair<int,RCP<STIPEV> >(2, iparm_2_validator) );
    
    Teuchos::RCP<EnhancedNumberValidator<int> > iparm_4_validator
      = Teuchos::rcp( new EnhancedNumberValidator<int>() );
    iparm_4_validator->setMin(0);

    RCP<STIPEV> iparm_24_validator
      = stringToIntegralParameterEntryValidator<int>(tuple<string>("0", "1"),
						     tuple<string>("PARDISO uses the previous algorithm for factorization",
								   "PARDISO uses the new two-level factorization algorithm"),
						     tuple<int>(0, 1),
						     toString(iparm_temp[23]));
    validators.insert( std::pair<int,RCP<STIPEV> >(24, iparm_24_validator) );

    RCP<STIPEV> iparm_25_validator
      = stringToIntegralParameterEntryValidator<int>(tuple<string>("0", "1"),
						     tuple<string>("PARDISO uses the parallel algorithm for the solve step",
								   "PARDISO uses the sequential forward and backward solve"),
						     tuple<int>(0, 1),
						     toString(iparm_temp[24]));
    validators.insert( std::pair<int,RCP<STIPEV> >(25, iparm_25_validator) );

    RCP<STIPEV> iparm_60_validator
      = stringToIntegralParameterEntryValidator<int>(tuple<string>("0", "2"),
						     tuple<string>("In-core PARDISO",
								   "Out-of-core PARDISO.  The OOC PARDISO can solve very "
								   "large problems by holding the matrix factors in files "
								   "on the disk. Hence the amount of RAM required by OOC "
								   "PARDISO is significantly reduced."),
						     tuple<int>(0, 2),
						     toString(iparm_temp[59]));
    validators.insert( std::pair<int,RCP<STIPEV> >(60, iparm_60_validator) );

    Teuchos::AnyNumberParameterEntryValidator::AcceptedTypes accept_int( false );
    accept_int.allowInt( true );

    pl->set("IPARM(2)" , validators[2]->getDefaultParameterName(),
	    "Fill-in reducing ordering for the input matrix", validators[2]);

    pl->set("IPARM(4)" , as<int>(iparm_temp[3]) , "Preconditioned CGS/CG",
            iparm_4_validator);

    pl->set("IPARM(8)" , as<int>(iparm_temp[8]) , "Iterative refinement step",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IPARM(10)", as<int>(iparm_temp[9]) , "Pivoting perturbation",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IPARM(18)", as<int>(iparm_temp[17]), "Report the number of non-zero elements in the factors",
            anyNumberParameterEntryValidator(preferred_int, accept_int));

    pl->set("IPARM(24)", validators[24]->getDefaultParameterName(),
	    "Parallel factorization control", validators[24]);
    
    pl->set("IPARM(25)", validators[25]->getDefaultParameterName(),
	    "Parallel forward/backward solve control", validators[25]);

    pl->set("IPARM(60)", validators[60]->getDefaultParameterName(),
	    "PARDISO mode (OOC mode)", validators[60]);

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
    nzvals_.resize(this->globalNumNonZeros_);
    colind_.resize(this->globalNumNonZeros_);
    rowptr_.resize(this->globalNumRows_ + 1);
  }

  int_t nnz_ret = 0;
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

    Util::get_crs_helper<
    MatrixAdapter<Matrix>,
      solver_scalar_type,
      int_t,int_t>::do_get(this->matrixA_.ptr(),
                           nzvals_(), colind_(), rowptr_(),
                           nnz_ret, ROOTED, SORTED_INDICES);
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

  TEST_FOR_EXCEPTION( true, std::runtime_error, errmsg );
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
      TEST_FOR_EXCEPTION( complex_,
                          std::invalid_argument,
                          "Cannot set a real Pardiso matrix type with scalar type complex" );
      mtype_ = 11; break;
    case 13:
      TEST_FOR_EXCEPTION( !complex_,
                          std::invalid_argument,
                          "Cannot set a complex Pardiso matrix type with non-complex scalars" );
      mtype_ = 13; break;
    default:
      TEST_FOR_EXCEPTION( true,
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
