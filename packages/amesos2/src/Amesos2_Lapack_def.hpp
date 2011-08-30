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
   \file   Amesos2_Lapack_def.hpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Fri Jul 29 18:54:39 MDT 2011
   
   \brief  Definitions for the Amesos2 Lapack interface.
*/

#ifndef AMESOS2_LAPACK_DEF_HPP
#define AMESOS2_LAPACK_DEF_HPP

#include <Teuchos_RCP.hpp>

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_Lapack_decl.hpp"

namespace Amesos2 {


  template <class Matrix, class Vector>
  Lapack<Matrix,Vector>::Lapack(Teuchos::RCP<const Matrix> A,
				Teuchos::RCP<Vector>       X,
				Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::Lapack,Matrix,Vector>(A, X, B) // instantiate superclass
    , nzvals_()
    , rowind_()
    , colptr_()
  {
    // Set default parameters
    Teuchos::RCP<Teuchos::ParameterList> default_params
      = Teuchos::parameterList( *(this->getValidParameters()) );
    this->setParameters(default_params);
  }


  template <class Matrix, class Vector>
  Lapack<Matrix,Vector>::~Lapack( )
  {
    /*
     * Free any memory allocated by the Lapack library functions (i.e. none)
     */
  }


  template<class Matrix, class Vector>
  int
  Lapack<Matrix,Vector>::preOrdering_impl()
  {
    return(0);
  }


  template <class Matrix, class Vector>
  int
  Lapack<Matrix,Vector>::symbolicFactorization_impl()
  {
    return(0);
  }


  template <class Matrix, class Vector>
  int
  Lapack<Matrix,Vector>::numericFactorization_impl()
  {
    // Set here so that solver_ can refresh it's internal state
    solver_.setMatrix( Teuchos::rcpFromRef(lu_) );

    int factor_ierr = solver_.factor();

    TEST_FOR_EXCEPTION( factor_ierr != 0,
			std::runtime_error,
			"Lapack factor routine returned error code "
			<< factor_ierr );

    return( 0 );
  }


  template <class Matrix, class Vector>
  int
  Lapack<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
				    const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    using Teuchos::as;
    
    // Convert X and B to SerialDenseMatrix's
    const global_size_type ld_rhs = this->root_ ? X->getGlobalLength() : 0;
    const size_t nrhs = X->getGlobalNumVectors();

    const size_t val_store_size = as<size_t>(ld_rhs * nrhs);
    if( this->root_ ){
      rhsvals_.resize(val_store_size);
    }

    {                             // Get values from RHS B
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mvConvTimer( this->timers_.vecConvTime_ );
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif
      Util::get_1d_copy_helper<MultiVecAdapter<Vector>, scalar_type>::do_get(B, rhsvals_(),
									     as<size_t>(ld_rhs),
									     ROOTED);
    }

    if( this->root_ ){
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor solveTimer( this->timers_.solveTime_ );
#endif
      
      using Teuchos::rcpFromRef;
      typedef Teuchos::SerialDenseMatrix<int,scalar_type> DenseMat;
      
      DenseMat rhs_dense_mat(Teuchos::View, rhsvals_.getRawPtr(),
				  as<int>(ld_rhs), as<int>(ld_rhs), as<int>(nrhs));

      solver_.setVectors( rcpFromRef(rhs_dense_mat),
			  rcpFromRef(rhs_dense_mat) );
    
      int solve_ierr = solver_.solve();
      TEST_FOR_EXCEPTION( solve_ierr != 0,
			  std::runtime_error,
			  "Lapack solver solve method returned with error code "
			  << solve_ierr );

      /* Work around a bug in Teuchos::SerialDenseSolver that does not
       * perform the unequilibration of the solution.  If the matrix
       * and RHS were not equilibrated, then this is a no-op
       */
      int unequilibrate_ierr = solver_.unequilibrateLHS();
      TEST_FOR_EXCEPTION( unequilibrate_ierr != 0,
			  std::runtime_error,
			  "Lapack solver returned with error code "
			  << unequilibrate_ierr
			  << " when unequilibrating the system solution." );

      // Solution is found in rhsvals_
    }

    /* Update X's global values */
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor redistTimer( this->timers_.vecRedistTime_ );
#endif

      Util::put_1d_data_helper<
      MultiVecAdapter<Vector>,scalar_type>::do_put(X, rhsvals_(),
						   as<size_t>(ld_rhs),
						   ROOTED);
    }

    return( 0 );
  }


  template <class Matrix, class Vector>
  bool
  Lapack<Matrix,Vector>::matrixShapeOK_impl() const
  {
    // Factorization of rectangular matrices is supported, but not
    // their solution.  For solution we can have square matrices.
    
    return( this->globalNumCols_ == this->globalNumRows_ );
  }


  template <class Matrix, class Vector>
  void
  Lapack<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    solver_.solveWithTranspose( parameterList->get<bool>("Transpose",
							 this->control_.useTranspose_) );

    solver_.factorWithEquilibration( parameterList->get<bool>("Equilibrate", true) );

    // solver_.solveToRefinedSolution( parameterList->get<bool>("Refine", false) );
  }

  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  Lapack<Matrix,Vector>::getValidParameters_impl() const
  {
    using Teuchos::ParameterList;

    static Teuchos::RCP<const Teuchos::ParameterList> valid_params;
    
    if( is_null(valid_params) ){
      Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

      pl->set("Equilibrate", true, "Whether to equilibrate the input matrix");

      // TODO: Refinement does not seem to be working with the SerialDenseSolver.  Not sure why.
      // pl->set("Refine", false, "Whether to apply iterative refinement");

      valid_params = pl;
    }

    return valid_params;
  }

  template <class Matrix, class Vector>
  bool
  Lapack<Matrix,Vector>::loadA_impl(EPhase current_phase)
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // We only load the matrix when numericFactorization is called
    if( current_phase < NUMFACT ) return( false );

    if( this->root_ ){
      nzvals_.resize(this->globalNumNonZeros_);
      rowind_.resize(this->globalNumNonZeros_);
      colptr_.resize(this->globalNumCols_ + 1);
    }

    // global_size_type nnz_ret = 0;
    int nnz_ret = 0;
    {
#ifdef HAVE_AMESOS2_TIMERS
      Teuchos::TimeMonitor mtxRedistTimer( this->timers_.mtxRedistTime_ );
#endif

      // typedef Util::get_ccs_helper<MatrixAdapter<Matrix>,
      // 	scalar_type, global_ordinal_type, global_size_type> ccs_helper;
      typedef Util::get_ccs_helper<MatrixAdapter<Matrix>,
	scalar_type, int, int> ccs_helper;
      ccs_helper::do_get(this->matrixA_.ptr(),
			 nzvals_(), rowind_(), colptr_(),
			 nnz_ret, ROOTED, SORTED_INDICES);
  }

    if( this->root_ ){
      // entries are initialized to zero in here:
      lu_.shape(this->globalNumRows_, this->globalNumCols_);

      // Put entries of ccs representation into the dense matrix
      global_size_type end_col = this->globalNumCols_;
      for( global_size_type col = 0; col < end_col; ++col ){
	global_ordinal_type ptr = colptr_[col];
	global_ordinal_type end_ptr = colptr_[col+1];
	for( ; ptr < end_ptr; ++ptr ){
	  lu_(rowind_[ptr], col) = nzvals_[ptr];
	}
      }

      // lu_.print(std::cout);
    }
    
  return( true );
}


  template<class Matrix, class Vector>
  const char* Lapack<Matrix,Vector>::name = "LAPACK";


} // end namespace Amesos2

#endif	// AMESOS2_LAPACK_DEF_HPP
