// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_NewSolver_def.hpp
   \author John Doe <jd@sandia.gov>
   \date   Thu Jul  8 16:08:10 2010
   
   \brief  Definitions for the Amesos2 NewSolver interface.
*/

#ifndef AMESOS2_NEWSOLVER_DEF_HPP
#define AMESOS2_NEWSOLVER_DEF_HPP

#include "Amesos2_SolverCore_def.hpp"
#include "Amesos2_NewSolver_decl.hpp"

namespace Amesos2 {


  template <class Matrix, class Vector>
  NewSolver<Matrix,Vector>::NewSolver(
				      Teuchos::RCP<const Matrix> A,
				      Teuchos::RCP<Vector>       X,
				      Teuchos::RCP<const Vector> B)
    : SolverCore<Amesos2::NewSolver,Matrix,Vector>(A, X, B) // instantiate superclass
    , nzvals_()
    , colind_()
    , rowptr_()
  {
    /* Assign private variables necessary for interfacing with the NewSolver TPL.
     *
     * Also set default solver options and paramters.
     */
  }


  template <class Matrix, class Vector>
  NewSolver<Matrix,Vector>::~NewSolver( )
  {
    /*
     * Free any memory allocated by the NewSolver library functions
     */
  }


  template<class Matrix, class Vector>
  int
  NewSolver<Matrix,Vector>::preOrdering_impl()
  {
    /*
     * TODO: implement for the NewSolver library
     *
     * Can assume that the internal matrix data is current
     */

    return(0);
  }


  template <class Matrix, class Vector>
  int
  NewSolver<Matrix,Vector>::symbolicFactorization_impl()
  {
    /* TODO: implement for the NewSolver library
     *
     * Can assume that the internal matrix data is current
     */

    return(0);
  }


  template <class Matrix, class Vector>
  int
  NewSolver<Matrix,Vector>::numericFactorization_impl()
  {
    /* TODO: implement for NewSolver TPL
     *
     * Use function_map to link to the appropriate NewSolver library
     * functions in the case that the TPL functions are not overloaded
     * by type.
     *
     * Can assume that the internal matrix data is current
     */

    return( 0 );
  }


  template <class Matrix, class Vector>
  int
  NewSolver<Matrix,Vector>::solve_impl(const Teuchos::Ptr<MultiVecAdapter<Vector> >       X,
				       const Teuchos::Ptr<const MultiVecAdapter<Vector> > B) const
  {
    /* TODO: implement for the NewSolver library
     *
     * Use function_map to link to the appropriate TPL functions in
     * the case that the TPL functions are not overloaded by type.
     * 
     * Multivector data should be extracted from B and solution placed
     * in X before return.  Can use utilities
     * Amesos2::Util::get_1d_copy_helper and
     * Amesos2::Util::put_1d_data_helper if the solver manipulates
     * RHS's as 1D arrays.
     */

    return( 0 );
  }


  template <class Matrix, class Vector>
  bool
  NewSolver<Matrix,Vector>::matrixShapeOK_impl() const
  {
    /*
     * Test shape of Matrix and return appropriate boolean
     */
  
    return( true );
  }


  template <class Matrix, class Vector>
  void
  NewSolver<Matrix,Vector>::setParameters_impl(const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
  {
    /*
     * Set solver-specific parameters
     *
     * Say, for example, the solver supports the options
     *
     * - Equilibrate
     * - PivotThreshold
     * - Transpose
     */

    if( parameterList->isParameter("Equilibrate") ){
      bool equil = parameterList->template get<bool>("Equil");
      if( equil ){
	// Set flag for solver to do equilibration
      } else {
	// Set flag for solver *not* to do equilibration
      }
    }

    if( parameterList->isParameter("DiagPivotThresh") ){
      double threshold = parameterList->template get<double>("DiagPivotThresh");
      // Set solver threshold parameter
    }

    if ( parameterList->isParameter("Trans") ){
      std::string fact = parameterList->template get<std::string>("Trans");
      if( fact == "Tranpose" ){
	// Set flag for solver to use transpose
      } else if ( fact == "NoTranspose" ){
	// Set flag for solver to *not* use transpose
      } else if ( fact == "Conjugate" ) {
	// Set flag for solver to use the conjugate transpose
      }
    } else {                      
      // Set default value (if not set in constructor)
    }
  }

  template <class Matrix, class Vector>
  Teuchos::RCP<const Teuchos::ParameterList>
  NewSolver<Matrix,Vector>::getValidParameters_impl() const
  {
    using Teuchos::ParameterList;

    ParameterList valid_params;

    // valid_params.set("Parameter", default_value, doc_string);

    return Teuchos::rcpFromRef( valid_params );
  }

  template <class Matrix, class Vector>
  bool
  NewSolver<Matrix,Vector>::loadA_impl(EPhase current_phase)
  {
#ifdef HAVE_AMESOS2_TIMERS
    Teuchos::TimeMonitor convTimer(this->timers_.mtxConvTime_);
#endif

    // resize storage arrays if necessary

    // use Amesos2::Util::get_crs_helper or
    // Amesos2::Util::get_ccs_helper to get a CRS or CCS
    // representation of the matrix, if appropriate.  Otherwise do
    // whatever conversion is necessary to get the matrix data ready
    // for the other library functions.
    
    return( true );
  }


  template<class Matrix, class Vector>
  const char* NewSolver<Matrix,Vector>::name = "NewSolver";


} // end namespace Amesos2

#endif	// AMESOS2_NEWSOLVER_DEF_HPP
