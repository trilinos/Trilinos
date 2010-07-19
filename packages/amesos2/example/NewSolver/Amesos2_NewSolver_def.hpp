// @HEADER
// ***********************************************************************
//
//                Amesos2: Direct Sparse Solver Package
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


#ifndef AMESOS2_NEWSOLVER_DEF_HPP
#define AMESOS2_NEWSOLVER_DEF_HPP

#include "Amesos2_Solver.hpp"
#include "Amesos2_NewSolver_MatrixHelper.hpp"
#include "Amesos2_NewSolver_TypeMap.hpp"
#include "Amesos2_NewSolver_FunctionMap.hpp"


namespace Amesos {


template <class Matrix, class Vector>
NewSolver<Matrix,Vector>::NewSolver(
  Teuchos::RCP<Matrix> A,
  Teuchos::RCP<Vector> X,
  Teuchos::RCP<Vector> B)
  : Solver<Amesos::NewSolver,Matrix,Vector>(A, X, B) // instantiate superclass
  , nzvals_(this->globalNumNonZeros_)                // initialize private members
  , rowind_(this->globalNumNonZeros_)
  , colptr_(this->globalNumRows_ + 1)
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
   * Free any memory allocated by the TPL
   */
}


template<class Matrix, class Vector>
int
NewSolver<Matrix,Vector>::preOrdering_impl()
{
  /* TODO: implement for NewSolver TPL */

  return(0);
}


template <class Matrix, class Vector>
int
NewSolver<Matrix,Vector>::symbolicFactorization_impl()
{
  /* TODO: implement for NewSolver TPL
   *
   * The Matrix A should be refreshed each time.
   */

  return(0);
}


template <class Matrix, class Vector>
int
NewSolver<Matrix,Vector>::numericFactorization_impl()
{
  /* TODO: implement for NewSolver TPL
   *
   * Use Amesos::TypeMap<Amesos::NewSolver,scalar_type> to link to the
   * appropriate TPL functions in the case that the TPL functions are not
   * overloaded by type.
   *
   * The Matrix A should be refreshed each time.
   */

  return(0);
}


template <class Matrix, class Vector>
int
NewSolver<Matrix,Vector>::solve_impl()
{
  /* TODO: implement for NewSolver TPL
   *
   * Use Amesos::TypeMap<Amesos::NewSolver,scalar_type> to link to the
   * appropriate TPL functions in the case that the TPL functions are not
   * overloaded by type.
   * 
   * The MultiVectors X and B should be refreshed each time.
   */

  return(0);
}


template <class Matrix, class Vector>
bool
NewSolver<Matrix,Vector>::matrixShapeOK_impl() const
{
  /*
   * Test shape of Matrix and return appropriate boolean
   */
  
  return(true);
}


template <class Matrix, class Vector>
void
NewSolver<Matrix,Vector>::setParameters_impl(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
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


template<class Matrix, class Vector>
const char* NewSolver<Matrix,Vector>::name = "Amesos2::NewSolver";


} // end namespace Amesos

#endif	// AMESOS2_NEWSOLVER_DEF_HPP
