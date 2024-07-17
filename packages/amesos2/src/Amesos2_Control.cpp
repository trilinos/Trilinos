// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
  \file   Amesos2_Control.cpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Jan 16 09:23:38 2010

  \brief  Implementation for Amesos2::Control
*/

#ifndef AMESOS2_CONTROL_CPP
#define AMESOS2_CONTROL_CPP

#include <sstream>
#include "Amesos2_Control.hpp"

namespace Amesos2 {


/**
 * \brief Sets various parameters concerning solver control.
 *
 * \param parameterList A Teuchos::ParameterList which contains the control
 * parameters to set.
 *
 * <br>
 * <b>Supported Parameters:</b>
 * <ul>
 * <li>"Transpose": { \c true | \c false }.  If \c true , the solver instead
 *   solves \f$ A^T X = B \f$.</li>
 * <li>"AddToDiag": \c std::string .  The \c string should be a conforming
 *   representation of the scalar type of the matrix being used.  The scalar
 *   will be added to non-zero diagonal entries of te matrix.  The non-zero
 *   structure will not be affected.</li>
 * <li>"AddZeroToDiag": { \c true | \c false }.  Zero will be added to the
 *   diagonal where diagonal elements are not present.  Not supported for
 *   matrices which are "locked up" by the time Amesos2 sees them.</li>
 * <li>"MatrixProperty": { 0 | 1 | 2 }, where 0 is a general sparse matrix, 1
 *   is a sparse diagonal matrix, and 2 denotes a symmetric sparse matrix.
 *   This parameter is only effective if the underlying solver supports the
 *   distinction.</li>
 * <li>"ScaleMethod": { 0 | 1 | 2 ... }, where 0 denotes no scaling, 1 is the
 *   underlying solver's first method, 2 is the solver's 2nd alternative, and
 *   so on.</li>
 * </ul>
 */
void Control::setControlParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  // Only check for boolean "Transpose" parameters.  If they are of some other
  // type, we will let the solver interfaces deal with them.
  if( parameterList->isType<bool>("Transpose") ){
    useTranspose_ = parameterList->get<bool>("Transpose");
  }

  if( parameterList->isType<bool>("Iterative refinement") ){
    useIterRefine_ = parameterList->get<bool>("Iterative refinement");
  }
  if( parameterList->isType<int>("Number of iterative refinements") ){
    maxNumIterRefines_ = parameterList->get<int>("Number of iterative refinements");
  }
  if( parameterList->isType<bool>("Verboes for iterative refinement") ){
    verboseIterRefine_ = parameterList->get<bool>("Verboes for iterative refinement");
  }

  // Add this value to all diagonal elements which are structurally
  // non-zero. No change is made to non-zero structure of the matrix.
  if( parameterList->isParameter("AddToDiag") ){
    addToDiag_ = parameterList->get<std::string>("AddToDiag");
  }

  // Add zero to diagonal if diagonal element is not present.
  // - not supported for matrices which are missing elements from the diagonal.  See bug #1928 for discussion
  if( parameterList->isParameter("AddZeroToDiag") ){
    addZeroToDiag_ = parameterList->get<bool>("AddZeroToDiag");
  }

  // Matrix property, defined internally in Amesos2::Mumps as an integer,
  // whose value can be:
  // - 0 : general unsymmetric matrix;
  // - 1 : SPD;
  // - 2 : general symmetric matrix.
  if( parameterList->isParameter("MatrixProperty") ) {
    std::string MatrixProperty;
    MatrixProperty = parameterList->get<std::string>("MatrixProperty");
    if( MatrixProperty == "general" )
      matrixProperty_ = 0;
    else if( MatrixProperty == "SPD" )
      matrixProperty_ = 1;
    else if( MatrixProperty == "symmetric" )
      matrixProperty_ = 2;
    else {
      std::ostringstream oss;
      oss << "Amesos2 : ERROR" << std::endl
          << "Amesos2 : matrixProperty value not recognized ("
          << MatrixProperty << ")" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(true,
        std::invalid_argument,
        oss.str());
    }

    // If true, the Amesos2 class should reindex the matrix to
    // standard indexing (i.e. 0-(n-1)).  At present, only
    // Amesos2::Klu supports this option.
    if( parameterList->isParameter("Reindex") ){
      reindex_ = parameterList->get<bool>("Reindex");
    }
  }
}


} // end namespace Amesos2

#endif	// AMESOS2_CONTROL_CPP
