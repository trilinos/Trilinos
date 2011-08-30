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
      TEST_FOR_EXCEPTION(true,
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
