/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// @HEADER
*/

// /////////////////////////////////////////////
// Teuchos_BLAS_types.hpp

#ifndef TEUCHOS_BLAS_TYPES_HPP
#define TEUCHOS_BLAS_TYPES_HPP

/*! \file Teuchos_BLAS_types.hpp
	\brief Enumerated types for BLAS input characters.
*/

/*! \defgroup BLASEnum_grp Enumerations for character inputs in Teuchos::BLAS methods

  \brief These enumerated lists are used in compile time checking of the input characters
  for BLAS methods.  

	\note Any other input other than those specified here will result
	in an error at compile time and are not supported by the templated BLAS/LAPACK interface.

	<ul>
	<li><b>Teuchos::ESide</b> : Enumerated list for BLAS character input "SIDE".
		<ul>
		<li>LEFT_SIDE : The matrix/std::vector is on, or applied to, the left side of the equation
		<li>RIGHT_SIDE : The matrix/std::vector is on, or applied to, the right side of the equation
		</ul><br>
	<li><b>Teuchos::ETransp</b> : Enumerated list for BLAS character input "TRANS".
		<ul>
		<li>NO_TRANS : The matrix/std::vector is not transposed
		<li>TRANS : The matrix/std::vector is transposed
		<li>CONJ_TRANS : The matrix/std::vector is conjugate transposed
		</ul><br>
	<li><b>Teuchos::EUplo</b> : Enumerated list for BLAS character input "UPLO".
		<ul>
		<li>UPPER_TRI : The matrix is upper triangular
		<li>LOWER_TRI : The matrix is lower triangular
		</ul><br>
	<li><b>Teuchos::EDiag</b> : Enumerated list for BLAS character input "DIAG".
		<ul>
		<li>UNIT_DIAG : The matrix has all ones on its diagonal
		<li>NON_UNIT_DIAG : The matrix does not have all ones on its diagonal
		</ul><br>
        </ul>
*/

namespace Teuchos {
  enum ESide { 	
    LEFT_SIDE,	/*!< Left side */ 
    RIGHT_SIDE 	/*!< Right side */
  };

  enum ETransp { 	
    NO_TRANS,	/*!< Not transposed */ 
    TRANS, 		/*!< Transposed */
    CONJ_TRANS 	/*!< Conjugate transposed */
  };
  
  enum EUplo { 	
    UPPER_TRI,	/*!< Upper triangular */ 
    LOWER_TRI,	/*!< Lower triangular */
    UNDEF_TRI   /*!< Unspeficied/undefined triangular structure */
  };
  
  enum EDiag { 	
    UNIT_DIAG,	/*!< Unit diagaonal */ 
    NON_UNIT_DIAG	/*!< Not unit diagonal */ 
  };

  enum EType {
    FULL,	/*!< Full matrix */
    LOWER,	/*!< Lower triangular */
    UPPER,	/*!< Upper triangular */
    HESSENBERG, /*!< Upper Hessenberg */
    SYM_BAND_L, /*!< Symmetric band, lower half stored */
    SYM_BAND_U, /*!< Symmetric band, upper half stored */
    BAND        /*!< General band */
  };
}

#endif // TEUCHOS_BLAS_TYPES_HPP
