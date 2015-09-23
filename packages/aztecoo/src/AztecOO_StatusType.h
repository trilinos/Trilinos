
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef AZTECOO_STATUSTYPE_H
#define AZTECOO_STATUSTYPE_H
/*! \file AztecOO_StatusType.h 
    \brief AztecOO StatusType: Used to return convergence status information for AztecOO_StatusTest objects.
 */

/*! \enum AztecOO_StatusType 
    When the CheckStatus and GetStatus methods of AztecOO_StatusTest objects are called a variable
    of type AztecOO_StatusType is returned.
*/



enum AztecOO_StatusType { Unchecked = 2,        /*!< Initial state of status */
			  Unconverged = 1,      /*!< Convergence is not reached. */
			  Converged = 0,        /*!< Convergence is reached. */
			  Failed = -1,          /*!< Some failure occured.  Should stop */
			  NaN = -2,              /*!< Result from test contains a NaN value.  Should stop */
			  PartialFailed = -3   /*!< Some success, e.g. implicit residual converged but explicit did not. */
			  
};
#endif /* AZTECOO_STATUSTYPE_H */
