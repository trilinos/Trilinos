// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTOR_STD_OPS_TESTER_DECL_HPP
#define THYRA_VECTOR_STD_OPS_TESTER_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Thyra {

/** \brief Testing class that tests all of the standard vector
 * operations defined in ??? using an arbitrary vector space.
 *
 * ToDo: Finish documentation!
 */
template <class Scalar>
class VectorStdOpsTester {
public:

  /** \brief . */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Set the maximum relative error before a warning is generated. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, warning_tol );

  /** \brief Set the maximum relative error before an error is generated. */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, error_tol );

  /** \brief . */
  VectorStdOpsTester(
    const ScalarMag    &warning_tol = 0
    ,const ScalarMag   &error_tol   = 0
    );

  /** \brief Run the tests using a vector space.
   *
   * @param  vecSpc   [in] VectorBase space used to generate vectors in tests.
   * @param  out      [in/out] If <tt>out!=NULL</tt> then <tt>*out</tt> will
   *                  receive output about the tests.
   * @param  dumpAll  [in] If <tt>true</tt> then vector elements will be printed after
   *                  each transformation operation.  Default is <tt>false</tt>.
   *
   * @return Returns <tt>true</tt> if all of the tests check out and
   * <tt>false</tt> otherwise.
   */
  bool checkStdOps(
    const VectorSpaceBase<Scalar>    &vecSpc
    ,std::ostream                    *out      = 0
    ,const bool                      &dumpAll  = false
    );

};

} // namespace Thyra

#endif // THYRA_VECTOR_STD_OPS_TESTER_DECL_HPP
