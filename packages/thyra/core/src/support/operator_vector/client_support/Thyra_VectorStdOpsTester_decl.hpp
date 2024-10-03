// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
