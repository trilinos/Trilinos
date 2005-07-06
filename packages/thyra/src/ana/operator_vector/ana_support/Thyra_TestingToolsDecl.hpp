// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_TESTING_TOOLS_DECL_HPP
#define THYRA_TESTING_TOOLS_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"

namespace Thyra {

/** \defgroup Thyra_Op_Vec_test_tools_code_grp Miscellaneous C++ utility code for testing and debugging.

\ingroup Thyra_Op_Vec_ANA_Development_grp

Here is some assorted C++ code to aid in testing and debugging
%Thyra code.

*/

/** \brief Return "passed" or "failed".
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
inline
const char* passfail(bool pass) { return pass ? "passed" : "failed"; }

/** \brief Return relative error of two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template <class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
relErr( const Scalar &s1, const Scalar &s2 );

/** \brief Compute, check and optionally print the relative error in two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template<class Scalar>
bool testRelErr(
  const std::string                                             &v1_name
  ,const Scalar                                                 &v1
  ,const std::string                                            &v2_name
  ,const Scalar                                                 &v2
  ,const std::string                                            &maxRelErr_error_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr_error
  ,const std::string                                            &maxRelErr_warning_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &maxRelErr_warning
  ,std::ostream                                                 *out
  ,const std::string                                            &leadingIndent = std::string("")
  );

/** \brief Print summary outputting for a test or just <tt>passed</tt> or
 * <tt>failed</tt>.
 *
 * @param  result          [in] Bool for of the test was successful or unsuccessful.
 * @param  test_summary    [in] The summary of the test that was just completed.
 * @param  show_all_tests  [in] Bool for if the test summary should be shown even if
 *                         the test passed.
 * @param  success         [out] Update of the success bool.
 * @param  out             [out] Stream where output will be sent if <tt>*out!=NULL</tt>.
 *
 * Preconditions:<ul>
 * <li><tt>success!=NULL</tt>
 * </ul>
 *
 * Preconditions:<ul>
 * <li><tt>*success==false</tt> if <tt>result==false</tt>
 * </ul>
 * 
 * Just the the definition of this function to see what it does.
 */
void printTestResults(
  const bool              result
  ,const std::string      &test_summary
  ,const bool             show_all_tests
  ,bool                   *success
  ,std::ostream           *out
  );

/** \brief Output operator to pretty print any <tt>Thyra::VectorBase</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const VectorBase<Scalar>& v );

/** \brief Output operator to pretty print any <tt>Thyra::LinearOpBase</tt> object.
 *
 * Calls <tt>M.describe(o,Teuchos::VERB_EXTREME);</tt>
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const LinearOpBase<Scalar>& M );

} // namespace Thyra

// //////////////////////////
// Inline functions                        

inline
void Thyra::printTestResults(
  const bool              result
  ,const std::string      &test_summary
  ,const bool             show_all_tests
  ,bool                   *success
  ,std::ostream           *out
  )
{
  if(!result) *success = false;
  if(out) {
    if( !result || show_all_tests )
      *out << std::endl << test_summary;
    else
      *out << "passed!\n";
  }
}

#endif // THYRA_TESTING_TOOLS_DECL_HPP
