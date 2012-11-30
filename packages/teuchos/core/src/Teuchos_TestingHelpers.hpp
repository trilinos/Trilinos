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

#ifndef TEUCHOS_TESTING_HELPERS_HPP
#define TEUCHOS_TESTING_HELPERS_HPP


/*! \file Teuchos_TestingHelpers.hpp
    \brief Utilities to make writing tests easier.
*/


#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_FancyOStream.hpp"


namespace Teuchos {


/** \brief Update the success bool flag.
 *
 * \ingroup teuchos_testing_grp
 */
inline void updateSuccess(const bool result, bool &success);


/** \brief Return "passed" or "failed".
 *
 * \ingroup teuchos_testing_grp
 */
inline const std::string passfail(const bool result);


/** \brief Helper function for TEUCHOS_PASS_FAIL(...).
 *
 * \ingroup teuchos_testing_grp
 */
TEUCHOSCORE_LIB_DLL_EXPORT const std::string passfail_with_location(const bool result, const std::string &file, const int lineNumber);

/** \brief Set if TEUCHOS_PASS_FAIL(...) should print test failure location.
 *
 * \ingroup teuchos_testing_grp
 */
void showTestFailureLocation(bool);


/** \brief Return if TEUCHOS_PASS_FAIL(...) should print test failure location.
 *
 * \ingroup teuchos_testing_grp
 */
bool showTestFailureLocation();


/** \brief .
 *
 * \ingroup teuchos_testing_grp
 */
template <bool hasMachineParameters, class Scalar>
class RelErrSmallNumber {
public:
  static Scalar smallNumber()
    {
      return ScalarTraits<Scalar>::ThisShouldNotCompile();
    }
};


/** \brief .
 *
 * \ingroup teuchos_testing_grp
 */
template <class Scalar>
class RelErrSmallNumber<false,Scalar> {
public:
  static Scalar smallNumber()
    {
      return Scalar(1e-8);
    }
};


/** \brief .
 *
 * \ingroup teuchos_testing_grp
 */
template <class Scalar>
class RelErrSmallNumber<true,Scalar> {
public:
  static Scalar smallNumber()
    {
      return Teuchos::ScalarTraits<Scalar>::eps();
    }
};


/** \brief . 
 *
 * \ingroup teuchos_testing_grp
 */
template <class Scalar>
Scalar defaultSmallNumber()
{
  const bool hasMachineParameters = ScalarTraits<Scalar>::hasMachineParameters;
  return RelErrSmallNumber<hasMachineParameters,Scalar>::smallNumber();
}


/** \brief Return relative error of two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup teuchos_testing_grp
 */
template <class Scalar>
typename ScalarTraits<Scalar>::magnitudeType
relErr( const Scalar &s1, const Scalar &s2 );


/** \brief Compute, check and optionally print the relative error in two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template<class Scalar>
bool testRelErr(
  const std::string &v1_name,
  const Scalar &v1,
  const std::string &v2_name,
  const Scalar &v2,
  const std::string &maxRelErr_error_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_error,
  const std::string &maxRelErr_warning_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_warning,
  const Ptr<std::ostream> &out
  );


/** \brief Compare if two array objects are the same or not.
 *
 * This function works with any two array objects are the same size and have
 * the same element value types.  The funtion is templated on the container
 * types and therefore can compare any two objects that have size() and
 * operator[](i) defined.
 *
 * \returns Returns <tt>true</tt> if the compare and <tt>false</tt> otherwise.
 *
 * \ingroup teuchos_testing_grp
 */
template<class Array1, class Array2>
bool compareArrays(
  const Array1 &a1, const std::string &a1_name,
  const Array2 &a2, const std::string &a2_name,
  Teuchos::FancyOStream &out
  );


/** \brief Compare if two array objects are the same or not up to a relative
 * floating point precision.
 *
 * This function works with any two array objects are the same size and have
 * the same element value types.  The funtion is templated on the container
 * types and therefore can compare any two objects that have size() and
 * operator[](i) defined.
 *
 * \returns Returns <tt>true</tt> if the compare and <tt>false</tt> otherwise.
 *
 * \ingroup teuchos_testing_grp
 */
template<class Array1, class Array2, class ScalarMag>
bool compareFloatingArrays(
  const Array1 &a1, const std::string &a1_name,
  const Array2 &a2, const std::string &a2_name,
  const ScalarMag &tol,
  Teuchos::FancyOStream &out
  );


} // namespace Teuchos


/** \brief Macro that prints "passed" or "failed" and optionally prints the
 * file name and line number as well.
 *
 * Only prints the file name and line number if
 * Teuchos::showTestFailureLocation() == true.
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_PASS_FAIL(RESULT) \
  Teuchos::passfail_with_location((RESULT), __FILE__, __LINE__)


/** \brief Echo a statement and then invoke it.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_ECHO( statement, out ) \
  (out) << #statement ";\n"; \
  statement;

/** \brief Test that an object is equal to a given constant.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success ) \
  { \
    (out) << #v1" = "<<Teuchos::toString(v1)<<" == "<<Teuchos::toString(v2)<<" : "; \
    const bool l_result = (v1) == (v2); \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    if (!l_result) (success) = false; \
  }

/** \brief Assert that a give object is true.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_ASSERT( v1, out, success ) \
  { \
    const bool l_result = v1; \
    (out) << #v1" = "<<l_result<<" == true : "; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    if (!l_result) (success) = false; \
  }

/** \brief Test that two values are equal.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_EQUALITY( v1, v2, out, success ) \
  { \
    (out) << #v1" = "<<Teuchos::toString(v1)<<" == "#v2" = "<<Teuchos::toString(v2)<<" : "; \
    const bool l_result = (v1) == (v2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test that an object is not equal to a given constant.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_INEQUALITY_CONST( v1, v2, out, success ) \
  { \
    (out) << #v1" = "<<Teuchos::toString(v1)<<" != "<<Teuchos::toString(v2)<<" : "; \
    const bool l_result = (v1) != (v2); \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    if (!l_result) (success) = false; \
  }


/** \brief Test that two values are not equal.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_INEQUALITY( v1, v2, out, success ) \
  { \
    (out) << #v1" = "<<Teuchos::toString(v1)<<" != "#v2" = "<<Teuchos::toString(v2)<<" : "; \
    const bool l_result = (v1) != (v2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test if two floating point values are equal to a given tolerance.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_FLOATING_EQUALITY( v1, v2, tol, out, success ) \
  { \
    const bool l_result = Teuchos::testRelErr( \
      #v1, v1, #v2, v2, "tol", tol, "tol", tol, Teuchos::outArg(out) ); \
    if (!l_result) (success) = false; \
  }


/** \brief Test if two iterators are equal.
 *
 * This macro does not try to print the iterators so it is more portable (in
 * terms of types).
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success ) \
  { \
    (out) << #iter1" == "#iter2" =  : "; \
    const bool l_result = (iter1) == (iter2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test if two iterators are NOT equal.
 *
 * This macro does not try to print the iterators so it is more portable (in
 * terms of types).
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_ITER_INEQUALITY( iter1, iter2, out, success ) \
  { \
    (out) << #iter1" != "#iter2" =  : "; \
    const bool l_result = (iter1) != (iter2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test that an array element value is equal to a given constant.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, printPass, out, success ) \
  { \
    const bool l_result = ( (a)[i] == (val) ); \
    if (!l_result) (success) = false; \
    if (printPass || !(l_result)) { \
      out << #a"["<<i<<"] = " << Teuchos::toString((a)[i]) << " == "#val" = " << Teuchos::toString(val) \
          << " : " << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    } \
  }


/** \brief Test that an array element value is not equal to a given constant.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_ARRAY_ELE_INEQUALITY( a, i, val, printPass, out, success ) \
  { \
    const bool l_result = ( (a)[i] != (val) ); \
    if (!l_result) (success) = false; \
    if (printPass || !(l_result)) { \
      out << #a"["<<i<<"] = " << Teuchos::toString((a)[i]) << " != "#val" = " << Teuchos::toString(val) \
          << " : " << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    } \
  }


/** \brief Test if a floating-point array element value is equal to a given
 * constant for a given tolerance.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_MATRIX_ELE_FLOATING_EQUALITY( a, i, j, val, tol, printPass, out, success ) \
  { \
    std::ostringstream a_i_str; \
    a_i_str <<#a<<"("<<i<<","<<j<<")"; \
    const bool l_result = Teuchos::testRelErr( \
      a_i_str.str(), (a)(i,j), #val, val, "tol", tol, "tol", tol, \
      (printPass) ? Teuchos::outArg(out) : Teuchos::null ); \
    if (!l_result) (success) = false; \
  }


/** \brief Test if a matrix element value is equal to a given constant.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_MATRIX_ELE_EQUALITY( a, i, j, val, printPass, out, success ) \
  { \
    const bool l_result = ( (a)(i,j) == (val) ); \
    if (!l_result) (success) = false; \
    if (printPass || !(l_result)) { \
      out << #a"("<<i<<","<<j<<") = " << (a)(i,j) << " == "#val" = " << (val) \
          << " : " << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
    } \
  }


/** \brief Compare two objects using an input comparion operator.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success ) \
  { \
    out << #v1" = "<<(v1)<<" "#comp" "#v2" = "<<(v2)<<" : "; \
    const bool l_result = (v1) comp (v2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test that a chunk of code throws an expected exception.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_THROW( code, ExceptType, out, success  ) \
  try { \
    (out) << "Test that code {"#code";} throws " \
          <<Teuchos::TypeNameTraits<ExceptType>::name()<<": "; \
    code; \
    (success) = false; \
    (out) << "failed\n"; \
  } \
  catch (const ExceptType& except) { \
    out << "passed\n"; \
    out << "\nException message for expected exception:\n\n"; \
    { \
      Teuchos::OSTab l_tab(out); \
      out << except.what() << "\n\n"; \
    } \
  }


/** \brief Test that a chunk of code does not throw any exceptions.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_NOTHROW( code, out, success  ) \
  try { \
    (out) << "Test that code {"#code";} does not throw : "; \
    code; \
    (out) << "passes\n"; \
  } \
  catch (...) { \
    (success) = false; \
    out << "failed\n"; \
  }


//
// Implementations
//


inline
void Teuchos::updateSuccess(const bool result, bool &success)
{
  if (!result) success = false;
}


inline
const std::string
Teuchos::passfail(const bool result)
{
  if (!result)
    return "FAILED";
  return "passed";
}


template <class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Teuchos::relErr( const Scalar &s1, const Scalar &s2 )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  return
    ST::magnitude( s1 - s2 )
    / (
      ST::magnitude(
        RelErrSmallNumber<ST::hasMachineParameters,Scalar>::smallNumber()
        )
      + std::max( ST::magnitude(s1), ST::magnitude(s2) )
      );
}


template<class Scalar>
bool Teuchos::testRelErr(
  const std::string &v1_name,
  const Scalar &v1,
  const std::string &v2_name,
  const Scalar &v2,
  const std::string &maxRelErr_error_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_error,
  const std::string &maxRelErr_warning_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_warning,
  const Ptr<std::ostream> &out
  )
{
  using std::endl;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const ScalarMag rel_err = relErr( v1, v2 );
  const bool success = ( !SMT::isnaninf(rel_err) && !SMT::isnaninf(maxRelErr_error)
    && rel_err <= maxRelErr_error );
  if (!is_null(out)) {
    *out
      << endl
      << "Check: rel_err(" << v1_name << ", " << v2_name << ")\n"
      << "       = rel_err(" << v1 << ", " << v2 << ") "
      << "= " << rel_err << endl
      << "         <= " << maxRelErr_error_name
      << " = " << maxRelErr_error << " : " << passfail(success) << endl;
    if( success && rel_err >= maxRelErr_warning ) {
      *out
        << "Warning! rel_err(" << v1_name << ", " << v2_name << ")\n"
        << "       = rel_err(" << v1 << ", " << v2 << ") "
        << "= " << rel_err << endl
        << "         >= " << maxRelErr_warning_name
        << " = " << maxRelErr_warning << "!\n";
    }
  }
  return success;
}


template<class Array1, class Array2>
bool Teuchos::compareArrays(
  const Array1 &a1, const std::string &a1_name,
  const Array2 &a2, const std::string &a2_name,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as;
  bool success = true;

  out << "Comparing " << a1_name << " == " << a2_name << " ... ";

  const int n = a1.size();

  // Compare sizes
  if (as<int>(a2.size()) != n) {
    out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == " 
        << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
    return false;
  }
  
  // Compare elements
  for( int i = 0; i < n; ++i ) {
    const bool result = ( a1[i] == a2[i] ); // Tests C::operator[](i) const
    if (!result) {
      out << "\nError, "<<a1_name<<"["<<i<<"] = "<<a1[i]<<" == "
          << a2_name<<"["<<i<<"] = "<<a2[i]<<": failed!\n";
      success = false;
    }
  }
  if (success) {
    out << "passed\n";
  }

  return success;

}


template<class Array1, class Array2, class ScalarMag>
bool Teuchos::compareFloatingArrays(
  const Array1 &a1, const std::string &a1_name,
  const Array2 &a2, const std::string &a2_name,
  const ScalarMag &tol,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as;
  bool success = true;

  out << "Comparing " << a1_name << " == " << a2_name << " ... ";

  const int n = a1.size();

  // Compare sizes
  if (as<int>(a2.size()) != n) {
    out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == " 
        << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
    return false;
  }
  
  // Compare elements
  for( int i = 0; i < n; ++i ) {
    const ScalarMag err = relErr( a1[i], a2[i] );
    if ( err > tol ) {
      out
        <<"\nError, relErr("<<a1_name<<"["<<i<<"],"
        <<a2_name<<"["<<i<<"]) = relErr("<<a1[i]<<","<<a2[i]<<") = "
        <<err<<" <= tol = "<<tol<<": failed!\n";
      success = false;
    }
  }
  if (success) {
    out << "passed\n";
  }

  return success;

}


#endif  // TEUCHOS_TESTING_HELPERS_HPP
