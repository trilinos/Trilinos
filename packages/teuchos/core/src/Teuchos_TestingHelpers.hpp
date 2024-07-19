// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
template <class Scalar1, class Scalar2>
typename Teuchos::ScalarTraits< typename std::common_type<Scalar1,Scalar2>::type >::magnitudeType
relErr( const Scalar1 &s1, const Scalar2 &s2 );


/** \brief Compute, check and optionally print the relative error in two scalars.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_test_tools_code_grp
 */
template <typename T1, typename T2, typename Enabled = void>
struct TestRelErr {
  typedef typename Teuchos::ScalarTraits<T1>::magnitudeType magType1;
  typedef typename Teuchos::ScalarTraits<T2>::magnitudeType magType2;
  typedef typename std::common_type<magType1,magType2>::type magnitudeType;
  static bool eval(
    const std::string &v1_name,
    const T1 &v1,
    const std::string &v2_name,
    const T2 &v2,
    const std::string &maxRelErr_error_name,
    const magnitudeType &maxRelErr_error,
    const std::string &maxRelErr_warning_name,
    const magnitudeType &maxRelErr_warning,
    const Ptr<std::ostream> &out
  )
  {
    using std::endl;
    typedef ScalarTraits<magnitudeType> SMT;
    const magnitudeType rel_err = relErr( v1, v2 );
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
};

template <typename T1, typename T2>
bool testRelErr(
    const std::string &v1_name,
    const T1 &v1,
    const std::string &v2_name,
    const T2 &v2,
    const std::string &maxRelErr_error_name,
    const typename TestRelErr<T1,T2>::magnitudeType &maxRelErr_error,
    const std::string &maxRelErr_warning_name,
    const typename TestRelErr<T1,T2>::magnitudeType &maxRelErr_warning,
    const Ptr<std::ostream> &out) {
  return TestRelErr<T1,T2>::eval(v1_name, v1, v2_name,v2, maxRelErr_error_name, maxRelErr_error, maxRelErr_warning_name, maxRelErr_warning, out);
}

/** \brief Compare if two array objects are the same or not.
 *
 * This function works with any two array-like objects:
 * <tt>Array1</tt> and <tt>Array2</tt> must have <tt>size()</tt>
 * and <tt>operator[](i)</tt> methods defined. Their element types
 * may be different, but it must be possible to compare them with
 * <tt>==</tt>.
 *
 * \returns Returns <tt>true</tt> if <tt>a1.size() == a2.size()</tt> and
 *   their elements are equal. Otherwise returns <tt>false</tt>.
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
 * This function works with any two array-like objects:
 * <tt>Array1</tt> and <tt>Array2</tt> must have <tt>size()</tt>
 * and <tt>operator[](i)</tt> methods defined. Their element
 * types may be different, as long as <tt>Teuchos::relErr(a1[i], a2[i])</tt>
 * can be called to determine the relative difference between them.
 *
 * \returns Returns <tt>true</tt> if <tt>a1.size() == a2.size()</tt> and
 *   <tt>relErr(a1[i], a2[i]) <= tol</tt> for all <tt>0 <= i < a1.size()</tt>.
 *   Otherwise returns <tt>false</tt>.
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

/** \brief Compare if two array objects are the same up to an absolute
 * tolerance: elements <tt>a1[i]</tt> and <tt>a2[i]</tt> are considered
 * the same if <tt>|a1[i]-a2[i]| <= tol</tt>.
 *
 * This function works with any two array-like objects
 * with the same element types. <tt>Array1</tt> and <tt>Array2</tt> must have
 * <tt>size()</tt> and <tt>operator[](i)</tt> methods defined.
 * <tt>Teuchos::ScalarTraits</tt> must also have a specialization
 * for the element type.
 *
 * \returns Returns <tt>true</tt> if <tt>a1.size() == a2.size()</tt> and
 *   <tt>|a1[i]-a2[i]| <= tol</tt> for all <tt>0 <= i < a1.size()</tt>.
 *   Otherwise returns <tt>false</tt>.
 *
 * \ingroup teuchos_testing_grp
 */
template<class Array1, class Array2, class ScalarMag>
bool compareFloatingArraysAbsolute(
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


/** \brief Compare two objects using an input comparsion operator.
 *
 * The test succeeds (passes) if and only if "(v1) comp (v2)".
 * For example, TEUCHOS_TEST_COMPARE( 2, <, 3, out, success )
 * succeeds, but TEUCHOS_TEST_COMPARE( 2, >, 3, out, success )
 * and TEUCHOS_TEST_COMPARE( 3, <, 2, out, success ) both fail.
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


/** \brief Compare an object and a constant using an input comparsion
 * operator.
 *
 * The test succeeds (passes) if and only if "(v1) comp (v2)".
 * For example, TEUCHOS_TEST_COMPARE( 2, <, 3, out, success )
 * succeeds, but TEUCHOS_TEST_COMPARE( 2, >, 3, out, success )
 * and TEUCHOS_TEST_COMPARE( 3, <, 2, out, success ) both fail.
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_COMPARE_CONST( v1, comp, v2, out, success ) \
  { \
    out << #v1" = "<<(v1)<<" "#comp" "<<(v2)<<" : "; \
    const bool l_result = (v1) comp (v2); \
    if (!l_result) (success) = false; \
    (out) << TEUCHOS_PASS_FAIL(l_result) << "\n"; \
  }


/** \brief Test that the chunk of code 'code' throws an expected exception.
 *
 * 'code' is a chunk of code to execute.  It will be executed exactly
 * once.  If it throws an exception of type ExceptType, this test
 * passes (and prints "passed").  Otherwise, it prints "failed" with
 * an informative message.  The macro prints all messages to the given
 * output stream (std::ostream&) out.  Furthermore, if the test
 * passes, it assigns true to success; if the test fails, it assigns
 * false to success.
 *
 * The macro's implementation does not evaluate 'out' more than once.
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_THROW( code, ExceptType, out, success  ) \
  { \
    std::ostream& l_out = (out); \
    try { \
      l_out << "Test that code {"#code";} throws " \
            << Teuchos::TypeNameTraits<ExceptType>::name () << ": "; \
      code; \
      (success) = false; \
      l_out << "failed (code did not throw an exception at all)\n"; \
    } \
    catch (const ExceptType& except) { \
      l_out << "passed\n";                                        \
      l_out << "\nException message for expected exception:\n\n";   \
      { \
        Teuchos::OSTab l_tab(out); \
        l_out << except.what () << "\n\n"; \
        (void)l_tab; \
      } \
    } \
    catch (std::exception& except) { \
      l_out << "The code was supposed to throw an exception of type "   \
            << Teuchos::TypeNameTraits<ExceptType>::name () << ", but " \
            << "instead threw an exception of type " \
            << typeid (except).name () << ", which is a subclass of " \
            << "std::exception.  The exception's message is:\n\n"; \
      { \
        Teuchos::OSTab l_tab(out); \
        l_out << except.what () << "\n\n"; \
        (void)l_tab; \
      } \
      l_out << "failed\n"; \
    } \
    catch (...) { \
      l_out << "The code was supposed to throw an exception of type "   \
            << Teuchos::TypeNameTraits<ExceptType>::name () << ", but " \
            << "instead threw an exception of some unknown type, which is " \
            << "not a subclass of std::exception.  This means we cannot " \
            << "show you the exception's message, if it even has one.\n\n"; \
      l_out << "failed\n"; \
    } \
  }


/** \brief Test that a chunk of code does not throw any exceptions.
 *
 * This macro is not complicated so take a look for yourself!
 *
 * \ingroup teuchos_testing_grp
 */
#define TEUCHOS_TEST_NOTHROW( code, out, success  ) \
  { \
    std::ostream& l_out = (out); \
    try { \
      l_out << "Test that code {"#code";} does not throw : "; \
      code; \
      l_out << "passed\n"; \
    } \
    catch (std::exception& except) { \
      (success) = false; \
      l_out << "The code was not supposed to throw an exception, but " \
            << "instead threw an exception of type " \
            << typeid (except).name () << ", which is a subclass of " \
            << "std::exception.  The exception's message is:\n\n"; \
      { \
        Teuchos::OSTab l_tab(out); \
        l_out << except.what () << "\n\n"; \
        (void)l_tab; \
      } \
      l_out << "failed\n"; \
    } \
    catch (...) { \
      (success) = false; \
      l_out << "The code was not supposed to throw an exception, but " \
            << "instead threw an exception of some unknown type, which is " \
            << "not a subclass of std::exception.  This means we cannot " \
            << "show you the exception's message, if it even has one.\n\n"; \
      l_out << "failed\n"; \
    } \
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


template <class Scalar1, class Scalar2>
typename Teuchos::ScalarTraits< typename std::common_type<Scalar1,Scalar2>::type >::magnitudeType
Teuchos::relErr( const Scalar1 &s1, const Scalar2 &s2 )
{
  typedef typename std::common_type<Scalar1,Scalar2>::type Scalar;
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

  // Determine the element types of Array1 and Array2
  using Elem1 = std::decay_t<decltype(std::declval<Array1>()[0])>;
  using Elem2 = std::decay_t<decltype(std::declval<Array2>()[0])>;

  static_assert(std::is_same_v<Elem1, Elem2>,
      "Teuchos::compareFloatingArrays: element types of Array1 and Array2 must be the same.");

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
    if ( !(err <= tol) ) {
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

template<class Array1, class Array2, class ScalarMag>
bool Teuchos::compareFloatingArraysAbsolute(
  const Array1 &a1, const std::string &a1_name,
  const Array2 &a2, const std::string &a2_name,
  const ScalarMag &tol,
  Teuchos::FancyOStream &out
  )
{
  using Teuchos::as;
  using Teuchos::ScalarTraits;

  // Determine the element types of Array1 and Array2
  using Elem1 = std::decay_t<decltype(std::declval<Array1>()[0])>;
  using Elem2 = std::decay_t<decltype(std::declval<Array2>()[0])>;

  static_assert(std::is_same_v<Elem1, Elem2>,
      "Teuchos::compareFloatingArraysAbsolute: element types of Array1 and Array2 must be the same.");

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
    const ScalarMag err = ScalarTraits<Elem1>::magnitude( a1[i] - a2[i] );
    if ( !(err <= tol) ) {
      out
        <<"\nError, ||"<<a1_name<<"["<<i<<"] - " << a2_name<<"["<<i<<"]|| = "
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
