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

#ifndef TEUCHOS_UNIT_TEST_HELPERS_HPP
#define TEUCHOS_UNIT_TEST_HELPERS_HPP


/*! \file Teuchos_UnitTestHelpers.hpp
\brief Macros for defining unit tests.

The macros in this file are for naming and defining unit tests.  They
give your unit test a group and name, so that you can identify it in
the test output.  You are responsible for filling in the actual test.

For macros (like TEST_NOTHROW) to help you write the actual unit test,
see Teuchos_LocalTestingHelpers.hpp and Teuchos_TestingHelpers.hpp.
*/

#include "Teuchos_UnitTestBase.hpp"
#include "Teuchos_StaticSetupMacro.hpp"


/** \brief Macro for defining a (non-templated) unit test.
 *
 * The macro parameters TEST_GROUP and TEST_NAME must each be a valid
 * part of a C++ identifier.  In particular, they should not contain
 * any spaces or other characters that are not allowed in the name of
 * a class.
 *
 * Here is a brief example of how to declare a unit test using this macro.
 * \code
 * TEUCHOS_UNIT_TEST( myTestGroup, testDouble ) {
 *   double x, y;
 *   // Make sure that assigning 42 to y doesn't throw an exception.
 *   TEST_NOTHROW( x = 42 );
 *   // Make sure that assigning 42 to y doesn't throw an exception.
 *   TEST_NOTHROW( y = 42 );
 *   // Make sure that x and y are now equal.
 *   TEST_EQUALITY_CONST( x, y );
 * }
 * \endcode
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST(TEST_GROUP, TEST_NAME) \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest() \
      : Teuchos::UnitTestBase( #TEST_GROUP, #TEST_NAME ) \
    {} \
    virtual void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  TEST_GROUP##_##TEST_NAME##_UnitTest \
    instance_##TEST_GROUP##_##TEST_NAME##_UnitTest; \
  \
        void TEST_GROUP##_##TEST_NAME##_UnitTest::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \


/** \brief Macro for defining a templated unit test with one template parameter.
 *
 * Use this macro to <i>define</i> the templated unit test.  To
 * <i>instantiate</i> the unit test for a particular template
 * parameter, use the TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT macro.

 * \note If you don't instantiate the unit test, it won't be run.  In
 *   fact, in that case, the unit test will only be compiled as a
 *   template class.  This may not necessarily catch all compilation
 *   errors, since
 *   <a href="http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error">invalid
 *   substitution of template parameters is not an error</a>.
 *
 * The macro parameters TEST_GROUP and TEST_NAME must each be a valid
 * part of a C++ identifier.  In particular, they should not contain
 * any spaces or other characters that are not allowed in the name of
 * a class.
 *
 * On instantiation, the macro parameter TYPE must be the name of a
 * valid C++ type.  The type will be treated as a template parameter.
 * Thus, if your unit test references typedefs in TYPE, it should use
 * \c typename to access them.  For example:
 * \code
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( myTestGroup, myTestName, StlContainerType ) {
 *   typedef typename StlContainerType::value_type value_type;
 *   // ... the rest of the unit test ...
 * }
 * \endcode
 *
 * When instantiating the unit test for a particular template
 * parameter using TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT, the macro
 * parameter TYPE must not contain any spaces (Bug 5757).  If you want
 * TYPE to be a type with spaces, like <tt>unsigned int</tt> or
 * <tt>long double</tt>, you must first use a typedef to alias the
 * original type to a name without spaces, and then use the name
 * without spaces.  For example:
 *
 * \code
 * // Define the templated unit test.
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( myTestGroup, myTestName, RealType ) {
 *   RealType x = 42, y = 42;
 *   TEST_EQUALITY_CONST( x, y );
 *   // ... the rest of the unit test ...
 * }
 *
 * //
 * // Instantiate the unit test for various values of RealType.
 * //
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( myTestGroup, myTestName, float )
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( myTestGroup, myTestName, double )
 *
 * // Typedef to work around Bug 5757 (TYPE values cannot have spaces).
 * typedef long double long_double_type;
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( myTestGroup, myTestName, long_double_type )
 * \endcode
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(TEST_GROUP, TEST_NAME, TYPE) \
  template<class TYPE> \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest(const std::string& typeName) \
      : Teuchos::UnitTestBase( std::string(#TEST_GROUP)+"_"+typeName, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE> \
        void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \

/** \brief Instantiate a templated unit test with one template parameter.
 *
 * Use this macro to <i>instantiate</i> for a particular template
 * parameter value the templated unit test defined by
 * TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL with the same TEST_GROUP and
 * TEST_NAME values.  For more details and examples, see the
 * documentation of TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, TYPE) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE>; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE> \
  instance_##TEST_GROUP##_##TYPE##_##TEST_NAME##_UnitTest(#TYPE);


#ifdef HAVE_TEUCHOS_FLOAT
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME)\
    TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, float)
#else
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME)
#endif

#define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME)\
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, double)

#if defined(HAVE_TEUCHOS_COMPLEX) && defined(HAVE_TEUCHOS_FLOAT)
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME)\
     typedef std::complex<float> ComplexFloat; \
     TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, ComplexFloat)
#else
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME)
#endif

#ifdef HAVE_TEUCHOS_COMPLEX
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)\
     typedef std::complex<double> ComplexDouble; \
     TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(TEST_GROUP, TEST_NAME, ComplexDouble)
#else
#  define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)
#endif


/** \brief Instantiate a whole group of tests for supported real Scalar
 * types.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_REAL_SCALAR_TYPES(TEST_GROUP, TEST_NAME)\
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME)


/** \brief Instantiate a whole group of tests for supported Scalar types.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(TEST_GROUP, TEST_NAME)\
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_FLOAT(TEST_GROUP, TEST_NAME) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_DOUBLE(TEST_GROUP, TEST_NAME) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_FLOAT(TEST_GROUP, TEST_NAME) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_COMPLEX_DOUBLE(TEST_GROUP, TEST_NAME)


/** \brief Macro for defining a templated unit test with two template parameters.
 *
 * Use this macro to <i>define</i> the templated unit test.  To
 * <i>instantiate</i> the unit test for particular template parameter
 * values, use the TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT macro.

 * \note If you don't instantiate the unit test, it won't be run.  In
 *   fact, in that case, the unit test will only be compiled as a
 *   template class.  This may not necessarily catch all compilation
 *   errors, since
 *   <a href="http://en.wikipedia.org/wiki/Substitution_failure_is_not_an_error">invalid
 *   substitution of template parameters is not an error</a>.
 *
 * The macro parameters TEST_GROUP and TEST_NAME must each be a valid
 * part of a C++ identifier.  In particular, they should not contain
 * any spaces or other characters that are not allowed in the name of
 * a class.
 *
 * On instantiation, the macro parameters TYPE1 and TYPE2 must be the
 * name of a valid C++ type.  The types will be treated as template
 * parameters.  Thus, if your unit test references typedefs in TYPE1
 * or TYPE2, it should use \c typename to access them.  For example:
 * \code
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( myTestGroup, myTestName, StlContainer1Type, StlContainer2Type ) {
 *   typedef typename StlContainer1Type::value_type value1_type;
 *   typedef typename StlContainer2Type::value_type value2_type;
 *   // ... the rest of the unit test ...
 * }
 * \endcode
 *
 * When instantiating the unit test for a particular template
 * parameter using TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT, the macro
 * parameter values TYPE1 and TYPE2 must not contain any spaces (Bug
 * 5757).  If you want TYPE1 or TYPE2 to be a type with spaces, like
 * <tt>unsigned int</tt> or <tt>long double</tt>, you must first use a
 * typedef to alias the original type to a name without spaces, and
 * then use the name without spaces.  For example:
 *
 * \code
 * // Define the templated unit test.
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( myTestGroup, myTestName, IntegerType, RealType ) {
 *   IntegerType x = 42;
 *   RealType y = 42;
 *   TEST_EQUALITY_CONST( x, static_cast<IntegerType> (y) );
 *   // ... the rest of the unit test ...
 * }
 *
 * //
 * // Instantiate the unit test for various values of IntegerType and RealType.
 * //
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( myTestGroup, myTestName, int, float )
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( myTestGroup, myTestName, long, double )
 *
 * // Typedef to work around Bug 5757 (TYPE values cannot have spaces).
 * typedef long double long_double_type;
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( myTestGroup, myTestName, int, long_double_type )
 * \endcode
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL(TEST_GROUP, TEST_NAME, TYPE1, TYPE2) \
  template<class TYPE1, class TYPE2> \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest( \
      const std::string& type1Name, \
      const std::string& type2Name \
       ) \
      :Teuchos::UnitTestBase( \
         std::string(#TEST_GROUP)+"_"+type1Name+"_"+type2Name, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE1, class TYPE2> \
        void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1,TYPE2>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \

/** \brief Instantiate a templated unit test with two template parameters.
 *
 * Use this macro to <i>instantiate</i> for particular template
 * parameter values the templated unit test defined by
 * TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL with the same TEST_GROUP and
 * TEST_NAME values.  For more details and examples, see the
 * documentation of TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2);


/** \brief Macro for defining a templated unit test with three template parameters.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3) \
  template<class TYPE1, class TYPE2, class TYPE3> \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest( \
      const std::string& type1Name, \
      const std::string& type2Name, \
      const std::string& type3Name \
       ) \
      :Teuchos::UnitTestBase( \
         std::string(#TEST_GROUP)+"_"+type1Name+"_"+type2Name+"_"+type3Name, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE1, class TYPE2, class TYPE3> \
        void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1,TYPE2,TYPE3>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \


/** \brief Instantiate a templated unit test with three template parameters.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3);


/** \brief Macro for defining a templated unit test with four template parameters.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4) \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4> \
  class TEST_GROUP##_##TEST_NAME##_UnitTest : public Teuchos::UnitTestBase \
        { \
  public: \
    TEST_GROUP##_##TEST_NAME##_UnitTest( \
      const std::string& type1Name, \
      const std::string& type2Name, \
      const std::string& type3Name, \
      const std::string& type4Name \
       ) \
      :Teuchos::UnitTestBase( \
         std::string(#TEST_GROUP)+"_"+type1Name+"_"+type2Name+"_"+type3Name+"_"+type4Name, #TEST_NAME ) \
    {} \
    void runUnitTestImpl( Teuchos::FancyOStream &out, bool &success ) const; \
    virtual std::string unitTestFile() const { return __FILE__; } \
    virtual long int unitTestFileLineNumber() const { return __LINE__; } \
  }; \
  \
  template<class TYPE1, class TYPE2, class TYPE3, class TYPE4> \
        void TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1,TYPE2,TYPE3,TYPE4>::runUnitTestImpl( \
    Teuchos::FancyOStream &out, bool &success ) const \


/** \brief Instantiate a templated unit test with four template parameters.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TYPE4##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3,#TYPE4);


#endif  // TEUCHOS_UNIT_TEST_HELPERS_HPP
