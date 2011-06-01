// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_UNIT_TEST_HELPERS_HPP
#define TEUCHOS_UNIT_TEST_HELPERS_HPP


/*! \file Teuchos_UnitTestHelpers.hpp

\brief Macros for helping to create concrete unit tests.
*/


#include "Teuchos_UnitTestBase.hpp"
#include "Teuchos_StaticSetupMacro.hpp"


/** \brief Basic unit test creation macro for non-templated code.
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


/** \brief Basic unit test creation macro for templated code on one template
 * parameter.
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


/** \brief Template instantiation for a single templated type.
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


/** \brief Basic unit test creation macro for templated code on two template
 * parameters.
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

/** \brief Template instantiation for two templated types.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2);


/** \brief Basic unit test creation macro for templated code on three template
 * parameters.
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


/** \brief Template instantiation for three templated types.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3);


/** \brief Basic unit test creation macro for templated code on four template
 * parameters.
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


/** \brief Template instantiation for four templated types.
 *
 * \ingroup Teuchos_UnitTestDefinitionMacros_grp
 */
#define TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(TEST_GROUP, TEST_NAME, TYPE1, TYPE2, TYPE3, TYPE4) \
  \
  template class TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4 >; \
  TEST_GROUP##_##TEST_NAME##_UnitTest<TYPE1, TYPE2, TYPE3, TYPE4 > \
  instance_##TEST_GROUP##_##TYPE1##_##TYPE2##_##TYPE3##_##TYPE4##_##TEST_NAME##_UnitTest(#TYPE1,#TYPE2,#TYPE3,#TYPE4);


#endif  // TEUCHOS_UNIT_TEST_HELPERS_HPP
