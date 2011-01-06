/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_percept_unit_tests_stk_utest_macros_hpp
#define stk_percept_unit_tests_stk_utest_macros_hpp

#ifndef STK_BUILT_IN_SIERRA
#include <STK_config.h>
#endif

//
//This file is kind of like a unit-test abstraction layer:
//A series of STKUNIT_* macros are defined in terms of either
//cppunit macros, or trilinos/teuchos unit-test macros, depending
//on whether stk_mesh is being built as a sierra product or Trilinos package.
//

#ifdef HAVE_STK_Trilinos
//If we're building as a Trilinos package, then we'll use the Teuchos unit-test macros.

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#define STKUNIT_UNIT_TEST(testclass,testmethod) TEUCHOS_UNIT_TEST(testclass,testmethod)

#define STKUNIT_ASSERT(A) \
    {bool success = true; TEUCHOS_TEST_ASSERT(A,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_EQUAL(A,B) \
    {bool success = true; TEUCHOS_TEST_EQUALITY(B,A,std::cout,success); if (!success) throw 1;}
#define STKUNIT_EXPECT_EQUAL(A,B) STKUNIT_ASSERT_EQUAL(A,B)
#define STKUNIT_ASSERT_THROW(A,B) \
    {bool success = true; TEUCHOS_TEST_THROW(A,B,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_NO_THROW(A) \
    {TEUCHOS_TEST_NOTHROW(A,out,success)}
#define STKUNIT_EXPECT_TRUE(A) \
    {bool success = true; TEUCHOS_TEST_ASSERT(A,std::cout,success); if (!success) throw 1;}
#define STKUNIT_EXPECT_FALSE(A) \
    {bool success = true; TEUCHOS_TEST_ASSERT(!(A),std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_TRUE(A) STKUNIT_EXPECT_TRUE(A)
#define STKUNIT_ASSERT_FALSE(A) STKUNIT_EXPECT_FALSE(A)


#define STKUNIT_MAIN(argc,argv) \
int main(int argc,char**argv) {\
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); \
  bool result = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv); \
  std::cout << std::endl ; \
  return result; \
}


#elif HAVE_STK_CPPUNIT

//If we're not in Trilinos, we're a sierra product, so we're using cppunit:

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#define STKUNIT_ASSERT(A) CPPUNIT_ASSERT(A)
#define STKUNIT_ASSERT_EQUAL(A,B) CPPUNIT_ASSERT_EQUAL(A,B)
#define STKUNIT_EXPECT_EQUAL(A,B) CPPUNIT_ASSERT_EQUAL(A,B)
#define STKUNIT_ASSERT_THROW(A,B) CPPUNIT_ASSERT_THROW(A,B)

#define STKUNIT_UNIT_TEST(testclass,testmethod) \
class testclass##_##testmethod : public CppUnit::TestCase { \
private: \
  CPPUNIT_TEST_SUITE( testclass##_##testmethod ); \
  CPPUNIT_TEST( testmethod ); \
  CPPUNIT_TEST_SUITE_END(); \
public: \
  testclass##_##testmethod() : CppUnit::TestCase() {} \
  void setUp() {} \
  void tearDown() {} \
  void testmethod(); \
}; \
CPPUNIT_TEST_SUITE_REGISTRATION( testclass##_##testmethod ); \
void testclass##_##testmethod::testmethod() \

#define STKUNIT_MAIN(argc,argv) \
int main(int argc, char ** argv) \
{ \
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) { \
    std::cerr << "MPI_Init FAILED" << std::endl ; \
    std::abort(); \
  } \
 \
  std::string test_names; \
 \
  { \
    for (int i = 0; i < argc; ++i) { \
      if (std::string(argv[i]) == "-test") { \
        if (i + 1 < argc) \
          test_names = argv[i + 1]; \
        else \
          std::cout << "Running all tests" << std::endl; \
      } \
    } \
 \
    CppUnit::TextUi::TestRunner runner; \
    CppUnit::TestFactoryRegistry & registry = \
      CppUnit::TestFactoryRegistry::getRegistry(); \
    runner.addTest( registry.makeTest() ); \
    runner.run(test_names); \
  } \
 \
  MPI_Finalize(); \
 \
  return 0; \
}

#else // HAVE_STK_GTEST

#include <gtest/gtest.h>

#define STKUNIT_ASSERT(A) ASSERT_TRUE(A)
#define STKUNIT_ASSERT_EQUAL(A,B) ASSERT_EQ(A,B)
#define STKUNIT_EXPECT_EQUAL(A,B) EXPECT_EQ(A,B)
#define STKUNIT_ASSERT_THROW(A,B) ASSERT_THROW(A,B)
#define STKUNIT_ASSERT_NO_THROW(A) ASSERT_NO_THROW(A)
#define STKUNIT_EXPECT_TRUE(A) EXPECT_TRUE(A)
#define STKUNIT_EXPECT_FALSE(A) EXPECT_FALSE(A)
#define STKUNIT_ASSERT_TRUE(A) ASSERT_TRUE(A)
#define STKUNIT_ASSERT_FALSE(A) ASSERT_FALSE(A)

#define STKUNIT_UNIT_TEST(testclass,testmethod) TEST(testclass,testmethod)

#define STKUNIT_MAIN(argc,argv) \
int main(int argc, char **argv) { \
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) { \
    std::cerr << "MPI_Init FAILED" << std::endl ; \
    std::abort(); \
  } \
  std::cout << "Running main() from gtest_main.cc\n"; \
  testing::InitGoogleTest(&argc, argv); \
  bool result = RUN_ALL_TESTS(); \
  MPI_Finalize(); \
  return result; \
}


#endif // HAVE_STK_Trilinos

#endif

