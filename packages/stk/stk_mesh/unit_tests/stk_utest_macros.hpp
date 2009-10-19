#ifndef stk_mesh_unit_tests_stk_utest_macros_hpp
#define stk_mesh_unit_tests_stk_utest_macros_hpp

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

#define STKUNIT_ASSERT(A) \
    {bool success; TEUCHOS_TEST_ASSERT(A,std::cout,success)}
#define STKUNIT_ASSERT_EQUAL(A,B) \
    {bool success; TEUCHOS_TEST_EQUALITY(A,B,std::cout,success)}
#define STKUNIT_ASSERT_THROW(A,B) \
    {bool success; TEUCHOS_TEST_THROW(A,B,std::cout,success)}

#define STKUNIT_MAIN(argc,argv) \
int main(int argc,char**argv) {\
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); \
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv); \
}

#define STKUNIT_UNIT_TEST(testclass,testmethod) TEUCHOS_UNIT_TEST(testclass,testmethod)

#else
//If we're not in Trilinos, we're a sierra product, so we're using cppunit:

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#define STKUNIT_ASSERT(A) CPPUNIT_ASSERT(A)
#define STKUNIT_ASSERT_EQUAL(A,B) CPPUNIT_ASSERT_EQUAL(A,B)
#define STKUNIT_ASSERT_THROW(A,B) CPPUNIT_ASSERT_THROW(A,B)

#define STKUNIT_UNIT_TEST(testclass,testmethod) \
class testclass : public CppUnit::TestCase { \
private: \
  CPPUNIT_TEST_SUITE( testclass ); \
  CPPUNIT_TEST( testmethod ); \
  CPPUNIT_TEST_SUITE_END(); \
public: \
  testclass() : CppUnit::TestCase() {} \
  void setUp() {} \
  void tearDown() {} \
  void testmethod(); \
}; \
CPPUNIT_TEST_SUITE_REGISTRATION( testclass ); \
void testclass::testmethod() \

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

#endif

#endif

