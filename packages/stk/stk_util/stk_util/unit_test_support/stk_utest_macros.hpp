/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_unit_tests_stk_utest_macros_hpp
#define stk_mesh_unit_tests_stk_utest_macros_hpp

#ifndef STK_BUILT_IN_SIERRA
#include <STK_config.h>
#else
#define HAVE_MPI
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

//
//This file is kind of like a unit-test abstraction layer:
//A series of STKUNIT_* macros are defined in terms of either
//gtest macros, or trilinos/teuchos unit-test macros, depending
//on whether stk_mesh is being built as a sierra product or Trilinos package.
//
#ifdef HAVE_MPI
#define RUN_TEST_REDUCE(error) \
  int tmp_error = error; \
  int reduce_result = MPI_Allreduce ( &tmp_error, &error, 1 /*count*/, \
                                      MPI_INT, MPI_MAX, MPI_COMM_WORLD ); \
  if (reduce_result != MPI_SUCCESS) { \
    std::cerr << "MPI_Allreduce FAILED" << std::endl; \
    error = true; \
  } \
  if ( !error ) \
    std::cout << "STKUNIT_ALL_PASS" << std::endl;
#else
#define RUN_TEST_REDUCE(error) \
  if (error != 0) { \
    std::cerr << "Test FAILED" << std::endl; \
    error = true; \
  } \
  if ( !error ) \
    std::cout << "STKUNIT_ALL_PASS" << std::endl;
#endif


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
#define STKUNIT_ASSERT_EQ(A,B) STKUNIT_ASSERT_EQUAL(A,B)
#define STKUNIT_ASSERT_NE(A,B) \
    {bool success = true; TEUCHOS_TEST_INEQUALITY(B,A,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_LE(A,B) \
    {bool success = true; TEUCHOS_TEST_COMPARE(A,<=,B,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_LT(A,B) \
    {bool success = true; TEUCHOS_TEST_COMPARE(A,<,B,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_GE(A,B) \
    {bool success = true; TEUCHOS_TEST_COMPARE(A,>=,B,std::cout,success); if (!success) throw 1;}
#define STKUNIT_ASSERT_GT(A,B) \
    {bool success = true; TEUCHOS_TEST_COMPARE(A,>,B,std::cout,success); if (!success) throw 1;}
#define STKUNIT_EXPECT_EQ(A,B) STKUNIT_ASSERT_EQ(A,B)
#define STKUNIT_EXPECT_NE(A,B) STKUNIT_ASSERT_NE(A,B)
#define STKUNIT_EXPECT_LE(A,B) STKUNIT_ASSERT_LE(A,B)
#define STKUNIT_EXPECT_LT(A,B) STKUNIT_ASSERT_LT(A,B)
#define STKUNIT_EXPECT_GE(A,B) STKUNIT_ASSERT_GE(A,B)
#define STKUNIT_EXPECT_GT(A,B) STKUNIT_ASSERT_GT(A,B)
#define STKUNIT_ASSERT_STREQ(A,B) STKUNIT_ASSERT_EQ(A,B)
#define STKUNIT_ASSERT_STRNE(A,B) STKUNIT_ASSERT_NE(A,B)
#define STKUNIT_ASSERT_STRCASEEQ(A,B) STKUNIT_ASSERT_STREQ(toupper(A),toupper(B))
#define STKUNIT_ASSERT_STRCASENE(A,B) STKUNIT_ASSERT_STRNE(toupper(A),toupper(B))
#define STKUNIT_EXPECT_STREQ(A,B) STKUNIT_EXPECT_EQ(A,B)
#define STKUNIT_EXPECT_STRNE(A,B) STKUNIT_EXPECT_NE(A,B)
#define STKUNIT_EXPECT_STRCASEEQ(A,B) STKUNIT_EXPECT_STREQ(toupper(A),toupper(B))
#define STKUNIT_EXPECT_STRCASENE(A,B) STKUNIT_EXPECT_STRNE(toupper(A),toupper(B))
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
#define STKUNIT_ASSERT_DOUBLE_EQ(A,B) \
    {bool success = true; TEUCHOS_TEST_EQUALITY(B,A,std::cout,success); if (!success) throw 1;}
#define STKUNIT_EXPECT_DOUBLE_EQ(A,B) STKUNIT_ASSERT_DOUBLE_EQ(A,B)
#define STKUNIT_ASSERT_NEAR(A,B, tol) \
    {bool success = true; TEUCHOS_TEST_FLOATING_EQUALITY(B,A,tol,std::cout,success); if (!success) throw 1;}
#define STKUNIT_EXPECT_NEAR(A,B, tol) STKUNIT_ASSERT_NEAR(A, B, tol)

#define STKUNIT_MAIN(argc,argv) \
int* STKUNIT_ARGC;                                                 \
char** STKUNIT_ARGV;                                               \
int main(int argc,char**argv) {\
STKUNIT_ARGC = &argc; \
STKUNIT_ARGV = argv; \
  Teuchos::GlobalMPISession mpiSession(&argc, &argv); \
  int error = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv); \
  RUN_TEST_REDUCE(error); \
  return error; \
}

#else // HAVE_STK_GTEST

#include <gtest/gtest.h>

#define STKUNIT_ASSERT(A) ASSERT_TRUE(A)
#define STKUNIT_ASSERT_EQUAL(A,B) ASSERT_EQ(A,B)
#define STKUNIT_EXPECT_EQUAL(A,B) EXPECT_EQ(A,B)
#define STKUNIT_ASSERT_EQ(A,B) ASSERT_EQ(A,B)
#define STKUNIT_ASSERT_NE(A,B) ASSERT_NE(A,B)
#define STKUNIT_ASSERT_LE(A,B) ASSERT_LE(A,B)
#define STKUNIT_ASSERT_LT(A,B) ASSERT_LT(A,B)
#define STKUNIT_ASSERT_GE(A,B) ASSERT_GE(A,B)
#define STKUNIT_ASSERT_GT(A,B) ASSERT_GT(A,B)
#define STKUNIT_EXPECT_EQ(A,B) EXPECT_EQ(A,B)
#define STKUNIT_EXPECT_NE(A,B) EXPECT_NE(A,B)
#define STKUNIT_EXPECT_LE(A,B) EXPECT_LE(A,B)
#define STKUNIT_EXPECT_LT(A,B) EXPECT_LT(A,B)
#define STKUNIT_EXPECT_GE(A,B) EXPECT_GE(A,B)
#define STKUNIT_EXPECT_GT(A,B) EXPECT_GT(A,B)
#define STKUNIT_ASSERT_STREQ(A,B) ASSERT_STREQ(A,B)
#define STKUNIT_ASSERT_STRNE(A,B) ASSERT_STRNE(A,B)
#define STKUNIT_ASSERT_STRCASEEQ(A,B) ASSERT_STRCASEEQ(A,B)
#define STKUNIT_ASSERT_STRCASENE(A,B) ASSERT_STRCASENE(A,B)
#define STKUNIT_EXPECT_STREQ(A,B) EXPECT_STREQ(A,B)
#define STKUNIT_EXPECT_STRNE(A,B) EXPECT_STRNE(A,B)
#define STKUNIT_EXPECT_STRCASEEQ(A,B) EXPECT_STRCASEEQ(A,B)
#define STKUNIT_EXPECT_STRCASENE(A,B) EXPECT_STRCASENE(A,B)
#define STKUNIT_ASSERT_THROW(A,B) ASSERT_THROW(A,B)
#define STKUNIT_ASSERT_NO_THROW(A) ASSERT_NO_THROW(A)
#define STKUNIT_EXPECT_TRUE(A) EXPECT_TRUE(A)
#define STKUNIT_EXPECT_FALSE(A) EXPECT_FALSE(A)
#define STKUNIT_ASSERT_TRUE(A) ASSERT_TRUE(A)
#define STKUNIT_ASSERT_FALSE(A) ASSERT_FALSE(A)
#define STKUNIT_ASSERT_DOUBLE_EQ(A,B) ASSERT_DOUBLE_EQ(A,B)
#define STKUNIT_EXPECT_DOUBLE_EQ(A,B) EXPECT_DOUBLE_EQ(A,B)
#define STKUNIT_ASSERT_NEAR(A,B,tol) ASSERT_NEAR(A,B,tol)
#define STKUNIT_EXPECT_NEAR(A,B,tol) EXPECT_NEAR(A,B,tol)

#define STKUNIT_UNIT_TEST(testclass,testmethod) TEST(testclass,testmethod)

#define STKUNIT_MAIN(argc,argv)                                    \
int* STKUNIT_ARGC;                                                 \
char** STKUNIT_ARGV;                                               \
int main(int argc, char **argv) {                                  \
  if ( MPI_SUCCESS != MPI_Init( & argc , & argv ) ) {              \
    std::cerr << "MPI_Init FAILED" << std::endl ;                  \
    std::abort();                                                  \
  }                                                                \
  STKUNIT_ARGC = &argc;                                            \
  STKUNIT_ARGV = argv;                                             \
  std::cout << "Running main() from gtest_main.cc\n";              \
  testing::InitGoogleTest(&argc, argv);                            \
  int error = RUN_ALL_TESTS();                                     \
  RUN_TEST_REDUCE(error);                                          \
  MPI_Finalize();                                                  \
  return error;                                                    \
}

#define STKUNIT_WITH_TRACING_MAIN(argc, argv)                      \
int* STKUNIT_ARGC;                                                 \
char** STKUNIT_ARGV;                                               \
int main(int argc, char **argv) {                                  \
  use_case::UseCaseEnvironment use_case_environment(&argc, &argv); \
  std::cout << "Running main() from gtest_main.cc\n";              \
  testing::InitGoogleTest(&argc, argv);                            \
  STKUNIT_ARGC = &argc;                                            \
  STKUNIT_ARGV = argv;                                             \
  int error = RUN_ALL_TESTS();                                     \
  RUN_TEST_REDUCE(error);                                          \
  return error;                                                    \
}

#define STKUNIT_WITH_SIERRA_MAIN(argc,argv,prod)    \
int main(int argc, char **argv) { \
  sierra::Env::set_input_file_required(false); \
  testing::InitGoogleTest(&argc, argv); \
  sierra::Env::Startup startup__(&argc, &argv, sierra::prod::get_product_name(), __DATE__ " " __TIME__); \
  int error = RUN_ALL_TESTS();                                          \
  RUN_TEST_REDUCE(error);                                          \
  return error; \
}

#endif // HAVE_STK_Trilinos

#endif // stk_mesh_unit_tests_stk_utest_macros_hpp

