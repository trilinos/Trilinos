// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TestHelpers.hpp
    \brief common code used by tests
*/
#ifndef ZOLTAN2_TESTHELPERS_HPP
#define ZOLTAN2_TESTHELPERS_HPP

#include <Zoltan2_Util.hpp>
#include <iostream>

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

// The path to the directory of test data

#define STR_VALUE(path) #path
#define PATH_NAME(path) STR_VALUE(path)

#ifdef Z2_DATA_DIR
  std::string testDataFilePath(PATH_NAME(Z2_DATA_DIR));
#else
  std::string testDataFilePath(".");
#endif

// The path to the Zoltan1 test directory.  We use
// some of their data for testing.

#ifdef Z1_TEST_DIR
  std::string zoltanTestDirectory(PATH_NAME(Z1_TEST_DIR));
#else
  std::string zoltanTestDirectory(".");
#endif

//
// If Tpetra is compiled with explicit instantiation,
// then we will use these types in our tests.
//
// Epetra uses int, int, double data types.  If we
// are using these data types, then we can test 
// cases of Epetra user input.
//


#if defined HAVE_ZOLTAN2_INST_FLOAT_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef float scalar_t;

#elif defined HAVE_ZOLTAN2_INST_DOUBLE_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef double scalar_t;

#elif defined HAVE_ZOLTAN2_INST_FLOAT_INT_INT

typedef int lno_t;
typedef int gno_t;
typedef float scalar_t;

#elif defined HAVE_ZOLTAN2_INST_DOUBLE_INT_INT

typedef int lno_t;
typedef int gno_t;
typedef double scalar_t;
#define HAVE_EPETRA_DATA_TYPES

#else

typedef int lno_t;
typedef int gno_t;
typedef double scalar_t;
#define HAVE_EPETRA_DATA_TYPES

#endif

#define MEMORY_CHECK(iPrint, msg) \
  if (iPrint){ \
    long kb = Zoltan2::getProcessKilobytes(); \
    std::cout.width(10); \
    std::cout.fill('*'); \
    std::cout << kb << " KB, " << msg << std::endl; \
    std::cout.width(0); \
    std::cout.fill(' '); \
  }

#include <ErrorHandlingForTests.hpp>  
#include <UserInputForTests.hpp>
#include <PrintData.hpp>

#endif
