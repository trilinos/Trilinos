#include <Tpetra_config.h>
#include <Kokkos_DefaultNode.hpp>

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

// The path to the directory of test data

#define STR_VALUE(path) #path
#define PATH_NAME(path) STR_VALUE(path)

#ifdef Z2_DATA_DIR
  std::string testDataFilePath(PATH_NAME(Z2_DATA_DIR));
#else
  std::string testDataFilePath(".");
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

#include <ErrorHandlingForTests.hpp>  
#include <UserInputForTests.hpp>

