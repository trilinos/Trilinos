#include <Tpetra_config.h>
#include <Kokkos_DefaultNode.hpp>

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

//
// If Tpetra is compiled with explicit instantiation,
// then we will use these types in our tests.
//
// Epetra uses int, int, double data types.  If we
// are using these data types, then we can test 
// cases of Epetra user input.
//


#ifdef HAVE_ZOLTAN2_INST_FLOAT_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef float scalar_t;

#elif HAVE_ZOLTAN2_INST_DOUBLE_INT_LONG

typedef int lno_t;
typedef long gno_t;
typedef double scalar_t;

#elif HAVE_ZOLTAN2_INST_FLOAT_INT_INT

typedef int lno_t;
typedef int gno_t;
typedef float scalar_t;

#elif HAVE_ZOLTAN2_INST_DOUBLE_INT_INT

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

