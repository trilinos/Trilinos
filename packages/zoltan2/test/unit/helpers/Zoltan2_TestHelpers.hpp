#include <Tpetra_config.h>
#include <Kokkos_DefaultNode.hpp>

typedef Kokkos::DefaultNode::DefaultNodeType node_t;

//
// Epetra uses int, int, double data types.  If we
// are using these data types, then we can test 
// cases of Epetra user input.
//

#define HAVE_EPETRA_DATA_TYPES

#ifdef HAVE_TPETRA_EXPLICIT_INSTANTIATION

#ifdef HAVE_TPETRA_INST_INT_LONG
typedef int lno_t;
typedef long gno_t;
#undef HAVE_EPETRA_DATA_TYPES
#else
typedef int lno_t;
typedef int gno_t;
#endif

#ifdef HAVE_TPETRA_INST_FLOAT
typedef float scalar_t;
#undef HAVE_EPETRA_DATA_TYPES
#else
typedef double scalar_t;
#endif

#else   // no tpetra explicit instantiation

typedef int lno_t;
typedef int gno_t;
typedef double scalar_t;


#endif

#include <ErrorHandlingForTests.hpp>  
#include <UserInputForTests.hpp>

