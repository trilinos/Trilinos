#ifndef XPETRA_USEDEFAULTTYPES_HPP
#define XPETRA_USEDEFAULTTYPES_HPP

#include <Kokkos_DefaultNode.hpp> // Note: we should not need this header for Epetra
#include <Kokkos_DefaultKernels.hpp>

#include "Xpetra_ConfigDefs.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Define default data types
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif
