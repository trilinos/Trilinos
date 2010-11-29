#ifndef CTHULHU_USEDEFAULTTYPES_HPP
#define CTHULHU_USEDEFAULTTYPES_HPP

#include <Kokkos_DefaultNode.hpp> // Note: we should not need this header for Epetra
#include <Kokkos_DefaultKernels.hpp>

#include "Cthulhu_ConfigDefs.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// Define default data types
typedef double ScalarType;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<ScalarType,LocalOrdinal,Node>::SparseOps LocalMatOps;

// TODO remove from this file
// Define convenient shortcut for data types
typedef ScalarType    SC;
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif
