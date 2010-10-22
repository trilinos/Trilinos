#ifndef CTHULHU_DEFAULTTYPES_HPP
#define CTHULHU_DEFAULTTYPES_HPP

#include <Kokkos_DefaultNode.hpp> // Note: we should not need this header for Epetra
#include <Kokkos_DefaultKernels.hpp>

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

#endif
