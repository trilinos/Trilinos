// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef FAKEKOKKOS_SERIALNODE_HPP_
#define FAKEKOKKOS_SERIALNODE_HPP_

#include <Teuchos_ParameterList.hpp>

// This is the node definition used if Epetra is enabled only
// Epetra can be compiled using the SerialNode OR the OpenMP Node
// If there is no Kokkos provide dummy classes for all Kokkos node
// types that Epetra might be compiled for.
namespace Kokkos {
namespace Compat {
class KokkosSerialWrapperNode {
 public:
  KokkosSerialWrapperNode(Teuchos::ParameterList &pl) {}
  KokkosSerialWrapperNode() {}
};

class KokkosOpenMPWrapperNode {
 public:
  KokkosOpenMPWrapperNode(Teuchos::ParameterList &pl) {}
  KokkosOpenMPWrapperNode() {}
};
}  // namespace Compat
}  // namespace Kokkos

#endif
