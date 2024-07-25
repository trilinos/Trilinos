// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef FAKEKOKKOS_DEFAULT_NODE_HPP_
#define FAKEKOKKOS_DEFAULT_NODE_HPP_

#include <Teuchos_RCP.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

// forward declaration for (fake) KokkosSerialWrapperNode
// This is the node definition used if Epetra is enabled only
/*namespace Kokkos {
namespace Compat {
  class KokkosSerialWrapperNode;
}
}*/

// This KokkosClassic namespace is used for getting the DefaultNode in some classes
namespace KokkosClassic {

namespace Details {
}  // namespace Details

class DefaultNode {
 public:
#ifdef EPETRA_HAVE_OMP
  typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode DefaultNodeType;
#else
  typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode DefaultNodeType;
#endif
};

}  // namespace KokkosClassic

#endif
