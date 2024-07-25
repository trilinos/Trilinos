// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOSCOMPAT_DEFAULTNODE_HPP
#define TPETRA_KOKKOSCOMPAT_DEFAULTNODE_HPP

#include "TpetraCore_config.h"
#include "Tpetra_KokkosClassic_DefaultNode_config.h"
#include "Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Teuchos_RCP.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Dear users: This is just a forward declaration.
  // Please skip over it.
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace Tpetra {
namespace KokkosClassic {


  /// \brief Specify Tpetra's default Node type.
  ///
  /// Tpetra::Map uses this class to get Tpetra's default Node type.
  /// <i>This is an implementation detail of Tpetra</i>.  If you want
  /// to know the default Node type, just ask Tpetra::Map, like this:
  /// \code
  /// typedef Tpetra::Map<>::node_type default_node_type;
  /// \endcode
  class DefaultNode {
  public:
#if defined(HAVE_TPETRA_DEFAULTNODE_SYCLWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosSYCLWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_HIPWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosHIPWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosCudaWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosOpenMPWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosThreadsWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
    typedef Tpetra::KokkosCompat::KokkosSerialWrapperNode DefaultNodeType;
#else
#    error "No default Kokkos Node type specified.  Please set the CMake option Tpetra_DefaultNode to a valid Node type."
#endif

  };

} // namespace KokkosClassic
} // namespace Tpetra

#endif // TPETRA_KOKKOSCOMPAT_DEFAULTNODE_HPP
