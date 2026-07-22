// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MPIPLATFORM_HPP
#define XPETRA_MPIPLATFORM_HPP

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"

namespace Xpetra {

//! \brief A implementation of the Platform class for MPI-based platforms.
/*!
  This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal.
  The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal
  type, if omitted, defaults to the \c LocalOrdinal type.
*/
template <class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MpiPlatform : public Teuchos::Describable {
 public:
  //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
  typedef Node NodeType;
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  explicit MpiPlatform(Teuchos::RCP<Node> node);

  //! Constructor
  MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

  //! Destructor
  ~MpiPlatform();

  //@}

  //! @name Class Creation and Accessor Methods
  //@{

  //! Comm Instance
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //@}

 private:
  Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
  MpiPlatform(const MpiPlatform<Node> &platform);
};

template <class Node>
MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> /* node */, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
  : comm_(Teuchos::createMpiComm<int>(rawMpiComm)) {}

template <class Node>
MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> /* node */)
  : comm_(Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD))) {}  // CHECK: ALLOW MPI_COMM_WORLD

template <class Node>
MpiPlatform<Node>::~MpiPlatform() {}

template <class Node>
MpiPlatform<Node>::MpiPlatform(const MpiPlatform<Node> &platform) {
  comm_ = platform.comm_;
}

template <class Node>
Teuchos::RCP<const Teuchos::Comm<int> >
MpiPlatform<Node>::getComm() const {
  return comm_;
}

}  // namespace Xpetra

#endif  // XPETRA_MPIPLATFORM_HPP
