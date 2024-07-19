// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_SERIALPLATFORM_HPP
#define XPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>

namespace Xpetra {

//! \brief A implementation of the Platform class for serial platforms.
template <class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class SerialPlatform : public Teuchos::Describable {
 public:
  //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
  typedef Node NodeType;
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  explicit SerialPlatform(const Teuchos::RCP<Node> &node);

  //! Destructor
  ~SerialPlatform();

  //@}

  //! @name Class Creation and Accessor Methods
  //@{

  //! Comm Instance
  const Teuchos::RCP<const Teuchos::SerialComm<int> > getComm() const;

  //@}
 private:
  SerialPlatform(const SerialPlatform<Node> &platform);

 protected:
  //! Teuchos::Comm object instantiated for the platform.
  Teuchos::RCP<const Teuchos::SerialComm<int> > comm_;
};

template <class Node>
SerialPlatform<Node>::SerialPlatform(const Teuchos::RCP<Node> & /* node */)
  : comm_(Teuchos::rcp(new Teuchos::SerialComm<int>())) {}

template <class Node>
SerialPlatform<Node>::~SerialPlatform() {}

template <class Node>
const Teuchos::RCP<const Teuchos::SerialComm<int> >
SerialPlatform<Node>::getComm() const {
  return comm_;
}

}  // namespace Xpetra

#endif  // XPETRA_SERIALPLATFORM_HPP
