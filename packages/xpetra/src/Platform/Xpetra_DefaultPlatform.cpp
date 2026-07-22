// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

Teuchos::RCP<Xpetra::DefaultPlatform::DefaultPlatformType> Xpetra::DefaultPlatform::platform_ = Teuchos::null;

namespace Xpetra {

DefaultPlatform::DefaultPlatformType &DefaultPlatform::getDefaultPlatform() {
  XPETRA_MONITOR("DefaultPlatform::getDefaultPlatform");

  if (!platform_.get()) {
    using node_type = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;
    Teuchos::RCP<node_type> node;  // null
#ifdef HAVE_MPI
    platform_ = Teuchos::rcp(new MpiPlatform<node_type>(node));
#else
    platform_ = Teuchos::rcp(new SerialPlatform<node_type>(node));
#endif
  }
  return *platform_;
}

}  // namespace Xpetra
