// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DEF_HPP_
#define PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DEF_HPP_

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"

#include "MueLu_FacadeClassBase_decl.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FacadeClassBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>::FacadeClassBase() {
}
}  // namespace MueLu

#endif /* PACKAGES_MUELU_SRC_INTERFACE_FACADECLASSES_MUELU_FACADECLASSBASE_DEF_HPP_ */
