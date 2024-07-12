// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEKOSMOOTHER_DEF_HPP_
#define MUELU_TEKOSMOOTHER_DEF_HPP_

#ifdef HAVE_MUELU_TEKO

#include "MueLu_TekoSmoother_decl.hpp"

namespace MueLu {
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNodeSmootherComplexity() const {
  // FIXME: This is a placeholder
  return Teuchos::OrdinalTraits<size_t>::invalid();
}

}  // namespace MueLu

#endif  // HAVE_MUELU_TEKO

#endif /* MUELU_TEKOSMOOTHER_DEF_HPP_ */
