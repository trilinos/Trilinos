// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SMOOTHERPROTOTYPE_DEF_HPP
#define MUELU_SMOOTHERPROTOTYPE_DEF_HPP

#include "MueLu_SmootherPrototype_decl.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SmootherPrototype()
  : isSetup_(false) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SmootherPrototype() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsSetup() const { return isSetup_; }

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node>::IsSetup(bool const &ToF) { isSetup_ = ToF; }

}  // namespace MueLu

// TODO: private copy constructor
// TODO: update comments
#endif  // MUELU_SMOOTHERPROTOTYPE_DEF_HPP
