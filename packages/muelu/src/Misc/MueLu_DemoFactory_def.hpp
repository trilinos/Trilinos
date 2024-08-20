// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DEMOFACTORY_DEF_HPP
#define MUELU_DEMOFACTORY_DEF_HPP

#include "MueLu_DemoFactory_decl.hpp"

// #include <Xpetra_Matrix.hpp>

#include "MueLu_Level.hpp"
// #include "MueLu_Monitor.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DemoFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~DemoFactory() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* currentLevel */) const {
  // TODO: declare input for factory
  // Input(currentLevel, varName_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void DemoFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* currentLevel */) const {
  // TODO: implement factory
}

}  // namespace MueLu

#define MUELU_DEMOFACTORY_SHORT
#endif  // MUELU_DEMOFACTORY_DEF_HPP
