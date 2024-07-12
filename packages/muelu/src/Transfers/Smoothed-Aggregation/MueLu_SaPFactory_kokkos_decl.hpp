// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SAPFACTORY_KOKKOS_DECL_HPP
#define MUELU_SAPFACTORY_KOKKOS_DECL_HPP

#include "MueLu_SaPFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class [[deprecated]] SaPFactory_kokkos : public SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> {};

}  // namespace MueLu

#define MUELU_SAPFACTORY_KOKKOS_SHORT
#endif  // MUELU_SAPFACTORY_KOKKOS_DECL_HPP
