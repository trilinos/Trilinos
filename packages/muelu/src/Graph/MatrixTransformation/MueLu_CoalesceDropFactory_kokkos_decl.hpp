// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP
#define MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP

#include "MueLu_CoalesceDropFactory.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class [[deprecated]] CoalesceDropFactory_kokkos : public CoalesceDropFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> {};

}  // namespace MueLu

#define MUELU_COALESCEDROPFACTORY_KOKKOS_SHORT
#endif  // MUELU_COALESCEDROPFACTORY_KOKKOS_DECL_HPP
