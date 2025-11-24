// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REPARTITIONUTILITIES_DECL_HPP
#define MUELU_REPARTITIONUTILITIES_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MPI

// Some classes are only used in the definition (_def.hpp) of this class
// but forward declarations are needed here to enable the UseShortNames mechanism.
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_CloneRepartitionInterface_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
struct RepartitionUtilities {
#undef MUELU_REPARTITIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  std::tuple<Array<GO>, GO> static ConstructGIDs(RCP<GOVector> decomposition);
};

}  // namespace MueLu

#define MUELU_REPARTITIONUTILITIES_SHORT

#endif  // ifdef HAVE_MPI
#endif  // MUELU_REPARTITIONUTILITIES_DECL_HPP
