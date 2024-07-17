// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_KEEPTYPE_HPP
#define MUELU_KEEPTYPE_HPP

namespace MueLu {

//! Keep status of a variable of Level.
// Several keep status can be set at the same time
enum KeepEnum {
  UserData = 0x1,  //!< User data are always kept. This flag is set automatically when Level::Set("data", data) is used. The keep status of the variable is not propagated to coarser level (if you use Level::Build()).
  Keep     = 0x2,  //!< Always keep data, even accross run. This flag is set by Level::Keep(). This flag is propagated to coarser level by Level::Build().
  Final    = 0x4,  //!< Keep data only for this run. Used to keep data useful for Hierarchy::Iterate(). Data will be deleted if setup phase is re-run. This flag is set by default for A, P, R, PreSmoother and PostSmoother of NoFactory by Hierarchy::Setup(). Not propagated by Level::Build().

  NextRun = UserData | Keep,  //!< Both UserData and Keep flags force data to be kept and reused for the next run. Do not use MueLu::NextRun in AddKeepFlag. Use it only for testing keep == UserData || keep == Keep.
  All     = UserData | Keep | Final
};

//!
typedef short KeepType;  // TODO: name it KeepFlag?

}  // namespace MueLu

#endif  // MUELU_KEEPTYPE_HPP
