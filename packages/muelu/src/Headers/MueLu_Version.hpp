// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_VERSION_HPP
#define MUELU_VERSION_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

inline std::string const Version() {
  return ("MueLu development");
}

}  // namespace MueLu

#endif  // ifndef MUELU_VERSION_HPP
