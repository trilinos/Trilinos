// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MEMORY_HPP
#define MUELU_MEMORY_HPP

#include <string>
#include "Teuchos_Time.hpp"
#include "Teuchos_Comm.hpp"

namespace Teuchos {
class Time;
}
namespace Teuchos {
template <typename Ordinal>
class Comm;
}

namespace MueLu {

namespace MemUtils {

std::string PrintMemoryUsage();
std::string PrintMemoryInfo();
void ReportTimeAndMemory(Teuchos::Time const &timer, Teuchos::Comm<int> const &Comm);

}  // namespace MemUtils

}  // namespace MueLu

#endif  // ifndef MUELU_MEMORY_HPP
