// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MEMORY_PROFILER_HPP
#define MUELU_MEMORY_PROFILER_HPP

#include <string>

std::string GetMemoryUsage();
void PrintMemoryUsage(const std::string& description = "Memory Usage:", const std::string& filename = "");

void MemoryUsageStart(const std::string& autoLogPrefix = "memorylog-");
void MemoryUsageStop();

#endif  // MUELU_MEMORY_PROFILER_HPP
