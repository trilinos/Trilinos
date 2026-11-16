// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ArgvBuilder.h"

namespace utest_util {

  void ArgvBuilder::addArgument(const std::string &arg) { args_.push_back(arg); }

  void ArgvBuilder::addArguments(const std::vector<std::string> &args)
  {
    for (const auto &arg : args) {
      args_.push_back(arg);
    }
  }

  void ArgvBuilder::clear()
  {
    args_.clear();
    argv_pointers_.clear();
  }

  int ArgvBuilder::argc() const { return static_cast<int>(args_.size()); }

  char **ArgvBuilder::argv() const
  {
    argv_pointers_.clear();
    argv_pointers_.reserve(args_.size() + 1);

    for (const auto &arg : args_) {

      argv_pointers_.push_back(const_cast<char *>(arg.c_str()));
    }

    argv_pointers_.push_back(nullptr);

    return argv_pointers_.data();
  }

} // namespace utest_util