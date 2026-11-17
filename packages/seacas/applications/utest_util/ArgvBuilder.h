// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <iostream>
#include <memory>
#include <vector>

namespace utest_util {

  class ArgvBuilder
  {

  public:
    ArgvBuilder() = default;

    void addArgument(const std::string &arg);

    void addArguments(const std::vector<std::string> &args);

    void clear();

    int argc() const;

    char **argv() const;

  private:
    std::vector<std::string> args_;

    mutable std::vector<char *> argv_pointers_;
  };

} // namespace utest_util
