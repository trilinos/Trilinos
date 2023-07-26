// Copyright(C) 2022, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iotr_export.h"

#include <Ioss_CodeTypes.h>
#include <functional>
#include <map>
#include <string>

namespace Iotr {
  class Factory;
  using FactoryMap = std::map<std::string, Factory *, std::less<>>;

  class IOTR_EXPORT Factory
  {
  public:
    virtual ~Factory() = default;
    static Ioss::Transform *create(const std::string &type);

    static int            describe(Ioss::NameList *names);
    static Ioss::NameList describe();

  protected:
    explicit Factory(const std::string &type);
    virtual Ioss::Transform *make(const std::string &) const = 0;
    static void              alias(const std::string &base, const std::string &syn);

  private:
    static FactoryMap &registry();
  };
} // namespace Iotr
