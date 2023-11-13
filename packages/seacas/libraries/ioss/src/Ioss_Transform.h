// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ioss_export.h"

#include <Ioss_CodeTypes.h>

#include <string>
#include <vector>

namespace Ioss {
  class Field;
  class VariableType;

  class IOSS_EXPORT Transform
  {
  public:
    virtual ~Transform();
    virtual const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const = 0;
    virtual size_t                    output_count(size_t in) const                      = 0;

    bool execute(const Ioss::Field &field, void *data);

    virtual void set_property(const std::string &name, int value);
    virtual void set_property(const std::string &name, double value);
    virtual void set_properties(const std::string &name, const std::vector<int> &values);
    virtual void set_properties(const std::string &name, const std::vector<double> &values);

  protected:
    Transform();

    virtual bool internal_execute(const Ioss::Field &field, void *data) = 0;
  };
} // namespace Ioss
