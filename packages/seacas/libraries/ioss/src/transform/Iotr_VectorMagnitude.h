// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iotr_export.h"

#include <Ioss_Transform.h>    // for Transform, Factory
#include <Ioss_VariableType.h> // for VariableType
#include <transform/Iotr_Factory.h>

#include <string> // for string

namespace Ioss {
  class Field;
} // namespace Ioss

namespace Iotr {

  class IOTR_EXPORT VM_Factory : public Factory
  {
  public:
    static const VM_Factory *factory();

  private:
    VM_Factory();
    Ioss::Transform *make(const std::string & /*unused*/) const override;
  };

  class IOTR_EXPORT VectorMagnitude : public Ioss::Transform
  {
    friend class VM_Factory;

  public:
    const Ioss::VariableType *output_storage(const Ioss::VariableType *in) const override;
    size_t                    output_count(size_t in) const override;

  protected:
    VectorMagnitude();

    bool internal_execute(const Ioss::Field &field, void *data) override;
  };
} // namespace Iotr
