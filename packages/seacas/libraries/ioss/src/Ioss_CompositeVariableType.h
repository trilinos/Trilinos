// Copyright(C) 1999-2020, 2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "ioss_export.h"

#include "Ioss_CodeTypes.h"
#include "Ioss_VariableType.h" // for VariableType
#include <string>              // for string

namespace Ioss {
  class IOSS_EXPORT CompositeVariableType : public VariableType
  {
  public:
    IOSS_NODISCARD static std::string   composite_name(const std::string &base, int copies);
    IOSS_NODISCARD static VariableType *composite_variable_type(const VariableType *inst,
                                                                int                 copies);

    IOSS_NODISCARD std::string label(int which, char suffix_sep = '_') const override;
    CompositeVariableType(const std::string &my_name, int number_components, bool delete_me);
    CompositeVariableType(const VariableType *base_type, int copies, bool delete_me);
    CompositeVariableType(const CompositeVariableType &) = delete;

    IOSS_NODISCARD const VariableType *GetBaseType() const;
    IOSS_NODISCARD int                 GetNumCopies() const;

  private:
    const VariableType *baseType;
    int                 copies_;
  };
} // namespace Ioss
