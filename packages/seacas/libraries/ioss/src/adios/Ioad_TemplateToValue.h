// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Ioad_TemplateToValue_h
#define IOSS_Ioad_TemplateToValue_h

namespace Ioad {

  template <typename T> constexpr Ioss::Field::BasicType template_to_basic_type() noexcept;

  // Only implement for a few specialized classes on purpose. No generic implementation.
  template <typename T> constexpr char const *get_entity_type() noexcept;

  template <class T> inline std::string GetType() noexcept;

} // namespace Ioad

#include "adios/Ioad_TemplateToValue.hpp"

#endif
