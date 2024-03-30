// Copyright(C) 1999-2022, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "heartbeat/Iohb_Layout.h"
#include <string>

namespace Iohb {
  Layout::Layout(bool show_labels, int precision, std::string separator, int field_width)
      : layout_(), separator_(std::move(separator)), precision_(precision),
        fieldWidth_(field_width), showLabels(show_labels)
  {
  }

  void Layout::add_literal(const std::string &label) { fmt::print(layout_, "{}", label); }

  void Layout::add_legend(const std::string &label)
  {
    fmt::print(layout_, "{}{:>{}}", legendStarted ? separator_ : "", label, fieldWidth_);
    legendStarted = true;
  }
} // namespace Iohb
