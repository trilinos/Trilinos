// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <sstream>
#include <string>
#include <vector>

#include "iohb_export.h"

namespace Iohb {
  class IOHB_EXPORT Layout
  {
  public:
    Layout(bool show_labels, int precision, std::string separator, int field_width);
    Layout(const Layout &)            = delete;
    Layout &operator=(const Layout &) = delete;

    const std::string layout() const { return layout_.str(); }

    void add_literal(const std::string &label);
    void add_legend(const std::string &label);

    template <typename T> void add(const std::string &name, const T &value);
    template <typename T> void add(const std::string &name, const std::vector<T> &value);

  private:
    void               output_common(const std::string &name);
    std::ostringstream layout_{};
    std::string        separator_{", "};

    int  precision_{5};
    int  count_{0}; // Number of fields on current line...
    int  fieldWidth_{0};
    bool showLabels{true};
    bool legendStarted{false};
  };

  inline void Layout::output_common(const std::string &name)
  {
    if (count_++ > 0 && !separator_.empty()) {
      fmt::print(layout_, "{}", separator_);
    }

    if (showLabels && !name.empty()) {
      fmt::print(layout_, "{}=", name);
    }
  }

  template <typename T> inline void Layout::add(const std::string &name, const T &value)
  {
    output_common(name);
    if (!showLabels && fieldWidth_ > 0) {
      fmt::print(layout_, "{0:{1}}", value, fieldWidth_);
    }
    else {
      fmt::print(layout_, "{}", value);
    }
  }

  template <> inline void Layout::add(const std::string &name, const double &value)
  {
    output_common(name);
    if (precision_ == -1) {
      // Use lib::fmt full precision output -- as many digits as needed to fully represent the
      // double
      fmt::print(layout_, "{}", value);
    }
    else if (!showLabels && fieldWidth_ > 0) {
      fmt::print(layout_, "{0:{1}.{2}e}", value, fieldWidth_, precision_);
    }
    else {
      fmt::print(layout_, "{0:.{1}e}", value, precision_);
    }
  }

  template <typename T>
  inline void Layout::add(const std::string &name, const std::vector<T> &value)
  {
    if (value.size() == 1) {
      add(name, value[0]);
    }
    else {
      output_common(name);
      if (!showLabels && fieldWidth_ > 0) {
        fmt::print(layout_, "{0:{1}}", fmt::join(value, separator_), fieldWidth_);
      }
      else {
        fmt::print(layout_, "{}", fmt::join(value, separator_));
      }
    }
  }

  template <> inline void Layout::add(const std::string &name, const std::vector<double> &value)
  {
    if (value.size() == 1) {
      add(name, value[0]);
    }
    else {
      output_common(name);
      if (precision_ == -1) {
        // Use lib::fmt full precision output -- as many digits as needed to fully represent the
        // double
        fmt::print(layout_, "{}", fmt::join(value, separator_));
      }
      else if (!showLabels && fieldWidth_ > 0) {
        fmt::print(layout_, "{0:{2}.{1}e}", fmt::join(value, separator_), precision_, fieldWidth_);
      }
      else {
        fmt::print(layout_, "{0:.{1}e}", fmt::join(value, separator_), precision_);
      }
    }
  }

} // namespace Iohb
