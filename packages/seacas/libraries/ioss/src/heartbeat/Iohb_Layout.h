// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_Iohb_Layout_h
#define IOSS_Iohb_Layout_h

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace Iohb {
  class Layout
  {
  public:
    Layout(bool show_labels, int precision, std::string separator, int field_width);
    Layout(const Layout &) = delete;
    Layout &operator=(const Layout &) = delete;

    ~Layout();

    friend std::ostream &operator<<(std::ostream & /*o*/, Layout & /*lo*/);

    void add_literal(const std::string &label);
    void add_legend(const std::string &label);
    void add(const std::string &name, double value);
    void add(const std::string &name, int value);
    void add(const std::string &name, long value);
    void add(const std::string &name, const std::string &value);

    void add(const std::string &name, std::vector<double> &value);
    void add(const std::string &name, std::vector<int> &value);
    void add(const std::string &name, std::vector<long> &value);
    void add(const std::string &name, std::vector<std::string> &value);

  private:
    std::ostringstream layout_{};
    std::string        separator_{};

    int  precision_;
    int  count_; // Number of fields on current line...
    int  fieldWidth_;
    bool showLabels;
    bool legendStarted;
  };
} // namespace Iohb

#endif // IOSS_Iohb_Layout_h
