// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "fmt/color.h"
#include "fmt/ostream.h"
#include "util.h"
#include <cstring> // for nullptr, memset
#include <iostream>
#include <unistd.h>

#if defined(_MSC_VER)
#include <io.h>
#define isatty _isatty
#endif

char **get_name_array(int size, int length)
{
  char **names = nullptr;
  if (size > 0) {
    names = new char *[size];
    for (int i = 0; i < size; i++) {
      names[i] = new char[length + 1];
      std::memset(names[i], '\0', length + 1);
    }
  }
  return names;
}

void free_name_array(char **names, int size)
{
  for (int i = 0; i < size; i++) {
    delete[] names[i];
  }
  delete[] names;
}

namespace {
  bool term_out()
  {
    static bool is_term = isatty(fileno(stdout));
    return is_term;
  }

  bool cerr_out()
  {
    static bool is_term = isatty(fileno(stderr));
    return is_term;
  }
} // namespace

void Error(const std::string &x)
{
  std::ostringstream out;
  fmt::print(out, "exodiff: ERROR: {}", x);
  ERR_OUT(out);
}

void ERR_OUT(std::ostringstream &buf)
{
  if (cerr_out()) {
    fmt::print(stderr, fmt::fg(fmt::color::red), "{}", buf.str());
  }
  else {
    fmt::print(stderr, "{}", buf.str());
  }
}

void DIFF_OUT(std::ostringstream &buf, fmt::internal::color_type color)
{
  if (term_out()) {
    fmt::print(fmt::fg(color), "{}\n", buf.str());
  }
  else {
    fmt::print("{}\n", buf.str());
  }
}

void DIFF_OUT(const std::string &buf, fmt::internal::color_type color)
{
  if (term_out()) {
    fmt::print(fmt::fg(color), "{}\n", buf);
  }
  else {
    fmt::print("{}\n", buf);
  }
}
