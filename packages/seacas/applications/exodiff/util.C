// Copyright(C) 1999-2021, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "util.h"
#include <cstring> // for nullptr, memset
#include <fmt/color.h>
#include <fmt/ostream.h>
#include <iostream>
#include <stringx.h>
#include <unistd.h>

#include "ED_SystemInterface.h" // for SystemInterface, interFace

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include <io.h>
#define isatty _isatty
#endif

int name_length()
{
  static int max_name_length = -1;
  if (max_name_length < 0) {
    max_name_length = std::max(max_name_length, max_string_length(interFace.glob_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.node_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.elmt_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.elmt_att_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.ns_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.ss_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.eb_var_names));
    max_name_length = std::max(max_name_length, max_string_length(interFace.fb_var_names));
    max_name_length++;
  }
  return max_name_length;
}

char **get_name_array(size_t size, size_t length)
{
  char **names = nullptr;
  if (size > 0) {
    names = new char *[size];
    for (size_t i = 0; i < size; i++) {
      names[i] = new char[length + 1];
      std::memset(names[i], '\0', length + 1);
    }
  }
  return names;
}

void free_name_array(char **names, size_t size)
{
  for (size_t i = 0; i < size; i++) {
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

void Error(std::ostringstream &x)
{
  std::ostringstream out;
  fmt::print(out, "exodiff: ERROR: {}", x.str());
  ERR_OUT(out);
  exit(EXIT_FAILURE);
}

void Error(const std::string &x)
{
  std::ostringstream out;
  fmt::print(out, "exodiff: ERROR: {}", x);
  ERR_OUT(out);
  exit(EXIT_FAILURE);
}

void Warning(const std::string &x)
{
  std::ostringstream out;
  fmt::print(out, "exodiff: WARNING: {}", x);
  WARN_OUT(out);
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

void WARN_OUT(std::ostringstream &buf)
{
  if (cerr_out()) {
    fmt::print(stderr, fmt::fg(fmt::color::yellow), "{}", buf.str());
  }
  else {
    fmt::print(stderr, "{}", buf.str());
  }
}

void DIFF_OUT(std::ostringstream &buf, fmt::detail::color_type color)
{
  if (term_out()) {
    fmt::print(fmt::fg(color), "{}\n", buf.str());
  }
  else {
    fmt::print("{}\n", buf.str());
  }
}

void DIFF_OUT(const std::string &buf, fmt::detail::color_type color)
{
  if (term_out()) {
    fmt::print(fmt::fg(color), "{}\n", buf);
  }
  else {
    fmt::print("{}\n", buf);
  }
}
