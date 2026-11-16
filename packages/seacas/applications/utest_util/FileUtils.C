// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <memory>
#include <vector>

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif
#include <unordered_map>

#include "gtest/gtest.h"

#include "FileUtils.h"

namespace {
  std::string sanitize_filename(const std::string &test_name)
  {
    std::string sanitized_name = test_name;

    // Replace invalid characters with underscores
    for (char &c : sanitized_name) {
      if (c == '/' || c == '\\' || c == ':' || c == '*' || c == '?' || c == '"' || c == '<' ||
          c == '>' || c == '|' || c == ' ') {
        c = '_'; // Replace with underscore
      }
    }

    // Optionally, limit the length of the filename
    const size_t max_length = 255; // Common max length for filenames
    if (sanitized_name.length() > max_length) {
      sanitized_name = sanitized_name.substr(0, max_length);
    }

    // Trim leading and trailing underscores (optional)
    sanitized_name.erase(0, sanitized_name.find_first_not_of('_'));
    sanitized_name.erase(sanitized_name.find_last_not_of('_') + 1);

    return sanitized_name;
  }
} // namespace

namespace utest_util {

  std::string unique_basename(const std::string &base)
  {
    std::string                basename  = base;
    const ::testing::TestInfo *test_info = ::testing::UnitTest::GetInstance()->current_test_info();
    if (test_info) {
      auto test_suffix = std::string(test_info->test_case_name()) + "." + test_info->name();
      basename += "." + sanitize_filename(test_suffix);
    }
    return basename;
  }

  std::string unique_filename(const std::string &base, const std::string &extension)
  {
    std::string filename = unique_basename(base);
    filename += "." + extension;
    return filename;
  }

} // namespace utest_util
