/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <string>
#include <vector>

namespace SEAMS {
  /**
   * Take the 'str' argument and split it using the list of characters
   * in separators as separators. Use tokens to return the result.
   */
  std::vector<std::string> tokenize(const std::string &str, const std::string &separators);
} // namespace SEAMS
