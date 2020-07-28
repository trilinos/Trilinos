/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef TOKENIZE_H
#define TOKENIZE_H

#include <string>
#include <vector>

/**
 * Take the 'str' argument and split it using the list of characters
 * in separators as separators. Return tokens as a vector of strings.
 */
namespace SLIB {
  /**
   * If `allow_empty_token` is false, then multiple sequential delimiters will not produce an empty
   * token,
   * If it is true, then there is a token between each and every delimiter even if empty.
   * If | is delimiter, then when false: a|||b is two tokens `a` and `b`.
   * When true, a|||b is 4 tokens "a" "" "" "b"
   */
  std::vector<std::string> tokenize(const std::string &str, const std::string &separators,
                                    bool allow_empty_token = false);
} // namespace SLIB

#endif
