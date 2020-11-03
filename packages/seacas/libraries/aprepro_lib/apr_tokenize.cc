// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "apr_tokenize.h"
#include <algorithm>

std::vector<std::string> SEAMS::tokenize(const std::string &str, const std::string &separators)
{
  std::vector<std::string> tokens;
  auto                     first = std::begin(str);
  while (first != std::end(str)) {
    const auto second =
        std::find_first_of(first, std::end(str), std::begin(separators), std::end(separators));
    if (first != second) {
      tokens.emplace_back(first, second);
    }
    if (second == std::end(str)) {
      break;
    }
    first = std::next(second);
  }
  return tokens;
}

#if 0
#include <iostream>
typedef std::vector<std::string> TokenList;

int main()
{
  char s[128];
  while(!std::cin.eof()) {
    std::cout << "Enter a string: ";
    std::cin.getline(s,128);
    std::string input_line(s);
    if (input_line != "quit") {
      std::vector<std::string> tokens = tokenize(input_line, ": \t\r\v\n");
      std::cout << "There were " << tokens.size() << " tokens in the line\n";
      TokenList::const_iterator I = tokens.begin();
      while (I != tokens.end()) {
        std::cout << "'" << *I++ << "'\t";
      }
      std::cout << '\n';
    } else {
      exit(0);
    }
  }
}
#endif
