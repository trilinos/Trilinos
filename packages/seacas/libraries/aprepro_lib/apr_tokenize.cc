// Copyright (c) 2014, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "apr_tokenize.h"
#include <algorithm>

void SEAMS::tokenize(const std::string& str, const std::string& separators,
              std::vector<std::string>& tokens)
{
  std::string curr_token = "";
  for (size_t i = 0; i < str.length(); ++i) {
    char curr_char = str[i];

    // determine if current character is a separator
    bool is_separator = std::find(separators.begin(), separators.end(), curr_char) != separators.end();

    if (is_separator && curr_token != "") {
      // we just completed a token
      tokens.push_back(curr_token);
      curr_token.clear();
    }
    else if (!is_separator) {
      curr_token += curr_char;
    }
  }
  if (curr_token != "") {
    tokens.push_back(curr_token);
  }
}

#if 0
#include <iostream>
using std::cout;
using std::cin;

typedef std::vector<std::string> TokenList;

int main()
{
  char s[128];
  while(!cin.eof()) {
    cout << "Enter a string: ";
    cin.getline(s,128);
    std::string input_line(s);
    if (input_line != "quit") {
      std::vector<std::string> tokens;
      tokenize(input_line, ": \t\r\v\n", tokens);
      cout << "There were " << tokens.size() << " tokens in the line\n";
      TokenList::const_iterator I = tokens.begin();
      while (I != tokens.end()) {
        cout << "'" << *I++ << "'\t";
      }
      cout << '\n';
    } else {
      exit(0);
    }
  }
}
#endif
