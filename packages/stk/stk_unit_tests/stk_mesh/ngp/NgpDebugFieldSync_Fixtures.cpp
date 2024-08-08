// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "NgpDebugFieldSync_Fixtures.hpp"
#include <gtest/gtest.h>
#include "stk_util/stk_kokkos_macros.h"
#include "stk_mesh/base/FieldSyncDebugging.hpp"
#include <string>
#include <sstream>
#include <vector>

template <typename Out>
void split_lines(const std::string &s, Out result) {
  std::istringstream iss(s);
  std::string item;
  while (std::getline(iss, item)) {
    *result++ = item;
  }
}

std::vector<std::string> split_lines(const std::string &s) {
  std::vector<std::string> elems;
  split_lines(s, std::back_inserter(elems));
  return elems;
}

void extract_warning(std::string & stdoutString, int numExpectedOccurrences, const std::string & warningString)
{
  std::vector<std::string> warningLines = split_lines(stdoutString);
  std::string newStdoutString;
#if defined(STK_USE_DEVICE_MESH)
  int numFound = 0;
#endif

  for (const std::string & line : warningLines) {
    const size_t loc = line.find(warningString);
    if (loc != std::string::npos) {
#if defined(STK_USE_DEVICE_MESH)
      ++numFound;
#endif
    }
    else {
      newStdoutString += line + "\n";
    }
  }

#if defined(STK_USE_DEVICE_MESH)
  if (numFound != numExpectedOccurrences) {
    std::cout << "Warning string found " << numFound << " times when expecting " << numExpectedOccurrences << " occurrences: \""
              << warningString << "\"" << std::endl;
    ADD_FAILURE();
  }
#endif

  stdoutString = newStdoutString;
}

void check_no_warnings(const std::string & stdoutString)
{
  std::vector<std::string> warningLines = split_lines(stdoutString);

  for (const std::string & line : warningLines) {
    if (!line.empty()) {
      std::cout << "Found unexpected warning: \"" << line << "\"" << std::endl;
      ADD_FAILURE();
    }
  }
}

void check_contains_file_name(const std::string & stdoutString, const std::string & fileName)
{
#if defined(STK_USE_DEVICE_MESH) && defined(HOST_USE_LOCATION_BUILTINS)
  const size_t fileNameLoc = stdoutString.find(fileName);
  EXPECT_NE(fileNameLoc, std::string::npos) << "Could not find fileName '" << fileName
                                            << "' in stdoutString '" << stdoutString << "'" << std::endl;
#endif
}

void check_contains_a_line_number(const std::string & stdoutString)
{
#if defined(STK_USE_DEVICE_MESH) && defined(HOST_USE_LOCATION_BUILTINS)
  const size_t colonLoc = stdoutString.find(":");
  ASSERT_NE(colonLoc, std::string::npos);
  const size_t lineNumberLoc = colonLoc + 1;
  int lineNumber = std::stoi(stdoutString.data()+lineNumberLoc);
  EXPECT_GT(lineNumber, 0);
#endif
}
