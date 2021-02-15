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

#include "stk_util/environment/FileUtils.hpp"
#include "stk_util/environment/EnvData.hpp"         // for EnvData
#include "stk_util/environment/ParsedOptions.hpp"   // for ParsedOptions
#include "stk_util/environment/ProgramOptions.hpp"  // for get_parsed_options
#include <cstddef>                                  // for size_t
#include <algorithm>                                // for max


namespace {
  std::string  get_input_file_basename()
  {
    std::string filename;
    
    std::string input_file_name = "stdin";
    if (stk::get_parsed_options().count("input-deck")) {
      input_file_name = stk::EnvData::instance().m_inputFile;
    }

    std::string base = stk::util::basename(input_file_name);

    // See if name contains ".aprepro"
    size_t pos = base.find(".aprepro");
    if (pos != std::string::npos) {
      // Strip if off...
      filename = base.substr(0,pos) + base.substr(pos+8);
    } else {
      filename = base;
    }
    return filename;
  }
}

namespace stk {
  namespace util {

    std::string tailname(const std::string &filename)
    {
      size_t ind = filename.find_last_of("/", filename.size());
      if (ind != std::string::npos) {
        return filename.substr(ind+1, filename.size());
      }
      return filename; // No path, just return the filename
    }

    std::string basename(const std::string &filename)
    {
      std::string tail = tailname(filename);

      // Strip off the extension
      size_t ind = tail.find_last_of('.', tail.size());
      if (ind != std::string::npos) {
        return tail.substr(0,ind);
      }
      return tail;
    }

    void filename_substitution(std::string &filename)
    {
      // See if filename contains "%P" which is replaced by the number of processors...
      // Assumes that %P only occurs once...
      // filename is changed.
      size_t pos = filename.find("%P");
      if (pos != std::string::npos) {
	// Found the characters...  Replace with the processor count...
	int num_proc = std::max(1, EnvData::instance().m_parallelSize);
	std::string tmp(filename, 0, pos);
	tmp += std::to_string(num_proc);
	tmp += filename.substr(pos+2);
	filename = tmp;
      }

      pos = filename.find("%B");
      if (pos != std::string::npos) {
	// Found the characters...  Replace with the input file basename...
	std::string basename = get_input_file_basename();
	if (!basename.empty()) {
	  std::string tmp(filename, 0, pos);
	  tmp += basename;
	  tmp += filename.substr(pos+2);
	  filename = tmp;
	}
      }
    }
  }
}

