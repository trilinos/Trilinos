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


#ifndef STK_GET_MESH_SPEC_H
#define STK_GET_MESH_SPEC_H

#include <stk_unit_test_utils/getOption.h>
#include <stk_util/stk_config.h>
#include <string>
#include <vector>

namespace stk
{
namespace unit_test_util
{

inline
std::string get_mesh_spec(unsigned dim)
{
  std::string dimString = std::to_string(dim);
  std::string meshSpec("generated:");
  meshSpec += dimString + "x" + dimString + "x" + dimString;
  return meshSpec;
}

inline
std::string get_mesh_spec(const std::string &optionName)
{
  return get_mesh_spec(stk::unit_test_util::get_command_line_option<unsigned>(optionName, 20));
}

inline std::vector<double> get_many_block_coordinates(unsigned numBlocks)
{
  std::vector<double> planeCoords = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };

  std::vector<double> coordinates;

  coordinates.insert(coordinates.end(), planeCoords.begin(), planeCoords.end());

  for(unsigned i = 1; i <= numBlocks; ++i) {
    for(unsigned point = 0; point < 4; ++point) {
      planeCoords[3 * point + 2] += 1;
    }

    coordinates.insert(coordinates.end(), planeCoords.begin(), planeCoords.end());
  }

  return coordinates;
}

inline void get_block_proc_distribution(unsigned numBlocks, unsigned numProc, std::vector<unsigned>& procs)
{
  procs.resize(numBlocks);
  unsigned unload = numBlocks % numProc;
  unsigned initVal = numBlocks / numProc;
  std::vector<unsigned> numBlocksPerProc(numProc);

  std::fill(numBlocksPerProc.begin(), numBlocksPerProc.end(), initVal);
  
  for(unsigned i = 0; i < unload; i++) {
    numBlocksPerProc[i]++;
  }
  for(unsigned i = 0; i < numBlocks; i++) {
    for(unsigned j = 0; j < numProc; j++) {
      if(numBlocksPerProc[j] > 0) {
        numBlocksPerProc[j]--;
        procs[i] = j;
        break;
      }
    }
  }
}

inline std::string get_many_block_mesh_desc(unsigned numBlocks, unsigned numProc = 1)
{
  std::ostringstream oss;
  std::vector<unsigned> procs(numBlocks);
  get_block_proc_distribution(numBlocks, numProc, procs);

  unsigned proc = 0;
  for(unsigned i = 0; i < numBlocks; ++i) {
    proc = procs[i];
    unsigned elemId = i + 1;
    unsigned firstNodeId = i * 4 + 1;
    oss << proc << "," << elemId << ",HEX_8,";
    for(unsigned node = firstNodeId; node < firstNodeId + 8; ++node) {
      oss << node << ",";
    }
    unsigned blockId = i + 1;
    oss << "block_" << blockId;

    if(i < numBlocks - 1) {
      oss << "\n";
    }
  }

  return oss.str();
}

namespace simple_fields {

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::string get_mesh_spec(unsigned dim);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::string get_mesh_spec(const std::string &optionName);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::vector<double> get_many_block_coordinates(unsigned numBlocks);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void get_block_proc_distribution(unsigned numBlocks, unsigned numProc, std::vector<unsigned>& procs);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::string get_many_block_mesh_desc(unsigned numBlocks, unsigned numProc = 1);

} // namespace simple_fields

}
}

#endif
