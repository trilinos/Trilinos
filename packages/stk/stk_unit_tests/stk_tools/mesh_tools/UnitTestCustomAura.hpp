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

#pragma once

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <string>

namespace aura_unit_tests {

class FourQuadShellsInSequenceFixture : public stk::unit_test_util::MeshFixture {
public:
  FourQuadShellsInSequenceFixture() {
    reset_mesh();
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
    const std::vector<double> coordinates {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      2.0, 0.0, 0.0,
      3.0, 0.0, 0.0,
      4.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      1.0, 1.0, 0.0,
      2.0, 1.0, 0.0,
      3.0, 1.0, 0.0,
      4.0, 1.0, 0.0
    };
    std::string mesh_description =
        "0,1,SHELL_QUAD_4,1,2,7,6,block_1\n"
        "0,2,SHELL_QUAD_4,2,3,8,7,block_1\n"
        "0,3,SHELL_QUAD_4,3,4,9,8,block_1\n"
        "0,4,SHELL_QUAD_4,4,5,10,9,block_1";
    if (this->get_parallel_size() == 2) {
      mesh_description =
          "0,1,SHELL_QUAD_4,1,2,7,6,block_1\n"
          "0,2,SHELL_QUAD_4,2,3,8,7,block_1\n"
          "0,3,SHELL_QUAD_4,3,4,9,8,block_1\n"
          "1,4,SHELL_QUAD_4,4,5,10,9,block_1";
    } else if (this->get_parallel_size() == 3) {
      mesh_description =
          "0,1,SHELL_QUAD_4,1,2,7,6,block_1\n"
          "0,2,SHELL_QUAD_4,2,3,8,7,block_1\n"
          "1,3,SHELL_QUAD_4,3,4,9,8,block_1\n"
          "2,4,SHELL_QUAD_4,4,5,10,9,block_1";
    } else if (this->get_parallel_size() == 4) {
      mesh_description =
          "0,1,SHELL_QUAD_4,1,2,7,6,block_1\n"
          "1,2,SHELL_QUAD_4,2,3,8,7,block_1\n"
          "2,3,SHELL_QUAD_4,3,4,9,8,block_1\n"
          "3,4,SHELL_QUAD_4,4,5,10,9,block_1";
    }
    stk::unit_test_util::setup_text_mesh(
          get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(mesh_description, coordinates));
  }

  void print_local_node_comm(const int rank);
};

}

