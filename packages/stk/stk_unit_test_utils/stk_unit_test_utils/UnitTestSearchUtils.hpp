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

#ifndef STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTSEARCHUTILS_HPP_
#define STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTSEARCHUTILS_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_mesh/base/Types.hpp"                    // for EntityId, PartV...
#include "stk_search_util/spmd/ElementSendMesh.hpp"
#include "stk_search_util/spmd/ElementRecvMesh.hpp"
#include "stk_search_util/spmd/NodeSendMesh.hpp"
#include "stk_search_util/spmd/NodeRecvMesh.hpp"
#include "stk_search/Box.hpp"                         // for Box
#include "stk_search/FilterCoarseSearch.hpp"          // for ObjectOutsideDo...
#include "stk_search/IdentProc.hpp"                   // for IdentProc
#include "stk_search/Point.hpp"                       // for Point
#include "stk_search/Sphere.hpp"                      // for Sphere
#include "stk_topology/topology.hpp"                  // for topology
#include "stk_unit_test_utils/BuildMesh.hpp"               // for build_mesh
#include "stk_unit_test_utils/MockMasterElementProvider.hpp"
#include "stk_unit_test_utils/MockSearchMesh.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"                // for get_full_t...
#include "stk_util/parallel/Parallel.hpp"             // for ParallelMachine

#include <memory>                                     // for shared_ptr
#include <set>                                        // for set
#include <string>                                     // for string, basic_s...
#include <utility>                                    // for pair
#include <vector>                                     // for vector

namespace stk::mesh { class BulkData; }
namespace stk::mesh { class MetaData; }
namespace stk::mesh { class FieldBase; }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace unit_test_util {

stk::mesh::PartVector get_parts_list(const stk::mesh::MetaData& meta, const std::vector<std::string>& blocks);
std::vector<std::string> get_all_block_names(const stk::mesh::MetaData& meta);

std::shared_ptr<stk::mesh::BulkData> build_slanted_single_hex_bulk(double slant = 0.5);

std::shared_ptr<stk::mesh::BulkData> build_single_point_bulk(double x, double y, double z);

std::shared_ptr<stk::search::spmd::ElementSendMesh>
construct_hex_send_mesh(stk::mesh::BulkData& bulk, double parametricTolerance,
                        const std::vector<std::string>& blocks = { "block_1" },
                        const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement());

std::shared_ptr<stk::search::spmd::NodeRecvMesh>
construct_node_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                        const std::vector<std::string>& blocks = { "block_1" },
                        const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement());

std::shared_ptr<stk::search::spmd::ElementRecvMesh>
construct_element_centroid_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                                     const std::vector<std::string>& blocks = { "block_1" },
                                     const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement());

std::shared_ptr<stk::search::spmd::ElementRecvMesh>
construct_hex_gauss_point_recv_mesh(stk::mesh::BulkData& bulk, double parametricTolerance, double geometricTolerance,
                                    unsigned integrationOrder = 0,
                                    const std::vector<std::string>& blocks = { "block_1" },
                                    const stk::mesh::Selector& activeSelector = stk::mesh::Selector().complement());

}
}

#endif /* STK_STK_UNIT_TEST_UTILS_STK_UNIT_TEST_UTILS_UNITTESTSEARCHUTILS_HPP_ */
