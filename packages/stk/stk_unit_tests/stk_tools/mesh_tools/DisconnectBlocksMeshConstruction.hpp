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

#ifndef DISCONNECT_BLOCKS_MESH_CONSTRUCTION_HPP
#define DISCONNECT_BLOCKS_MESH_CONSTRUCTION_HPP

#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/BulkData.hpp>
#include <stk_tools/mesh_tools/DetectHingesImpl.hpp>
#include <stk_tools/mesh_tools/DisconnectTypes.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <string>
#include <algorithm>

#define EPS 0.15

stk::mesh::Part & create_part(stk::mesh::MetaData& meta, const stk::topology topology, const std::string & blockName, int64_t blockId);

//void print_node_count(stk::mesh::BulkData& bulk, const std::string str);

bool is_new_owner(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& elem, const stk::mesh::Entity& node);

void distribute_mesh(stk::mesh::BulkData& bulk, const stk::mesh::EntityIdProcVec& idProcVec);

stk::mesh::PartVector setup_mesh_1block_1quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2quad_only_on_proc_0(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2quad_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2quad_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2quad_2hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2quad_2hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3quad_1hinge_linear_stack(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3quad_1hinge_linear_stack(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4quad_bowtie_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4quad_2hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4quad_2hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4quad_4hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4quad_4hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4quad_pacman(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4quad_pacman(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4quad_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4quad_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_3quad_2tri_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_5block_3quad_2tri_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_1hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2hex_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2hex_2node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2hex_face_test(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_8hex_flower_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2tet_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_2hex_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2hex_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3hex_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3hex_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge2(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge2(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_3hex_1node_hinge_1edge_hinge3(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3hex_1node_hinge_1edge_hinge3(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_4hex_bowtie_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4hex_bowtie_1edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_two_by_two_hex_2edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_two_by_two_hex_2edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_four_hex_one_edge_one_node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_four_hex_one_edge_one_node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_four_hex_2node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_one_edge_hinge_manual(stk::mesh::BulkData& bulk);

stk::mesh::ConnectivityOrdinal destroy_element_node_relation(stk::mesh::BulkData& bulk, stk::mesh::Entity element, stk::mesh::Entity node);

void setup_mesh_with_hinge_ring(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_four_hex_2node_one_edge_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_1block_eight_tri_1node_hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4quad_bowtie_1hinge(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_3quad_1hinge(stk::mesh::BulkData& bulk);

void print_hinge_info(const stk::mesh::BulkData& bulk,
                      const stk::tools::impl::HingeNodeVector& hingeNodes,
                      const stk::tools::impl::HingeEdgeVector& hingeEdges);

bool is_debug();

// Common Decompositions
void two_elements_decomposition(stk::mesh::BulkData& bulk);

void three_elements_decomposition(stk::mesh::BulkData& bulk);

void four_elements_decomposition(stk::mesh::BulkData& bulk);

void four_elements_decomposition2(stk::mesh::BulkData& bulk);

void five_elements_decomposition(stk::mesh::BulkData& bulk);

void verify_test_run(stk::ParallelMachine pm, int locallyRanTest);

template<typename T>
inline void test_two_element_one_hinge_grouping(const stk::mesh::BulkData& bulk, const T& hinge)
{
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEMENT_RANK, 1u);
  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEMENT_RANK, 2u);
  stk::tools::impl::HingeGroupVector groupings = stk::tools::impl::get_convex_groupings(bulk, hinge);
  int ranTest = 0;

  if(!groupings.empty()) {
    ranTest = 1;
    EXPECT_EQ(groupings.size(), 2u);
    EXPECT_EQ(groupings[0].size(), 1u);
    EXPECT_EQ(groupings[1].size(), 1u);
    if(groupings[0][0] == elem1) {
      EXPECT_EQ(groupings[1][0], elem2);
    }
    if(groupings[0][0] == elem2) {
      EXPECT_EQ(groupings[1][0], elem1);
    }
  }

  verify_test_run(bulk.parallel(), ranTest);
}

stk::mesh::EntityVector get_entities_from_id_range(const stk::mesh::BulkData& bulk, stk::topology::rank_t rank, unsigned count);

stk::mesh::EntityVector get_nodes_from_id_range(const stk::mesh::BulkData& bulk, unsigned count);

stk::mesh::EntityVector get_elements_from_id_range(const stk::mesh::BulkData& bulk, unsigned count);

//void snip_and_get_remaining_hinge_count(stk::mesh::BulkData& bulk,
//                                        stk::mesh::EntityVector& elemVec,
//                                        stk::mesh::EntityVector& nodeVec,
//                                        std::pair<unsigned,unsigned>& hingeCounts);

//std::pair<unsigned,unsigned> get_locally_owned_elem_node_pair(const stk::mesh::BulkData& bulk);

std::pair<unsigned,unsigned> get_reduced_entity_counts(const stk::mesh::BulkData& bulk);

struct BlockConnection {
  BlockConnection(const std::string& b1, const std::string& b2)
    : block1(b1), block2(b2), numExpectedIntersectingNodes(0) { }

  BlockConnection(const std::string& b1, const std::string& b2, unsigned expectedNumNodes)
    : block1(b1), block2(b2), numExpectedIntersectingNodes(expectedNumNodes) { }

  std::string block1;
  std::string block2 = 0;
  unsigned numExpectedIntersectingNodes = 0;
};

typedef std::vector<BlockConnection> BlockConnectionVector;

stk::tools::BlockPairVector convert_connection_vector_to_pair_vector(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectConnVector);

stk::tools::BlockPairVector get_local_reconnect_list(const stk::mesh::BulkData& bulk, const BlockConnectionVector& disconnectList);

void create_sides_between_blocks(stk::mesh::BulkData& bulk,
                                 const std::string& block1Name,
                                 const std::string& block2Name,
                                 const std::string& sidePartName);

void create_all_boundary_sides(stk::mesh::BulkData& bulk, const std::string& sidePartName);

unsigned get_num_surface_nodes(const stk::mesh::BulkData& bulk, const std::vector<std::string>& blockPartNames);

void create_sideset(stk::mesh::BulkData& bulk,
                    const std::string& surfacePartName,
                    const std::vector<std::string>& blockPartNames);

void move_elems_from_block_to_block(stk::mesh::BulkData& bulk,
                                    const std::vector<stk::mesh::EntityId>& elemIDs,
                                    const std::string& fromBlockName,
                                    const std::string& toBlockName);

unsigned get_num_common_entities(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks, const stk::mesh::EntityRank rank);

unsigned get_num_intersecting_nodes(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks);

unsigned get_num_common_sides(const stk::mesh::BulkData & bulk, const stk::mesh::PartVector & blocks);

unsigned get_num_total_entities(const stk::mesh::BulkData & bulk, const stk::mesh::EntityRank rank);

unsigned get_num_total_nodes(const stk::mesh::BulkData & bulk);

unsigned get_num_total_sides(const stk::mesh::BulkData & bulk);

bool check_orphaned_nodes(stk::mesh::BulkData & bulk);

bool verify_attached_faces(stk::mesh::BulkData & bulk);

void output_mesh(stk::mesh::BulkData & bulk, const std::string & fileName);

void output_mesh(stk::mesh::BulkData & bulk);

int get_debug_level();

stk::mesh::PartVector setup_mesh_2block_1quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2quad_reversed(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_4quad_corner(stk::mesh::BulkData& bulk, int decompPattern = 0);

stk::mesh::PartVector setup_mesh_2block_4quad_swappedCorner(stk::mesh::BulkData& bulk, int decompPattern = 0);

stk::mesh::PartVector create_3_blocks_order1(stk::mesh::BulkData& bulk);

stk::mesh::PartVector create_3_blocks_order2(stk::mesh::BulkData& bulk);

stk::mesh::PartVector create_3_blocks_order3(stk::mesh::BulkData& bulk);

stk::mesh::PartVector create_3_blocks_order4(stk::mesh::BulkData& bulk);

stk::mesh::PartVector create_3_blocks_order5(stk::mesh::BulkData& bulk);

stk::mesh::PartVector create_3_blocks_order6(stk::mesh::BulkData& bulk);

void setup_mesh_3block_4quad_base(stk::mesh::BulkData& bulk, stk::mesh::PartVector & blocks, unsigned decompPattern);

stk::mesh::PartVector setup_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern);

void test_mesh_3block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder, unsigned decompPattern);

stk::mesh::PartVector setup_mesh_3block_4quad_keepLowerRight(stk::mesh::BulkData& bulk, unsigned decompPattern);

stk::mesh::PartVector setup_mesh_2block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern);

stk::mesh::PartVector setup_mesh_3block_4quad_checkerboard(stk::mesh::BulkData& bulk, unsigned decompPattern);

stk::mesh::PartVector setup_mesh_2block_2quad_diagonal(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4quad_bowtie(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4quad_reverse_ordinal(stk::mesh::BulkData& bulk);

void fill_mesh_description_4block_4quad_np1(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates);

void fill_mesh_description_4block_4quad_np2(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates);

void fill_mesh_description_4block_4quad_np3(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates);

void fill_mesh_description_4block_4quad_np4(stk::mesh::BulkData& bulk, unsigned blockOrder,
                                            std::string& meshDesc, std::vector<double>& coordinates);

void fill_mesh_description_6block_6quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates);

void fill_mesh_description_9block_9quad_np1(stk::mesh::BulkData& bulk, std::string& meshDesc, std::vector<double>& coordinates);

stk::mesh::PartVector setup_mesh_4block_4quad(stk::mesh::BulkData& bulk, unsigned blockOrder);

stk::mesh::PartVector setup_mesh_6block_6quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_9block_9quad(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_1hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2hex_with_internal_sides(stk::mesh::BulkData& bulk, bool loadMeshFirst = false);

stk::mesh::PartVector setup_mesh_2block_2hex_with_internal_and_external_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2hex_with_dual_internal_and_external_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2hex_with_external_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4hex_with_internal_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_2block_2cubeOfTet_with_internal_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_3block_4cubeOfTet_with_internal_sides(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4hex(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_4hex_vertical_stack(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_4block_8hex_cube(stk::mesh::BulkData& bulk);

stk::mesh::PartVector setup_mesh_8block_8hex_cube(stk::mesh::BulkData& bulk);

std::vector<std::string> get_part_names(const stk::mesh::PartVector& parts);

stk::mesh::PartVector setup_mesh_8block_8hex_with_external_sides(stk::mesh::BulkData& bulk);

#endif
