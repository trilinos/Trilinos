// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include <gtest/gtest.h>
#include <vector>                         // for vector
#include <stk_mesh/base/BulkData.hpp>     // for BulkData
#include <stk_util/parallel/Parallel.hpp> // for ParallelMachine
#include "stk_mesh/base/Selector.hpp"     // for Selector
#include <stk_mesh/base/CreateEdges.hpp>  // for create_edges
#include <stk_mesh/base/CreateFaces.hpp>  // for create_edges
#include "stk_mesh/base/Entity.hpp"       // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include "stk_mesh/base/MetaData.hpp"     // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"        // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"      // for topology, etc
#include "stk_mesh/base/GetEntities.hpp"  // for count_entities
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Comm.hpp>
#include <array>

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;

TEST(FEMHelper, get_ordinal_and_permutation)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        unsigned gold_side_node_ids[4] = {5,6,8,7};
        unsigned gold_num_nodes = 4;

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        std::string name = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        unsigned elem_id = 0;
        unsigned gold_local_side_id = 0;
        unsigned perm_value = 0;

        if (mesh.parallel_rank()==0)
        {
            gold_local_side_id=5;
            elem_id = 1;
            perm_value = 0;
        }
        else
        {
            gold_local_side_id=4;
            elem_id = 2;
            perm_value = 4;
        }

        stk::mesh::Permutation gold_permutation = static_cast<stk::mesh::Permutation>(perm_value);

        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);
        EXPECT_TRUE(mesh.bucket(elem).owned());

        stk::mesh::EntityVector side_nodes(gold_num_nodes);
        for(unsigned i = 0; i < gold_num_nodes; ++i)
        {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i]);
            side_nodes[i] = node;
        }

        stk::mesh::EntityRank to_rank = stk::topology::FACE_RANK;
        std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(mesh, elem, to_rank, side_nodes);

        ASSERT_TRUE(ordinalAndPermutation.first != stk::mesh::INVALID_CONNECTIVITY_ORDINAL);
        ASSERT_TRUE(ordinalAndPermutation.second != stk::mesh::INVALID_PERMUTATION);

        EXPECT_EQ(gold_local_side_id, ordinalAndPermutation.first);
        EXPECT_EQ(gold_permutation, ordinalAndPermutation.second);
    }
}

TEST(FEMHelper, check_permutation_consistency_using_FEMHelper_parallel)
{
    if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
    {
        stk::mesh::EntityId global_side_id = 1;
        unsigned gold_side_node_ids[4] = {5,6,8,7};
        unsigned gold_num_nodes = 4;

        stk::io::StkMeshIoBroker stkMeshIoBroker(MPI_COMM_WORLD);
        std::string name = "generated:1x1x2";
        stkMeshIoBroker.add_mesh_database(name, stk::io::READ_MESH);
        stkMeshIoBroker.create_input_mesh();
        stkMeshIoBroker.populate_bulk_data();

        stk::mesh::BulkData &mesh = stkMeshIoBroker.bulk_data();

        unsigned elem_id = 0;

        if (mesh.parallel_rank()==0)
        {
            elem_id = 1;
        }
        else
        {
            elem_id = 2;
        }

        stk::mesh::Entity elem = mesh.get_entity(stk::topology::ELEM_RANK, elem_id);
        EXPECT_TRUE(mesh.bucket(elem).owned());

        stk::mesh::EntityVector side_nodes(gold_num_nodes);
        for(unsigned i = 0; i < gold_num_nodes; ++i) {
            stk::mesh::Entity node = mesh.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i]);
            side_nodes[i] = node;
        }

        stk::mesh::Part &part = mesh.mesh_meta_data().get_topology_root_part(stk::topology::QUAD_4);
        mesh.modification_begin();
        stk::mesh::Entity side = declare_element_to_sub_topology_with_nodes(mesh, elem, side_nodes, global_side_id, stk::topology::FACE_RANK, part);
        EXPECT_NO_THROW(mesh.modification_end());

        std::vector<size_t> mesh_counts;
        stk::mesh::comm_mesh_counts(mesh, mesh_counts);
        size_t numFacesGlobal = 1u;
        EXPECT_EQ(numFacesGlobal, mesh_counts[stk::topology::FACE_RANK]);
        EXPECT_TRUE(mesh.is_valid(side));
    }
}
void build_element_from_topology_verify_ordinals_and_permutations(stk::mesh::BulkData &bulk,
		const stk::topology topo,
		stk::mesh::EntityId * elem_node_ids,
		stk::mesh::EntityId * side_ids,
		stk::mesh::EntityId * edge_ids,
		const std::vector < std::vector < unsigned > > &gold_side_node_ids,
		const unsigned * gold_side_permutations,
		const std::vector < std::vector < unsigned > > &gold_edge_node_ids,
		const unsigned * gold_edge_permutations)
{
	stk::mesh::EntityId element_id[1] = {1};
	stk::mesh::MetaData &meta = bulk.mesh_meta_data();
	stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part", topo);

	meta.commit();
	bulk.modification_begin();

	stk::mesh::Entity elem = stk::mesh::declare_element(bulk, elem_part, element_id[0], elem_node_ids);

	stk::mesh::EntityVector side_nodes;
	uint num_sides = topo.num_sides();
	stk::topology::rank_t sub_topo_rank = topo.side_rank();

	for(uint i = 0; i < num_sides; ++i)
	{
		stk::topology sub_topo = topo.side_topology(i);
		side_nodes.clear();

		stk::mesh::Entity side = bulk.declare_entity(sub_topo_rank, side_ids[i], meta.get_topology_root_part(sub_topo));

		for (uint j = 0; j < sub_topo.num_nodes(); ++j)
		{
			stk::mesh::Entity side_node = bulk.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i][j]);
			side_nodes.push_back(side_node);
			bulk.declare_relation(side, side_node, j);
		}

		std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(bulk, elem, sub_topo_rank, side_nodes);

		EXPECT_EQ(ordinalAndPermutation.second, gold_side_permutations[i]);
		EXPECT_EQ(ordinalAndPermutation.first, i);
	}

	if (edge_ids == NULL) {
		bulk.modification_end();
		return;
	}

	stk::mesh::EntityVector edge_nodes;
	uint num_edges = topo.num_edges();

	for(uint i = 0; i < num_edges; ++i)
	{
		edge_nodes.clear();

		stk::mesh::Entity edge = bulk.declare_entity(stk::topology::EDGE_RANK, edge_ids[i], meta.get_topology_root_part(topo.edge_topology()));

		for (uint j = 0; j < topo.edge_topology().num_nodes(); ++j)
		{
			stk::mesh::Entity edge_node = bulk.get_entity(stk::topology::NODE_RANK, gold_edge_node_ids[i][j]);
			edge_nodes.push_back(edge_node);
			bulk.declare_relation(edge, edge_node, j);
		}

		std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(bulk, elem, stk::topology::EDGE_RANK, edge_nodes);

		EXPECT_EQ(ordinalAndPermutation.second, gold_edge_permutations[i]);
		EXPECT_EQ(ordinalAndPermutation.first, i);
	}

	bulk.modification_end();
}

void verify_unbuildable_element(stk::mesh::BulkData &bulk,
		const stk::topology topo,
		stk::mesh::EntityId * elem_node_ids,
		stk::mesh::EntityId * side_ids,
		const std::vector < std::vector < unsigned > > &gold_side_node_ids,
		bool *sides_connectibility_check,
		stk::mesh::EntityId * edge_ids,
		const std::vector < std::vector < unsigned > > &gold_edge_node_ids,
		bool *edges_connectibility_check)
{
	stk::mesh::EntityId element_id[1] = {1};
	stk::mesh::MetaData &meta = bulk.mesh_meta_data();
	stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part", topo);

	meta.commit();
	bulk.modification_begin();

	stk::mesh::Entity elem = stk::mesh::declare_element(bulk, elem_part, element_id[0], elem_node_ids);

	stk::mesh::EntityVector side_nodes;
	uint num_sides = topo.num_sides();
	stk::topology::rank_t sub_topo_rank = topo.side_rank();

	for(uint i = 0; i < num_sides; ++i)
	{
		stk::topology sub_topo = topo.side_topology(i);
		side_nodes.clear();

		stk::mesh::Entity side = bulk.declare_entity(sub_topo_rank, side_ids[i], meta.get_topology_root_part(sub_topo));

		for (uint j = 0; j < sub_topo.num_nodes(); ++j)
		{
			stk::mesh::Entity side_node = bulk.get_entity(stk::topology::NODE_RANK, gold_side_node_ids[i][j]);
			side_nodes.push_back(side_node);
			bulk.declare_relation(side, side_node, j);
		}

		std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(bulk, elem, sub_topo_rank, side_nodes);

		if (sides_connectibility_check[i])
		{
		  EXPECT_NE(ordinalAndPermutation.first, stk::mesh::ConnectivityOrdinal::INVALID_CONNECTIVITY_ORDINAL);
		  EXPECT_NE(ordinalAndPermutation.second, stk::mesh::Permutation::INVALID_PERMUTATION);
		}
		else
		{
		  EXPECT_EQ(ordinalAndPermutation.first, stk::mesh::ConnectivityOrdinal::INVALID_CONNECTIVITY_ORDINAL);
		  EXPECT_EQ(ordinalAndPermutation.second, stk::mesh::Permutation::INVALID_PERMUTATION);
		}
	}

	if (edge_ids == NULL) {
		bulk.modification_end();
		return;
	}

	stk::mesh::EntityVector edge_nodes;
	uint num_edges = topo.num_edges();

	for(uint i = 0; i < num_edges; ++i)
	{
		edge_nodes.clear();

		stk::mesh::Entity edge = bulk.declare_entity(stk::topology::EDGE_RANK, edge_ids[i], meta.get_topology_root_part(topo.edge_topology()));

		for (uint j = 0; j < topo.edge_topology().num_nodes(); ++j)
		{
			stk::mesh::Entity edge_node = bulk.get_entity(stk::topology::NODE_RANK, gold_edge_node_ids[i][j]);
			edge_nodes.push_back(edge_node);
			bulk.declare_relation(edge, edge_node, j);
		}

		std::pair<stk::mesh::ConnectivityOrdinal, stk::mesh::Permutation> ordinalAndPermutation = get_ordinal_and_permutation(bulk, elem, stk::topology::EDGE_RANK, edge_nodes);

		if (edges_connectibility_check[i])
		{
		  EXPECT_NE(ordinalAndPermutation.first, stk::mesh::ConnectivityOrdinal::INVALID_CONNECTIVITY_ORDINAL);
		  EXPECT_NE(ordinalAndPermutation.second, stk::mesh::Permutation::INVALID_PERMUTATION);
		}
		else
		{
		  EXPECT_EQ(ordinalAndPermutation.first, stk::mesh::ConnectivityOrdinal::INVALID_CONNECTIVITY_ORDINAL);
		  EXPECT_EQ(ordinalAndPermutation.second, stk::mesh::Permutation::INVALID_PERMUTATION);
		}
	}

	bulk.modification_end();
}

TEST(FEMHelper, test_permutations_for_key_topologies)
{
	stk::ParallelMachine pm = MPI_COMM_WORLD;
	const int p_size = stk::parallel_machine_size(pm);

	if (p_size != 1) {
		return;
	}

	int spatial_dimension;
	stk::topology topo;

	const unsigned num_test_topologies = 10;
	stk::topology test_topologies[num_test_topologies] = {stk::topology::TRI_3_2D, stk::topology::QUAD_4_2D, stk::topology::SHELL_TRI_3,
	                                                      stk::topology::SHELL_QUAD_4, stk::topology::TET_4, stk::topology::PYRAMID_5,
	                                                      stk::topology::WEDGE_6, stk::topology::HEX_8, stk::topology::TRI_6_2D,
														  stk::topology::TET_10};

	// check that the permutations define the same sides
	for (size_t index = 0; index < num_test_topologies; ++index)
	{
		topo = test_topologies[index];
		std::cout << "topology " << topo << " topology rank " << topo.rank() << " topology dimension " << topo.dimension() << std::endl;
		if(topo.rank() == stk::topology::ELEMENT_RANK)
		{
			spatial_dimension = topo.dimension();

			stk::mesh::MetaData meta(spatial_dimension);
			stk::mesh::BulkData bulk(meta, pm);

			//specific to topology
			switch (topo)
			{
				case stk::topology::TRI_3_2D:
				{
					std::array <stk::mesh::EntityId, 3> elem_node_ids = { {1, 2, 3}};
					std::array <stk::mesh::EntityId, 3> side_ids = { {1, 2, 3}};

					std::array < std::array <unsigned, 2>, 3 > gold_side_node_ids_data = {{{{1,2}},{{3,2}},{{3,1}}}};
					std::vector < std::vector < unsigned > > gold_side_node_ids(3);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(2);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j]; //note that 2nd is inverted from the normal way you'd set this side up
						}
					}
					unsigned gold_side_permutations[3] = { 0, 1, 0 };  //the permutation of 1 is consistent with inversion of side 2 from above
					stk::mesh::EntityId * edge_ids = NULL;
					std::vector < std::vector < unsigned > > gold_edge_node_ids;
					unsigned *  gold_edge_permutations = NULL;
					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], edge_ids, gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, gold_edge_permutations);
					break;
				}
				case stk::topology::QUAD_4_2D:
				{
			        std::array <stk::mesh::EntityId, 4> elem_node_ids = {{1, 2, 3, 4}};
					std::array <stk::mesh::EntityId, 4> side_ids = { {1, 2, 3, 4}};

					//note that 3rd is inverted from the normal way you'd set this side up
					std::array < std::array <unsigned, 2>, 4 > gold_side_node_ids_data = {{ {{1,2}}, {{2,3}}, {{4,3}}, {{4,1}} }};
					std::vector < std::vector < unsigned > > gold_side_node_ids(4);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(2);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j]; //note that 2nd is inverted from the normal way you'd set this side up
						}
					}
					unsigned gold_side_permutations[4] = { 0, 0, 1, 0 };  //the permutation of 1 is consistent with inversion of side 3 from above
					stk::mesh::EntityId * edge_ids = NULL;
					std::vector < std::vector < unsigned > > gold_edge_node_ids;
					unsigned *  gold_edge_permutations = NULL;
					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], edge_ids, gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, gold_edge_permutations);

					break;
				}
				case stk::topology::SHELL_TRI_3:
				{
				    std::array <stk::mesh::EntityId, 3> elem_node_ids = {{1, 2, 3}};
					std::array <stk::mesh::EntityId, 2> side_ids = {{1, 2}};
					std::array <stk::mesh::EntityId, 3> edge_ids = {{1, 2, 3}};

					std::array < std::array <unsigned, 3>, 2 > gold_side_node_ids_data = {{ {{1,2,3}}, {{3,2,1}} }};
					std::vector < std::vector < unsigned > > gold_side_node_ids(2);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(3);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
							gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j]; //note that 2nd is rotated from the normal way you'd set this side up
						}
					}
					unsigned gold_side_permutations[2] = { 0, 1 };

					std::array < std::array <unsigned, 2>, 3 > gold_edge_node_ids_data = {{ {{1,2}}, {{3,2}}, {{3,1}} }};
					std::vector < std::vector < unsigned > > gold_edge_node_ids(3);
					for (unsigned i = 0; i < topo.num_edges(); ++i) {
						gold_edge_node_ids[i].resize(2);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
						}
					}
					unsigned gold_edge_permutations[4] = { 0, 1, 0 };

					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);

          break;
				}
				case stk::topology::SHELL_QUAD_4:
				{
				    std::array <stk::mesh::EntityId, 4> elem_node_ids = {{1, 2, 3, 4}};
					std::array <stk::mesh::EntityId, 2> side_ids = {{1, 2}};
					std::array <stk::mesh::EntityId, 4> edge_ids = {{1, 2, 3, 4}};

					std::array < std::array <unsigned, 4>, 2 > gold_side_node_ids_data = {{ {{1,2,3,4}}, {{4,3,2,1}} }};
					std::vector < std::vector < unsigned > > gold_side_node_ids(2);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(4);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
							gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j]; //note that 2nd is rotated from the normal way you'd set this side up
						}
					}
					unsigned gold_side_permutations[2] = { 0, 1 };

					std::array < std::array <unsigned, 2>, 4 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{4,3}}, {{4,1}} }};
					std::vector < std::vector < unsigned > > gold_edge_node_ids(4);
					for (unsigned i = 0; i < topo.num_edges(); ++i) {
						gold_edge_node_ids[i].resize(2);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
						}
					}
					unsigned gold_edge_permutations[4] = { 0, 0, 1, 0 };

					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);

          break;
				}
				case stk::topology::TET_4:
				{
				    std::array <stk::mesh::EntityId, 4> elem_node_ids = {{1, 2, 3, 4}};
					std::array <stk::mesh::EntityId, 4> side_ids = {{1, 2, 3, 4}};
					std::array <stk::mesh::EntityId, 6> edge_ids = {{1, 2, 3, 4, 5, 6}};

					std::array < std::array <unsigned, 3>, 4 > gold_side_node_ids_data = {{ {{1,2,4}}, {{2,3,4}}, {{1,3,4}}, {{3,2,1}} }};
					std::vector < std::vector < unsigned > > gold_side_node_ids(4);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(3);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
							gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
						}
					}
					unsigned gold_side_permutations[4] = { 0, 0, 3, 1 };

					std::array < std::array <unsigned, 2>, 6 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{3,1}}, {{1,4}}, {{4,2}}, {{3,4}} }};
					std::vector < std::vector < unsigned > > gold_edge_node_ids(6);
					for (unsigned i = 0; i < topo.num_edges(); ++i) {
						gold_edge_node_ids[i].resize(2);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
						}
					}
					unsigned gold_edge_permutations[6] = { 0, 0, 0, 0, 1, 0 };

					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);
          break;
				}
				case stk::topology::PYRAMID_5:
				{
				  std::array <stk::mesh::EntityId, 5> elem_node_ids = {{1, 2, 3, 4, 5}};
				  std::array <stk::mesh::EntityId, 5> side_ids = {{1, 2, 3, 4, 5}};
				  std::array <stk::mesh::EntityId, 8> edge_ids = {{1, 2, 3, 4, 5, 6, 7, 8}};

				  std::array < std::array <unsigned, 4>, 5 > gold_side_node_ids_data = {{ {{1,2,5,0}}, {{2,3,5,0}}, {{3,4,5,0}}, {{4,5,1,0}}, {{1,4,3,2}} }};
				  std::vector < std::vector < unsigned > > gold_side_node_ids(5);
				  for (unsigned i = 0; i < topo.num_sides(); ++i) {
				    gold_side_node_ids[i].resize(4);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
				      gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_side_permutations[5] = { 0, 0, 0, 3, 0 };

				  std::array < std::array <unsigned, 2>, 8 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{3,4}}, {{4,1}}, {{1,5}}, {{2,5}}, {{5,3}}, {{4,5}} }};
				  std::vector < std::vector < unsigned > > gold_edge_node_ids(8);
				  for (unsigned i = 0; i < topo.num_edges(); ++i) {
				    gold_edge_node_ids[i].resize(2);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
				      gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_edge_permutations[8] = { 0, 0, 0, 0, 0, 0, 1, 0 };

				  build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
				      &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);

				  break;
				}
				case stk::topology::WEDGE_6:
				{
				  std::array <stk::mesh::EntityId, 6> elem_node_ids = {{1, 2, 3, 4, 5, 6}};
				  std::array <stk::mesh::EntityId, 5> side_ids = {{1, 2, 3, 4, 5}};
				  std::array <stk::mesh::EntityId, 9> edge_ids = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

				  std::array < std::array <unsigned, 4>, 5 > gold_side_node_ids_data = {{ {{1,2,5,4}}, {{3,6,5,2}}, {{1,4,6,3}}, {{1,2,3,0}}, {{4,5,6,0}} }};
				  std::vector < std::vector < unsigned > > gold_side_node_ids(5);
				  for (unsigned i = 0; i < topo.num_sides(); ++i) {
				    gold_side_node_ids[i].resize(4);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
				      gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_side_permutations[5] = { 0, 1, 0, 3, 0 };

				  std::array < std::array <unsigned, 2>, 9 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{3,1}}, {{4,5}}, {{5,6}}, {{6,4}}, {{1,4}}, {{5,2}}, {{3,6}} }};
				  std::vector < std::vector < unsigned > > gold_edge_node_ids(9);
				  for (unsigned i = 0; i < topo.num_edges(); ++i) {
				    gold_edge_node_ids[i].resize(2);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
				      gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_edge_permutations[9] = { 0, 0, 0, 0, 0, 0, 0, 1, 0 };

				  build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
				      &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);
				  break;
				}
				case stk::topology::HEX_8:
				{
				  std::array <stk::mesh::EntityId, 8> elem_node_ids = {{1, 2, 3, 4, 5, 6, 7, 8}};
				  std::array <stk::mesh::EntityId, 6> side_ids = {{1, 2, 3, 4, 5, 6}};
				  std::array <stk::mesh::EntityId, 12> edge_ids = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}};

				  std::array < std::array <unsigned, 4>, 6 > gold_side_node_ids_data = {{ {{1,2,6,5}}, {{2,3,7,6}}, {{3,4,8,7}}, {{1,5,8,4}}, {{1,2,3,4}}, {{5, 6, 7, 8}} }};
				  std::vector < std::vector < unsigned > > gold_side_node_ids(6);
				  for (unsigned i = 0; i < topo.num_sides(); ++i) {
				    gold_side_node_ids[i].resize(4);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
				      gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_side_permutations[6] = { 0, 0, 0, 0, 4, 0 };

				  std::array < std::array <unsigned, 2>, 12 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{3,4}}, {{4,1}}, {{5,6}}, {{6,7}}, {{7,8}}, {{8,5}}, {{1,5}}, {{2,6}}, {{7,3}}, {{4,8}} }};
				  std::vector < std::vector < unsigned > > gold_edge_node_ids(12);
				  for (unsigned i = 0; i < topo.num_edges(); ++i) {
				    gold_edge_node_ids[i].resize(2);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
				      gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_edge_permutations[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 };

				  build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
				      &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);
				  break;
				}
				case stk::topology::TRI_6_2D:
				{
				  std::array <stk::mesh::EntityId, 6> elem_node_ids = {{1, 2, 3, 4, 5, 6}};
				  std::array <stk::mesh::EntityId, 3> side_ids = {{1, 2, 3}};

				  std::array < std::array <unsigned, 3>, 3 > gold_side_node_ids_data = {{ {{1,2,4}}, {{3,2,5}}, {{3,1,6}} }};
				  std::vector < std::vector < unsigned > > gold_side_node_ids(3);
				  for (unsigned i = 0; i < topo.num_sides(); ++i) {
				    gold_side_node_ids[i].resize(3);
				    for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
				      gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
				    }
				  }
				  unsigned gold_side_permutations[3] = { 0, 1, 0};
				  stk::mesh::EntityId * edge_ids = NULL;
				  std::vector < std::vector < unsigned > > gold_edge_node_ids;
				  unsigned *  gold_edge_permutations = NULL;

				  build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], edge_ids, gold_side_node_ids,
				      &gold_side_permutations[0], gold_edge_node_ids, gold_edge_permutations);
				  break;
				}
				case stk::topology::TET_10:
				{
				    std::array <stk::mesh::EntityId, 10> elem_node_ids = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
					std::array <stk::mesh::EntityId, 4> side_ids = {{1, 2, 3, 4}};
					std::array <stk::mesh::EntityId, 6> edge_ids = {{1, 2, 3, 4, 5, 6}};

					std::array < std::array <unsigned, 6>, 4 > gold_side_node_ids_data = {{ {{1,2,4, 5,9,8}},
					                                                                        {{4,2,3, 9,6,10}},
					                                                                        {{1,4,3, 8,10,7}},
					                                                                        {{1,3,2, 7,6,5}} }};

					std::vector < std::vector < unsigned > > gold_side_node_ids(4);
					for (unsigned i = 0; i < topo.num_sides(); ++i) {
						gold_side_node_ids[i].resize(6);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
						  gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
						}
					}
					unsigned gold_side_permutations[4] = { 0, 2, 0, 0 };

					std::array < std::array <unsigned, 3>, 6 > gold_edge_node_ids_data = {{ {{1,2, 5}}, {{2,3, 6}}, {{3,1, 7}},
					                                                                        {{1,4, 8}}, {{4,2, 9}}, {{3,4, 10}} }};
					std::vector < std::vector < unsigned > > gold_edge_node_ids(6);
					for (unsigned i = 0; i < topo.num_edges(); ++i) {
						gold_edge_node_ids[i].resize(3);
						for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
							gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
						}
					}
					unsigned gold_edge_permutations[6] = { 0, 0, 0, 0, 1, 0 };

					build_element_from_topology_verify_ordinals_and_permutations(bulk, topo, &elem_node_ids[0], &side_ids[0], &edge_ids[0], gold_side_node_ids,
																				 &gold_side_permutations[0], gold_edge_node_ids, &gold_edge_permutations[0]);
					break;
			    }
				default :
				{
					throw std::runtime_error("Invalid Topology\n");
				}
			}
		}
  }
}

TEST(FEMHelper, verify_connectibility_failure)
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;
  const int p_size = stk::parallel_machine_size(pm);

  if (p_size != 1) {
    return;
  }

  int spatial_dimension;
  stk::topology topo;

  const unsigned num_test_topologies = 2;
  stk::topology test_topologies[num_test_topologies] = { stk::topology::PYRAMID_5, stk::topology::TET_10};

  // check that the permutations define the same sides
  for (size_t index = 0; index < num_test_topologies; ++index)
  {
    topo = test_topologies[index];
    std::cout << "topology " << topo << " topology rank " << topo.rank() << " topology dimension " << topo.dimension() << std::endl;
    if(topo.rank() == stk::topology::ELEMENT_RANK)
    {
      spatial_dimension = topo.dimension();

      stk::mesh::MetaData meta(spatial_dimension);
      stk::mesh::BulkData bulk(meta, pm);

      //specific to topology
      switch (topo)
      {
      case stk::topology::PYRAMID_5:
      {
        std::array <stk::mesh::EntityId, 5> elem_node_ids = {{1, 2, 3, 4, 5}};
        std::array <stk::mesh::EntityId, 5> side_ids = {{1, 2, 3, 4, 5}};
        std::array <stk::mesh::EntityId, 8> edge_ids = {{1, 2, 3, 4, 5, 6, 7, 8}};

        std::array < std::array <unsigned, 4>, 5 > gold_side_node_ids_data = {{ {{1,2,5,0}}, {{2,3,5,0}}, {{3,4,5,0}}, {{4,5,1,0}}, {{1,3,4,2}} }};
        std::vector < std::vector < unsigned > > gold_side_node_ids(5);
        for (unsigned i = 0; i < topo.num_sides(); ++i) {
          gold_side_node_ids[i].resize(4);
          for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
            gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
          }
        }

        std::array < std::array <unsigned, 2>, 8 > gold_edge_node_ids_data = {{ {{1,2}}, {{2,3}}, {{3,4}}, {{4,1}}, {{1,5}}, {{2,5}}, {{5,3}}, {{4,5}} }};
        std::vector < std::vector < unsigned > > gold_edge_node_ids(8);
        for (unsigned i = 0; i < topo.num_edges(); ++i) {
          gold_edge_node_ids[i].resize(2);
          for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
            gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
          }
        }

        bool sides_connectibility_check[5] = { 1, 1, 1, 1, 0 };
        bool edges_connectibility_check[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
        verify_unbuildable_element(bulk, topo, &elem_node_ids[0], &side_ids[0], gold_side_node_ids, sides_connectibility_check,
            &edge_ids[0], gold_edge_node_ids, edges_connectibility_check);
        break;
      }
      case stk::topology::TET_10:
      {
        std::array <stk::mesh::EntityId, 10> elem_node_ids = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}};
        std::array <stk::mesh::EntityId, 4> side_ids = {{1, 2, 3, 4}};
        std::array <stk::mesh::EntityId, 6> edge_ids = {{1, 2, 3, 4, 5, 6}};

        std::array < std::array <unsigned, 6>, 4 > gold_side_node_ids_data = {{ {{1,2,4, 5,9,8}},
            {{4,2,3, 9,6,10}},
            {{1,4,3, 8,10,7}},
            {{1,3,2, 7,6,5}} }};

        std::vector < std::vector < unsigned > > gold_side_node_ids(4);
        for (unsigned i = 0; i < topo.num_sides(); ++i) {
          gold_side_node_ids[i].resize(6);
          for (unsigned j = 0; j < topo.sub_topology(stk::topology::FACE_RANK, i).num_nodes(); ++j) {
            gold_side_node_ids[i][j] = gold_side_node_ids_data[i][j];
          }
        }

        std::array < std::array <unsigned, 3>, 6 > gold_edge_node_ids_data = {{ {{1,2, 5}}, {{2,3, 6}}, {{3,1, 7}},
            {{1,4, 8}}, {{4,9, 2}}, {{3,4, 10}} }};
        std::vector < std::vector < unsigned > > gold_edge_node_ids(6);
        for (unsigned i = 0; i < topo.num_edges(); ++i) {
          gold_edge_node_ids[i].resize(3);
          for (unsigned j = 0; j < topo.sub_topology(stk::topology::EDGE_RANK, i).num_nodes(); ++j) {
            gold_edge_node_ids[i][j] = gold_edge_node_ids_data[i][j];
          }
        }

        bool sides_connectibility_check[4] = { 1, 1, 1, 1 };
        bool edges_connectibility_check[6] = { 1, 1, 1, 1, 0, 1 };
        verify_unbuildable_element(bulk, topo, &elem_node_ids[0], &side_ids[0], gold_side_node_ids, sides_connectibility_check,
            &edge_ids[0], gold_edge_node_ids, edges_connectibility_check);
        break;
      }
      default:
      {
        throw std::runtime_error("Invalid Topology\n");
      }
      }
    }
  }
}
