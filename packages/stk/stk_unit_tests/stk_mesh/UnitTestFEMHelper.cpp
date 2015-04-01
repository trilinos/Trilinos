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

TEST(FEMHelper, test_permutations_for_all_topologies)
{
	stk::ParallelMachine pm = MPI_COMM_WORLD;
	//const int p_rank = stk::parallel_machine_rank(pm);
	const int p_size = stk::parallel_machine_size(pm);

	if (p_size != 1) {
		return;
	}

	int spatial_dimension;
	stk::mesh::EntityId element_id[1] = {1};

	// check that the permutations define the same sides
	for (stk::topology topo = stk::topology::BEGIN_TOPOLOGY; topo < stk::topology::END_TOPOLOGY; ++topo)
	{

		if (topo.defined_on_spatial_dimension(2u)) {spatial_dimension = 2;}
		else if (topo.defined_on_spatial_dimension(3u)) {spatial_dimension = 3;}

		stk::mesh::MetaData meta(spatial_dimension);
		stk::mesh::BulkData bulk(meta, pm);
		stk::mesh::Part &elem_part = meta.declare_part_with_topology("elem_part", topo);
		meta.commit();
		bulk.modification_begin();

		//specific to topology
		switch (topo) {
		case stk::topology::TRIANGLE_3: {
			stk::mesh::EntityId elem_node_ids[3] = { 1, 2, 3 };
			stk::mesh::declare_element(bulk, elem_part, element_id[0],
					elem_node_ids);
			stk::mesh::create_edges(bulk);
			stk::mesh::create_faces(bulk);
			break;
		}
		case stk::topology::TRIANGLE_4: {
			stk::mesh::EntityId elem_node_ids[4] = { 1, 2, 3, 4 };
			stk::mesh::declare_element(bulk, elem_part, element_id[0],
					elem_node_ids);
			stk::mesh::create_edges(bulk);
			stk::mesh::create_faces(bulk);
			break;
		}
		default :
		{

		}
		}
  //stk::mesh::EntityId element_ids [2] = {1, 2};
  //stk::mesh::EntityId elem_node_ids [][4] = {{1, 2, 3, 4}, {4, 3, 6, 5}};
  // Start with all entities on proc 0
  //std::vector<stk::mesh::Entity> elems;
  //if (p_rank == 0) {
  //stk::mesh::declare_element(bulk, elem_part , element_id[0], elem_node_ids );
  //elems.push_back(stk::mesh::declare_element(bulk, elem_part ,element_ids[1], elem_node_ids[1] ) );
  //}
     bulk.modification_end();

   }
}
