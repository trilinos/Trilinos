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

#include <gtest/gtest.h>                // for AssertHelper, ASSERT_EQ, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/GetEntities.hpp>  // for get_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <string>                       // for allocator, operator+, etc
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_NodeSet.h"               // for NodeSet
#include "Ioss_Region.h"                // for NodeSetContainer, Region
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, readInitialConditionNodalSubset)
{
  //-BEGIN
  std::string dbFieldNameShell = "NodeBlock_1";
  std::string appFieldName = "temperature";
  size_t num_elems_per_edge = 9;

  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  {
    // ============================================================
    // INITIALIZATION
    //+ Create a generated mesh containg hexes and shells with a
    //+ single nodal variable -- NodeBlock_1
    std::string s_elems_per_edge = std::to_string(num_elems_per_edge);

    //+ Create a generated mesh containg hexes and shells.
    std::string input_filename = s_elems_per_edge + "x" +
        s_elems_per_edge + "x" +
        s_elems_per_edge;
    input_filename += "|shell:xyzXYZ|variables:nodal,1|times:1";

    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(input_filename, "generated", stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta_data = stkIo.meta_data();

    // Declare the nodal "temperature" field. Exists on all nodes.
    stk::mesh::Field<double> &temperature =
        stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, appFieldName,1);

    // "NodeBlock_1" is the name of the node field on the input mesh.
    stk::io::MeshField mf(temperature, dbFieldNameShell);
    double time = stkIo.get_input_ioss_region()->get_state_time(1);
    mf.set_read_time(time);

    const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];

      //+ It exists on all nodes in the mesh...
      stk::mesh::put_field_on_mesh(temperature, *part, nullptr);

      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
        //+ Temperature field exists on all nodes in the mesh,
        //+ but only initialize it on the shell nodes.
        mf.add_subset(*part);
      }
    }

    stkIo.populate_bulk_data();


    //+ Populate the fields with data from the input mesh.
    stkIo.read_input_field(mf);

    //-END
    // ============================================================
    //+ VERIFICATION
    //+ The value of the field on the first 729 elements should be 0.0;
    //+ The value of the field on the remaining elements should be sqrt(i+1)
    size_t num_nodes_per_edge = num_elems_per_edge+1;
    std::vector<stk::mesh::Entity> nodes;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK,
                            nodes);
    size_t all_nodes = num_nodes_per_edge * num_nodes_per_edge * num_nodes_per_edge;
    EXPECT_EQ(all_nodes, nodes.size());

    // Create a selector for the nodes that are attached to the shells.
    stk::mesh::Selector shell_subset;
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];
      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
        shell_subset |= *part;
      }
    }

    //+ Get all nodes attached and not attached to shells.
    std::vector<stk::mesh::Entity> shell_nodes;
    std::vector<stk::mesh::Entity> other_nodes;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, shell_subset, shell_nodes);
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::NODE_RANK, !shell_subset, other_nodes);

    // Calculate number of nodes on surface of the mesh...
    size_t expected_nodes = 6 * num_nodes_per_edge*num_nodes_per_edge;
    expected_nodes -= 12 * num_nodes_per_edge; // Nodes on each edge were double-counted in previous calculation.
    expected_nodes += 8; // Nodes on each corner were removed in previous calculation; add them back.
    EXPECT_EQ(expected_nodes, shell_nodes.size());
    EXPECT_EQ(all_nodes-expected_nodes, other_nodes.size());

    for(size_t i=0; i<other_nodes.size(); i++) {
      double *fieldDataForNode = stk::mesh::field_data(temperature, other_nodes[i]);
      EXPECT_DOUBLE_EQ(0.0, *fieldDataForNode);
    }

    for(size_t i=0; i<shell_nodes.size(); i++) {
      double *fieldDataForNode = stk::mesh::field_data(temperature, shell_nodes[i]);
      size_t id = stkIo.bulk_data().identifier(shell_nodes[i]);
      EXPECT_DOUBLE_EQ(sqrt(id), *fieldDataForNode);
    }
  }
}
}
