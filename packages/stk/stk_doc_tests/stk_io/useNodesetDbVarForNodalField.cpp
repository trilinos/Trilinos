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
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeSet.h"               // for NodeSet
#include "Ioss_Region.h"                // for NodeSetContainer, Region
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/FieldBLAS.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, useNodesetDbVarForNodalFields)
{
  std::string resultsFilename = "nodeset_fields.results";
  std::string dbFieldName = "temp";
  std::string appFieldName = "temperature";

  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) { return; }

  size_t num_elems_per_edge = 9;
  {
    //-BEGIN
    // ============================================================
    // INITIALIZATION
    std::string s_elems_per_edge = std::to_string(num_elems_per_edge);

    //+ Create a generated mesh containg hexes and shells.
    std::string input_filename = s_elems_per_edge + "x" +
        s_elems_per_edge + "x" +
        s_elems_per_edge + "|shell:xyzXYZ";

    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(input_filename, "generated",
                            stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta_data = stkIo.meta_data();
    stk::mesh::Field<double> &temperature =
        meta_data.declare_field<double>(stk::topology::NODE_RANK, appFieldName, 1);

    // ============================================================
    //+ Put the temperature field on the nodes of the shell parts.
    const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();
    stk::mesh::Selector shell_subset;
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];
      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
        stk::mesh::put_field_on_mesh(temperature, *part, nullptr);
      }
    }

    stkIo.populate_bulk_data();

    // Create the output...
    size_t fh = stkIo.create_output_mesh(resultsFilename, stk::io::WRITE_RESULTS);

    //+ The "temperature" field will be output on nodesets consisting
    //+ of the nodes of each part the field is defined on.
    stkIo.use_nodeset_for_sideset_nodes_fields(fh, true);
    stkIo.use_nodeset_for_block_nodes_fields(fh, true);
    stkIo.add_field(fh, temperature, dbFieldName);

    // Add three steps to the database
    // For each step, the value of the field is the value 'time'
    for (size_t i=0; i < 3; i++) {
      double time = i;

      stk::mesh::field_fill(time, temperature);

      stkIo.begin_output_step(fh, time);
      stkIo.write_defined_output_fields(fh);
      stkIo.end_output_step(fh);
    }
    // Verification omitted...
    //-END
  }

  {
    // ==================================================
    // VERIFICATION
    // Verify that the output mesh has 6 nodesets (input has 0)
    // corresponding to the nodes of the 6 shell blocks.
    // Each nodeset should have a single variable named dbFieldName

    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", resultsFilename,
                                                       Ioss::READ_MODEL, communicator);
    Ioss::Region ioRegion(iossDb);

    // Six nodesets should have been created -- 1 for each shell block in the mesh.
    const Ioss::NodeSetContainer &nsets = ioRegion.get_nodesets();
    size_t nset_count = 6;
    ASSERT_EQ(nset_count, nsets.size());

    for (size_t i=0; i < nsets.size(); i++) {
      const Ioss::NodeSet *nset = nsets[i];
      // Each nodeset should have a field named 'temp'
      ASSERT_TRUE(nset->field_exists(dbFieldName));
    }
  }
  unlink(resultsFilename.c_str());
}

}
