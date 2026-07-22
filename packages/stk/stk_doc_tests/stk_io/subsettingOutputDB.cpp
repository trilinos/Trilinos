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

#include <gtest/gtest.h>                // for ASSERT_EQ, AssertHelper, etc
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <string>                       // for allocator, operator+, etc
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region
#include "Ioss_Utils.h"                 // for Utils
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"
namespace Ioss { class DatabaseIO; }

namespace {

TEST(StkMeshIoBrokerHowTo, subsetOutputDatabase)
{
  std::string resultsFilename = "subsetted.results";
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
    size_t index = stkIo.add_mesh_database(input_filename, "generated", stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();
    stkIo.populate_bulk_data();

    stk::mesh::MetaData &meta_data = stkIo.meta_data();
    const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();

    // ============================================================
    //+ EXAMPLE
    //+ Create a selector containing just the shell parts.
    stk::mesh::Selector shell_subset;
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];
      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
        shell_subset |= *part;
      }
    }

    // Create the output...
    size_t fh = stkIo.create_output_mesh(resultsFilename,
                                         stk::io::WRITE_RESULTS);

    //+ Specify that only the subset of parts selected by the
    //+ "shell_subset" selector will be on the output database.
    stkIo.set_subset_selector(fh, shell_subset);
    stkIo.write_output_mesh(fh);
    // Verification omitted...
    //-END
  }

  {
    // ==================================================
    // VERIFICATION
    // Verify output mesh has correct number of nodes and elements.
    // Note that the output mesh will contain all element blocks;
    // however, the non-shell element block will have zero elements.
    // This is due to the subset_selector subsetting the entities and
    // not the parts...
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", resultsFilename,
                                                       Ioss::READ_MODEL, communicator);
    Ioss::Region ioRegion(iossDb);


    // The output model should consist of the elements and nodes in the 6 shell blocks.
    size_t num_elements = ioRegion.get_property("element_count").get_int();
    size_t num_nodes    = ioRegion.get_property("node_count").get_int();

    // Calculate the expected number of nodes and elements.
    size_t expected_elements = 6 * num_elems_per_edge * num_elems_per_edge ;

    size_t num_nodes_per_edge = num_elems_per_edge+1;
    size_t expected_nodes = 6 * num_nodes_per_edge*num_nodes_per_edge;
    expected_nodes -= 12 * num_nodes_per_edge; // Nodes on each edge were double-counted in previous calculation.
    expected_nodes += 8; // Nodes on each corner were removed in previous calculation; add them back.

    ASSERT_EQ(num_elements, expected_elements);
    ASSERT_EQ(num_nodes, expected_nodes);
  }
  unlink(resultsFilename.c_str());
}

}
