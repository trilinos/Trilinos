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
#include <stk_mesh/base/MeshBuilder.hpp>
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

TEST(StkMeshIoBrokerHowTo, readInitialConditionMultiSubset)
{
  std::string dbFieldNameShell = "ElementBlock_1";
  std::string appFieldName = "pressure";

  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  {
    //-BEGIN
    // ============================================================
    // INITIALIZATION
    //+ Create a generated mesh containg hexes and shells with a
    //+ single element variable -- pressure
    std::string input_filename = "9x9x9|shell:xyzXYZ|variables:element,1|times:1";

    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(input_filename, "generated", stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta_data = stkIo.meta_data();

    // Declare the element "pressure" field...
    stk::mesh::Field<double> &pressure =
        stkIo.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, appFieldName, 1);

    const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();
    for (size_t i=0; i < all_parts.size(); i++) {
      //+ Put the field on all element block parts...
      stk::mesh::put_field_on_mesh(pressure, *all_parts[i], nullptr);
    }

    // This commits BulkData and populates the coordinates, connectivity, mesh...
    stkIo.populate_bulk_data();

    double time = stkIo.get_input_ioss_region()->get_state_time(1);

    //+ Initialize the "pressure" field from mesh data on the shell parts on demand...
    for (size_t i=0; i < all_parts.size(); i++) {
      stk::topology topo = all_parts[i]->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {

        stk::io::MeshField mf(pressure, dbFieldNameShell);
        mf.set_read_time(time);
        mf.add_subset(*all_parts[i]);
        stkIo.add_input_field(mf);
      }
    }

    //+ Populate any other fields with data from the input mesh.
    //+ This would *not* know about the MeshFields above since
    //+ "add_input_field()" was not called...
    stkIo.read_defined_input_fields(time);


    //-END
    // ============================================================
    //+ VERIFICATION
    //+ The value of the field on the first 729 elements should be 0.0;
    //+ The value of the field on the remaining elements should be sqrt(i+1)
    std::vector<stk::mesh::Entity> elements;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEMENT_RANK,
                            elements);
    EXPECT_TRUE(elements.size() >= 729);

    for(size_t i=0; i<729; i++) {
      double *fieldDataForElement = stk::mesh::field_data(pressure, elements[i]);
      EXPECT_DOUBLE_EQ(0.0, *fieldDataForElement);
    }

    for(size_t i=729; i<elements.size(); i++) {
      double *fieldDataForElement = stk::mesh::field_data(pressure, elements[i]);
      EXPECT_DOUBLE_EQ(sqrt(i+1), *fieldDataForElement);
    }
  }
}
}
