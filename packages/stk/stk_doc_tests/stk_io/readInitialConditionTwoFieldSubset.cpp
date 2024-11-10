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

TEST(StkMeshIoBrokerHowTo, readInitialConditionTwoFieldSubset)
{
  //-BEGIN
  std::string dbFieldNameShell = "ElementBlock_1";
  std::string dbFieldNameOther = "ElementBlock_2";
  std::string appFieldName = "pressure";

  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  {
    // ============================================================
    // INITIALIZATION
    //+ Create a generated mesh containg hexes and shells with two
    //+ element variables -- ElementBlock_1 and ElementBlock_2
    std::string input_filename = "9x9x9|shell:xyzXYZ|variables:element,2|times:1";

    stk::io::StkMeshIoBroker stkIo(communicator);
    stkIo.add_mesh_database(input_filename, "generated", stk::io::READ_MESH);
    stkIo.create_input_mesh();

    stk::mesh::MetaData &meta_data = stkIo.meta_data();

    // Declare the element "pressure" field...
    stk::mesh::Field<double> &pressure =
        stkIo.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, appFieldName,1);

    stk::io::MeshField mf_shell(pressure, dbFieldNameShell);
    stk::io::MeshField mf_other(pressure, dbFieldNameOther);

    const stk::mesh::PartVector &all_parts = meta_data.get_mesh_parts();
    for (size_t i=0; i < all_parts.size(); i++) {
      const stk::mesh::Part *part = all_parts[i];

      //+ Put the field on all element block parts...
      stk::mesh::put_field_on_mesh(pressure, *part, nullptr);

      stk::topology topo = part->topology();
      if (topo == stk::topology::SHELL_QUAD_4) {
        //+ The shell blocks will have the pressure field initialized
        //+ from the dbFieldNameShell database variable.
        mf_shell.add_subset(*part);
      }
      else {
        //+ The non-shell blocks will have the pressure field initialized
        //+ from the dbFieldNameOther database variable.
        mf_other.add_subset(*part);
      }
    }

    stkIo.add_input_field(mf_shell);
    stkIo.add_input_field(mf_other);
    stkIo.populate_bulk_data();

    double time = stkIo.get_input_ioss_region()->get_state_time(1);

    //+ Populate the fields with data from the input mesh.
    stkIo.read_defined_input_fields(time);

    //-END
    // ============================================================
    //+ VERIFICATION
    //+ The value of the field on all elements should be sqrt(i+1)
    std::vector<stk::mesh::Entity> elements;
    stk::mesh::get_entities(stkIo.bulk_data(), stk::topology::ELEMENT_RANK,
                            elements);
    EXPECT_TRUE(elements.size() >= 729);

    for(size_t i=0; i<elements.size(); i++) {
      double *fieldDataForElement = stk::mesh::field_data(pressure, elements[i]);
      EXPECT_DOUBLE_EQ(sqrt(i+1), *fieldDataForElement);
    }
  }
}
}
