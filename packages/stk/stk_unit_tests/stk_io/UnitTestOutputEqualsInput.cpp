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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for rand, srand, RAND_MAX
#include <stk_io/IossBridge.hpp>        // for is_part_io_part
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for get_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for string
#include <vector>                       // for vector, etc
#include "gtest/gtest.h"                // for AssertHelper, ASSERT_TRUE
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator&
#include "stk_mesh/base/Types.hpp"      // for PartVector, BucketVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_mesh/base/Comm.hpp"
#include "stk_unit_test_utils/ioUtils.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "Ioss_Region.h"
#include "Ioss_ElementBlock.h"

namespace {

TEST(StkMeshIoBroker, outputEqualsInput)
{
    // A simple test for reading and writing an exodus file, to make sure
    // that the output file is the same as the input file.

    std::string input_filename = unitTestUtils::getOption("--input-mesh", "no-mesh-specified");

    if (input_filename != "no-mesh-specified")
    {
        stk::ParallelMachine comm = MPI_COMM_WORLD;
        stk::mesh::MetaData meta;
        stk::mesh::BulkData bulk(meta, comm);

        stk::io::StkMeshIoBroker ioBroker(comm);
        ioBroker.set_bulk_data(bulk);
        ioBroker.add_mesh_database(input_filename, stk::io::READ_MESH);
        ioBroker.create_input_mesh();
        ioBroker.populate_bulk_data();

        std::vector<size_t> counts;
        stk::mesh::comm_mesh_counts(bulk, counts);

        if (bulk.parallel_rank() == 0)
        {
            std::cerr << "input-mesh = " << input_filename << std::endl;
            std::cerr << "      num-nodes = " << counts[stk::topology::NODE_RANK]
                      << ", num-elems = " << counts[stk::topology::ELEM_RANK] << std::endl;
        }

        std::cerr << "writing mesh with same io-broker as input" << std::endl;
        std::string same_iobroker_output_filename = std::string("output-same-iobroker-")+input_filename;
        size_t output_file_index = ioBroker.create_output_mesh(same_iobroker_output_filename, stk::io::WRITE_RESULTS);
        ioBroker.write_output_mesh(output_file_index);

        std::string new_iobroker_output_filename = std::string("output-")+input_filename;
        std::cerr << "writing mesh with new io-broker" << std::endl;

        stk::io::StkMeshIoBroker new_iobroker(comm);

        new_iobroker.set_bulk_data(bulk);
        size_t new_iobroker_output_file_index = new_iobroker.create_output_mesh(new_iobroker_output_filename, stk::io::WRITE_RESULTS);
        new_iobroker.write_output_mesh(new_iobroker_output_file_index);

        //Ideally we could use exodiff to compare the two output meshes to
        //verify they are the same.
        //But exodiff won't say two mesh files are different if the difference is only
        //a block-ordering difference. So let's instead verify that the two ioss-regions
        //have the same element-block grouping-entities.
        const Ioss::ElementBlockContainer& input_broker_elem_blocks = ioBroker.get_output_io_region(output_file_index)->get_element_blocks();
        const Ioss::ElementBlockContainer& new_iobroker_elem_blocks = new_iobroker.get_output_io_region(new_iobroker_output_file_index)->get_element_blocks();

        EXPECT_EQ(input_broker_elem_blocks.size(), new_iobroker_elem_blocks.size());
        for(size_t i=0; i<input_broker_elem_blocks.size(); ++i)
        {
            std::string input_name = input_broker_elem_blocks[i]->name();
            std::string new_name = new_iobroker_elem_blocks[i]->name();
            EXPECT_EQ(input_name, new_name);
        }
    }
    else
    {
        std::cerr << "test StkMeshIoBroker.outputEqualsInput, no mesh specified" << std::endl;
    }
}

}

