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
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <algorithm>                    // for sort
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_io/StkMeshIoBroker.hpp"

//====================
extern int gl_argc;
extern char** gl_argv;

namespace {
//helper-functions for the unit-test------------------------------------------------------------

void verify_mesh_is_empty(const stk::mesh::BulkData& mesh)
{
  std::vector<unsigned> counts(mesh.mesh_meta_data().entity_rank_count());
  stk::mesh::count_entities(mesh.mesh_meta_data().universal_part(), mesh, counts);
  for(size_t i=0; i<counts.size(); ++i) {
    EXPECT_EQ(0u, counts[i]);
  }
}

//void create_coord_field(stk::mesh::MetaData& meta)
//{
//  stk::mesh::Field<double,stk::mesh::Cartesian>& coord_field = meta.declare_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");
//  stk::mesh::put_field(coord_field, meta.universal_part(), 3);
//  meta.set_coordinate_field(&coord_field);
//}
//
//void create_hex_part(stk::mesh::MetaData& meta)
//{
//  meta.declare_part_with_topology("block_1", stk::topology::HEX_8);
//  meta.commit();
//}

void create_1_hex_element(stk::mesh::BulkData& mesh)
{
  verify_mesh_is_empty(mesh);

  stk::mesh::Part& hex_block = *mesh.mesh_meta_data().get_part("block_1");

  stk::mesh::EntityId elem_id = 1;
  stk::mesh::EntityId node_ids[] = {1, 2, 3, 4, 5, 6, 7, 8};

  mesh.modification_begin();
  stk::mesh::declare_element(mesh, hex_block, elem_id, node_ids);
  mesh.modification_end();
}

void read_mesh(const std::string& filename, stk::mesh::BulkData& mesh)
{
  stk::io::StkMeshIoBroker reader(mesh.parallel());
  reader.set_bulk_data(mesh);
  reader.add_mesh_database(filename, stk::io::READ_MESH);
  reader.create_input_mesh();
  reader.populate_bulk_data();
}

void write_mesh(const std::string& filename, stk::mesh::BulkData& mesh)
{
  stk::io::StkMeshIoBroker writer(mesh.parallel());
  writer.set_bulk_data(mesh);
  size_t output_handle = writer.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  writer.write_output_mesh(output_handle);
}

//end of helper-functions for the unit-test------------------------------------------------------------
}

TEST(ReplaceBulkData, read_replace_write)
{
    stk::ParallelMachine pm = MPI_COMM_WORLD;
    const int p_size = stk::parallel_machine_size(pm);
    if (p_size > 1) {
      return;
    }

//    const unsigned spatialDim = 3;
    stk::mesh::MetaData* meta = new stk::mesh::MetaData;
    stk::mesh::BulkData* bulk = new stk::mesh::BulkData(*meta, pm);

    std::string filename("generated:1x1x1");
    read_mesh(filename, *bulk);

    //replace BulkData:
    delete bulk;
    bulk = new stk::mesh::BulkData(*meta, pm);

    create_1_hex_element(*bulk);

    //write out to exodus file
    std::string out_filename("replace_out.exo");
    write_mesh(out_filename, *bulk);

    delete bulk;
    delete meta;
    unlink(out_filename.c_str());
}

