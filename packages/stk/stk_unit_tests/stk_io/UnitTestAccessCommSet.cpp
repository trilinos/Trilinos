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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include "stk_io/StkMeshIoBroker.hpp"
#include "Ioss_Region.h"
#include "Ioss_CommSet.h"

TEST(UnitTestAccessCommSet, basicNodeComm)
{
  stk::ParallelMachine communicator = MPI_COMM_WORLD;

  stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
  stkMeshIoBroker.add_mesh_database("generated:1x1x4", stk::io::READ_MESH);
  stkMeshIoBroker.create_input_mesh();
  stkMeshIoBroker.populate_bulk_data();

  const stk::mesh::MetaData& stkmeta = stkMeshIoBroker.meta_data();
  const stk::mesh::BulkData& stkmesh = stkMeshIoBroker.bulk_data();
  int myProc = stkmesh.parallel_rank();
  int numProcs = stkmesh.parallel_size();

  if (numProcs != 2) {
    std::cout<<"UnitTestAccessCommSet.basicNodeComm only runs on 2 procs"<<std::endl;
    return;
  }

  const unsigned expected_num_nodes = 16;
  unsigned num_nodes = stk::mesh::count_selected_entities(stkmeta.universal_part(), stkmesh.buckets(stk::topology::NODE_RANK));
  EXPECT_EQ(expected_num_nodes, num_nodes);

  std::cout<<"proc "<<myProc<<", num nodes: "<<num_nodes<<std::endl;

  Ioss::Region& ioss_region = *stkMeshIoBroker.get_input_ioss_region();
  Ioss::CommSet* io_cs = ioss_region.get_commset("commset_node");
  EXPECT_FALSE(NULL == io_cs);

  const int expected_num_shared_nodes = 4;
  int num_shared_nodes = io_cs->get_field("entity_processor").raw_count();
  EXPECT_EQ(expected_num_shared_nodes, num_shared_nodes);
  
  std::cout<<"proc "<<myProc<<", num shared nodes: "<<num_shared_nodes<<std::endl;

  std::vector<int> entity_proc;
  io_cs->get_field_data("entity_processor", entity_proc);

//  std::cout<<"proc "<<myProc<<", entity_processor field data size: "<<entity_proc.size()<<std::endl;
//  for(int i=0; i<num_shared_nodes; ++i) {
//     std::cout<<"on proc "<<myProc<<", node "<<entity_proc[i*2]<<" is shared with proc "<<entity_proc[i*2+1]<<std::endl;
//  }
}
