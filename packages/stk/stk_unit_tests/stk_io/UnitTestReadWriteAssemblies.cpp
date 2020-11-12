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
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <fstream>
#include "Assembly.hpp"

TEST_F(Assembly, readWriteAssembly_simple)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assemblyName("simpleAssembly");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assemblyPart = create_assembly(assemblyName);
  stk::mesh::Part& block1Part = create_io_part(partNames[0]);
  stk::mesh::Part& block2Part = create_io_part(partNames[1]);
  declare_subsets(assemblyPart, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);

  test_write_then_read_assemblies(1);
}

TEST_F(Assembly, readWriteAssembly_multiple)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
  const std::string assembly1Name("assembly1");
  const std::string assembly2Name("assembly2");
  const std::vector<std::string> partNames {"block_1", "block_2"};

  stk::mesh::Part& assembly1Part = create_assembly(assembly1Name);
  stk::mesh::Part& assembly2Part = create_assembly(assembly2Name);
  stk::mesh::Part& block1Part = create_io_part(partNames[0]);
  stk::mesh::Part& block2Part = create_io_part(partNames[1]);
  declare_subsets(assembly1Part, {&block1Part, &block2Part});
  declare_subsets(assembly2Part, {&block1Part, &block2Part});
  stk::io::fill_mesh("generated:2x2x2", get_bulk());
  move_element(1, block1Part, block2Part);

  test_write_then_read_assemblies(2);
}

TEST_F(Assembly, readWriteAssembly_hierarchy)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  std::string parentAssemblyName("myParentAssembly");
  std::string subAssembly1Name("mySubAssembly1");
  std::string subAssembly2Name("mySubAssembly2");
  std::vector<std::string> leafPartNames {"block_1", "block_2", "block_3"};

  stk::mesh::Part* parentAssemblyPart;
  stk::mesh::PartVector subAssemblyParts;
  std::tie(parentAssemblyPart, subAssemblyParts) = create_assembly_hierarchy(parentAssemblyName,
                                                     {subAssembly1Name, subAssembly2Name});
  stk::mesh::Part& block1Part = create_io_part(leafPartNames[0]);
  stk::mesh::Part& block2Part = create_io_part(leafPartNames[1]);
  stk::mesh::Part& block3Part = create_io_part(leafPartNames[2]);

  ASSERT_EQ(2u, subAssemblyParts.size());
  declare_subsets(*subAssemblyParts[0], {&block1Part, &block2Part});
  declare_subsets(*subAssemblyParts[1], {&block2Part, &block3Part});

  stk::io::fill_mesh("generated:2x2x2", get_bulk());

  move_element(1, block1Part, block2Part);
  move_element(2, block1Part, block3Part);

  test_write_then_read_assemblies(3);
}

TEST_F(Assembly, readWriteAssembly_deeperHierarchy)
{
  if(stk::parallel_machine_size(get_comm()) != 1) { return; }

  setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

  create_deep_assembly_hierarchy();

  stk::io::fill_mesh("generated:2x2x2", get_bulk());

  stk::mesh::Part& block1Part = *get_meta().get_part("block_1"); 
  stk::mesh::Part& block2Part = *get_meta().get_part("block_2"); 
  stk::mesh::Part& block3Part = *get_meta().get_part("block_3"); 

  move_element(1, block1Part, block2Part);
  move_element(2, block1Part, block3Part);

  test_write_then_read_assemblies(6);
}

