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
#include <stddef.h>                     // for size_t
#include <unistd.h>                     // for unlink
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <ostream>                      // for basic_ostream::operator<<
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker, etc
#include <stk_util/util/ParameterList.hpp>  // for ParameterList, etc
#include <stk_util/parallel/Parallel.hpp>
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_ElementTopology.h"       // for NameList
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_io/WriteMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace Ioss { class DatabaseIO; }

namespace
{

TEST(StkMeshIoBrokerHowTo, writeHistory)
{

  const std::string file_name = "History.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  stk::util::ParameterList params;

  {
    // ========================================================================
    // INITIALIZATION...
    // Add some params to write and read...
    params.set_param("PI", 3.14159);  // Double
    params.set_param("Answer", 42);   // Integer

    std::vector<double> my_vector;
    my_vector.push_back(2.78);
    my_vector.push_back(5.30);
    my_vector.push_back(6.21);
    params.set_param("some_doubles", my_vector);   // Vector of doubles...

    std::vector<int> ages;
    ages.push_back(55);
    ages.push_back(49);
    ages.push_back(21);
    ages.push_back(19);

    params.set_param("Ages", ages);   // Vector of integers...
  }

  {
    // ========================================================================
    // EXAMPLE USAGE...
    // Begin use of stk io history file...
    stk::io::StkMeshIoBroker stkIo(communicator);

    //-BEGIN
    //+ Define the heartbeat output and the format (BINARY)
    size_t hb = stkIo.add_heartbeat_output(file_name, stk::io::BINARY);
    //-END

    stk::util::ParameterMapType::const_iterator i = params.begin();
    stk::util::ParameterMapType::const_iterator iend = params.end();
    for (; i != iend; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &parameter = params.get_param(parameterName);

      // Tell history database which global variables should be output at each step...
      stkIo.add_heartbeat_global(hb, parameterName, parameter);
    }

    // Now output the global variables...
    int timestep_count = 2;
    double time = 0.0;
    for (int step=1; step <= timestep_count; step++) {
      time = step;
      stkIo.process_heartbeat_output(hb, step, time);
    }
  }

  {
    // ========================================================================
    // VERIFICATION:
    Ioss::DatabaseIO *iossDb = Ioss::IOFactory::create("exodus", file_name, Ioss::READ_MODEL, communicator);
    Ioss::Region region(iossDb);

    // Number steps
    EXPECT_EQ(region.get_property("state_count").get_int(), 2);
    region.begin_state(2);

    Ioss::NameList fields;
    region.field_describe(Ioss::Field::REDUCTION, &fields);
    EXPECT_EQ(fields.size(), 4u);

    std::vector<double> values;
    region.get_field_data("PI", values);
    EXPECT_NEAR(values[0], 3.14159, 1.0e-6);
    EXPECT_EQ(1u, values.size());

    std::vector<int> ages;
    region.get_field_data("Ages", ages);
    EXPECT_EQ(ages.size(), 4u);
    EXPECT_EQ(ages[0], 55);
    EXPECT_EQ(ages[1], 49);
    EXPECT_EQ(ages[2], 21);
    EXPECT_EQ(ages[3], 19);
  }

  // ========================================================================
  // CLEANUP:
  unlink(file_name.c_str());
}

std::string getElementFieldName()
{
  return "SomeElementField";
}

std::string getNodalFieldName()
{
  return "SomeNodalField";
}

double initialValue()
{
  return -1.0;
}

void setUpMeshWithFieldOnBlock1(stk::mesh::BulkData& bulk, stk::mesh::Field<double>& field1, stk::mesh::Field<double>& field2)
{
  stk::io::StkMeshIoBroker stkIo(bulk.parallel());
  stkIo.set_bulk_data(bulk);
  stkIo.add_mesh_database("generated:1x1x2", stk::io::READ_MESH);
  stkIo.create_input_mesh();

  double init = initialValue();
  stk::mesh::Part* block1 = bulk.mesh_meta_data().get_part("block_1");
  stk::mesh::put_field_on_mesh(field1, *block1, &init);
  stk::mesh::put_field_on_mesh(field2, *block1, &init);

  stkIo.add_all_mesh_fields_as_input_fields();
  stkIo.populate_bulk_data();
}

void putEntityIntoHistoryPart(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank, stk::mesh::EntityId id, stk::mesh::Part& historyPart)
{
  stk::mesh::Entity entity = bulk.get_entity(rank, id);
  bulk.modification_begin();
  if(bulk.is_valid(entity) && bulk.bucket(entity).owned())
  {
    bulk.change_entity_parts(entity, stk::mesh::ConstPartVector{&historyPart});
  }
  bulk.modification_end();
}

double getValue(stk::mesh::EntityId id, int stepNum)
{
  return 100*stepNum + id;
}

void writeHistoryFile(const std::string& historyFilename, stk::mesh::BulkData& bulk,
                      stk::mesh::Field<double>& elementField, stk::mesh::Field<double>& nodalField,
                      stk::mesh::Part& elementHistoryPart, stk::mesh::Part& nodeHistoryPart)
{
  stk::mesh::Selector subset = elementHistoryPart | nodeHistoryPart;
  stk::io::StkMeshIoBroker outStkIo;
  outStkIo.set_bulk_data(bulk);

  size_t outputFileIndex = outStkIo.create_output_mesh(historyFilename, stk::io::WRITE_RESULTS);
  outStkIo.set_subset_selector(outputFileIndex, subset);
  outStkIo.add_field(outputFileIndex, elementField);
  outStkIo.add_field(outputFileIndex, nodalField);
  outStkIo.write_output_mesh(outputFileIndex);

  stk::mesh::EntityVector elementEntities;
  stk::mesh::get_entities(bulk, stk::topology::ELEM_RANK, elementHistoryPart & bulk.mesh_meta_data().locally_owned_part(), elementEntities);

  stk::mesh::EntityVector nodeEntities;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, nodeHistoryPart & bulk.mesh_meta_data().locally_owned_part(), nodeEntities);

  for(int step = 1; step < 10; step++)
  {
    for(stk::mesh::Entity entity : elementEntities)
    {
      double *data = stk::mesh::field_data(elementField, entity);
      *data = getValue(bulk.identifier(entity), step);
    }

    for(stk::mesh::Entity entity : nodeEntities)
    {
      double *data = stk::mesh::field_data(nodalField, entity);
      *data = getValue(bulk.identifier(entity), step);
    }

    double time = 0.4*step;
    outStkIo.begin_output_step(outputFileIndex, time);
    outStkIo.write_defined_output_fields(outputFileIndex);
    outStkIo.end_output_step(outputFileIndex);
  }
}

void verify_data(stk::mesh::BulkData& bulk, stk::mesh::Field<double>& field, stk::mesh::EntityRank rank, int stepNum)
{
  stk::mesh::EntityVector entities;
  stk::mesh::get_entities(bulk, rank, bulk.mesh_meta_data().locally_owned_part(), entities);

  for(stk::mesh::Entity entity : entities)
  {
    double *data = stk::mesh::field_data(field, entity);
    double goldValue = getValue(bulk.identifier(entity), stepNum);
    if(*data != initialValue())
    {
      EXPECT_EQ(goldValue, *data);
    }
  }
}

void verifyHistoryFileOutput(const std::string& filename)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(MPI_COMM_WORLD).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  int numSteps = 0;
  double maxTime = 0;

  stk::io::fill_mesh_save_step_info(filename, *bulk, numSteps, maxTime);

  stk::mesh::Field<double>* field1 = static_cast<stk::mesh::Field<double>*>(meta.get_field(stk::topology::ELEM_RANK, getElementFieldName()));
  ASSERT_TRUE(field1!=nullptr);
  verify_data(*bulk, *field1, stk::topology::ELEM_RANK, numSteps);

  stk::mesh::Field<double>* field2 = static_cast<stk::mesh::Field<double>*>(meta.get_field(stk::topology::NODE_RANK, getNodalFieldName()));
  ASSERT_TRUE(field2!=nullptr);
  verify_data(*bulk, *field2, stk::topology::NODE_RANK, numSteps);
}

TEST(StkMeshIoBrokerHowTo, writeHistoryOfElementAndNode)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm)<=2)
  {
    stk::mesh::EntityId elementId = 1;
    stk::mesh::EntityId nodeId = 12;

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    std::shared_ptr<stk::mesh::BulkData> bulk = builder.create();
    stk::mesh::MetaData& meta = bulk->mesh_meta_data();

    stk::mesh::Field<double>& elemField = meta.declare_field<double>(stk::topology::ELEM_RANK, getElementFieldName());
    stk::mesh::Field<double>& nodalField = meta.declare_field<double>(stk::topology::NODE_RANK, getNodalFieldName());

    stk::mesh::Part& elementHistoryPart = meta.declare_part("ElementHistoryOutput", stk::topology::ELEM_RANK);
    stk::mesh::Part& nodeHistoryPart = meta.declare_part("NodeHistoryOutput", stk::topology::NODE_RANK);


    setUpMeshWithFieldOnBlock1(*bulk, elemField, nodalField);
    putEntityIntoHistoryPart(*bulk, stk::topology::ELEM_RANK, elementId, elementHistoryPart);
    putEntityIntoHistoryPart(*bulk, stk::topology::NODE_RANK, nodeId, nodeHistoryPart);

    const std::string historyFilename = "history.e";
    writeHistoryFile(historyFilename, *bulk, elemField, nodalField, elementHistoryPart, nodeHistoryPart);

    verifyHistoryFileOutput(historyFilename);
  }
}


TEST(StkMeshIoBrokerHowTo, writeEmptyHistory)
{

  const std::string file_name = "EmptyHistory.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(communicator);
  if (numProcs != 1) {
    return;
  }

  stk::util::ParameterList params;

  {
    // ========================================================================
    // EXAMPLE USAGE...
    // Begin use of stk io history file...
    stk::io::StkMeshIoBroker stkIo(communicator);

    //-BEGIN
    //+ Define the heartbeat output and the format (BINARY)
    size_t hb = stkIo.add_heartbeat_output(file_name, stk::io::BINARY);
    //-END

    stkIo.begin_define_transient_for_heartbeat(hb);

    stk::util::ParameterMapType::const_iterator i = params.begin();
    stk::util::ParameterMapType::const_iterator iend = params.end();
    for (; i != iend; ++i) {
      const std::string parameterName = (*i).first;
      stk::util::Parameter &parameter = params.get_param(parameterName);

      // Tell history database which global variables should be output at each step...
      stkIo.add_heartbeat_global(hb, parameterName, parameter);
    }

    stkIo.end_define_transient_for_heartbeat(hb);

    // Now output the global variables...
    int timestep_count = 2;
    double time = 0.0;
    for (int step=1; step <= timestep_count; step++) {
      time = step;
      EXPECT_NO_THROW(stkIo.process_heartbeat_output(hb, step, time));
    }
  }

  // ========================================================================
  // CLEANUP:
  unlink(file_name.c_str());
}

}
