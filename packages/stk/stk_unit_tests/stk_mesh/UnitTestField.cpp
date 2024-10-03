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

#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, etc
#include "mpi.h"                        // for ompi_communicator_t, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/MeshField.hpp"         // for MeshField
#include "stk_io/WriteMesh.hpp"
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/FieldDataManager.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector, etc
#include "stk_mesh/base/MeshBuilder.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/GenerateALefRAMesh.hpp"
#include <stk_unit_test_utils/GetMeshSpec.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <iostream>                     // for ostream, operator<<, etc
#include <stddef.h>                     // for size_t, NULL
#include <stdexcept>                    // for runtime_error
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_mesh/base/NgpUtils.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include "ngp/NgpUnitTestUtils.hpp"
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "stk_util/util/AdjustForAlignment.hpp"
#include <string>                       // for string, operator==, etc
#include <unistd.h>                     // for unlink
#include <vector>                       // for vector

namespace Ioss { class DatabaseIO; }

namespace {

const stk::topology::rank_t NODE_RANK = stk::topology::NODE_RANK;
using stk::unit_test_util::build_mesh;

TEST(UnitTestField, testFieldMaxSize)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for some test fields
  typedef stk::mesh::Field<double> rank_zero_field;
  typedef stk::mesh::Field<double> rank_one_field;
  typedef stk::mesh::Field<double> rank_two_field;

  const std::string name0("test_field_0");
  const std::string name1("test_field_1");
  const std::string name2("test_field_2");

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();

  rank_zero_field & f0 = meta_data.declare_field<double>(NODE_RANK, name0);
  rank_one_field  & f1 = meta_data.declare_field<double>(NODE_RANK, name1);
  rank_two_field  & f2 = meta_data.declare_field<double>(NODE_RANK, name2);

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK);
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK);
  stk::mesh::Part & p2 = meta_data.declare_part("P2", NODE_RANK);

  stk::mesh::put_field_on_mesh(f0, p0, nullptr);
  stk::mesh::put_field_on_mesh(f1, p1, 10, nullptr);
  stk::mesh::put_field_on_mesh(f2, p2, 10, 20, nullptr);

  meta_data.commit();

  EXPECT_EQ( f0.max_size(), 1u );

  EXPECT_EQ( f1.max_size(), 10u );

  EXPECT_EQ( f2.max_size(), 200u );
}

TEST(UnitTestField, fieldDataAccess_rankMustMatch)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { return; }

  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;
  using MyField = stk::mesh::Field<double>;
  MyField& nodalField = meta.declare_field<double>(NODE_RANK, "nodal_field");

  stk::mesh::put_field_on_mesh(nodalField, meta.universal_part(), nullptr);

  stk::io::fill_mesh("generated:2x2x2|sideset:xXyYzZ", bulk);

  stk::mesh::EntityVector nodes;
  stk::mesh::get_entities(bulk, stk::topology::NODE_RANK, meta.universal_part(), nodes);
  stk::mesh::EntityVector faces;
  stk::mesh::get_entities(bulk, stk::topology::FACE_RANK, meta.universal_part(), faces);

  ASSERT_TRUE(!nodes.empty());
  ASSERT_TRUE(!faces.empty());

  EXPECT_NO_THROW(stk::mesh::field_data(nodalField, nodes[0]));
#ifndef NDEBUG
  EXPECT_THROW(stk::mesh::field_data(nodalField, faces[0]), std::logic_error);
#endif
}

TEST(UnitTestField, testFieldWithSelector)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for test field
  typedef stk::mesh::Field<double> rank_zero_field;

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;

  rank_zero_field  & f0 = meta_data.declare_field<double>( NODE_RANK, name0 );

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK );

  stk::mesh::Selector select_p0 = p0;
  std::cout <<"select_p0: "<< select_p0 << std::endl;

  stk::mesh::put_field_on_mesh( f0 , select_p0 , nullptr);

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare 10 nodes on each part

  for ( unsigned i = 1 ; i < 11 ; ++i )
    bulk_data.declare_node(i, std::vector< stk::mesh::Part * >(1, &p0));

  for ( unsigned i = 11 ; i < 21 ; ++i )
    bulk_data.declare_node(i, std::vector< stk::mesh::Part * >(1, &p1));

  const stk::mesh::BucketVector & node_buckets =
      bulk_data.buckets( NODE_RANK );

  unsigned num = stk::mesh::count_selected_entities(select_p0, node_buckets);

  ASSERT_EQ( 10u, num );

  stk::mesh::Selector select_f0 = stk::mesh::selectField(f0);

  std::cout <<"select_f0: "<< select_f0 << std::endl;

  unsigned num_f0 = stk::mesh::count_selected_entities(select_f0, node_buckets);
  ASSERT_EQ(10u, num_f0);

  stk::mesh::BucketVector const& f0_buckets = bulk_data.get_buckets(NODE_RANK, select_p0);
  unsigned num_buckets = f0_buckets.size();
  ASSERT_EQ(1u, num_buckets);

  for(stk::mesh::Bucket* b : f0_buckets) {
    unsigned f0_size = field_bytes_per_entity(f0, *b);
    ASSERT_EQ(8u, f0_size);
  }
}

TEST(UnitTestField, testFieldWithSelectorAnd)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  typedef stk::mesh::Field<double> rank_one_field;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;

  rank_one_field & f0 = meta_data.declare_field<double>( stk::topology::ELEMENT_RANK, name0 );

  stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  stk::mesh::Part & elements = meta_data.declare_part("Elements", elem_rank);
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );
  stk::mesh::Part & tet4s = meta_data.declare_part("Tet4", elem_rank );

  stk::mesh::Selector elem_hex_selector = elements & hex8s;
  stk::mesh::Selector elem_tet_selector = elements & tet4s;
  std::cout <<"elem_hex_selector: "<< elem_hex_selector << std::endl;
  std::cout <<"elem_tet_selector: "<< elem_tet_selector << std::endl;

  stk::mesh::put_field_on_mesh( f0 , elem_hex_selector, 8u , nullptr);
  stk::mesh::put_field_on_mesh( f0 , elem_tet_selector, 4u , nullptr);

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  // Declare 10 elements on each part

  stk::mesh::PartVector parts;
  parts.push_back(&elements);
  parts.push_back(&hex8s);

  for ( unsigned i = 1 ; i < 11 ; ++i )
    bulk_data.declare_element(i, parts);

  parts.clear();
  parts.push_back(&elements);
  parts.push_back(&tet4s);

  for ( unsigned i = 11 ; i < 21 ; ++i )
    bulk_data.declare_element(i, parts);

  {
    stk::mesh::BucketVector const& f0_buckets = bulk_data.get_buckets(elem_rank, elem_hex_selector);

    for(stk::mesh::Bucket* b : f0_buckets) {
      unsigned f0_size = field_bytes_per_entity(f0, *b);
      ASSERT_EQ(64u, f0_size);
    }
  }

  {
    stk::mesh::BucketVector const& f0_buckets = bulk_data.get_buckets(elem_rank, elem_tet_selector);

    for(stk::mesh::Bucket* b : f0_buckets) {
      unsigned f0_size = field_bytes_per_entity(f0, *b);
      ASSERT_EQ(32u, f0_size);
    }
  }
}


TEST(UnitTestField, testFieldWithSelectorInvalid)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  typedef stk::mesh::Field<double> rank_one_field;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatial_dimension, pm);
  stk::mesh::MetaData& meta_data = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk_data = *bulkPtr;

  rank_one_field & f0 = meta_data.declare_field<double>( stk::topology::ELEMENT_RANK, name0 );

  stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );

  stk::mesh::Part & universal_part = meta_data.universal_part();
  stk::mesh::Selector elem_hexA_selector = hex8s;
  stk::mesh::Selector elem_hexB_selector = universal_part & hex8s;

  std::cout <<"elem_hexA_selector: "<< elem_hexA_selector << std::endl;
  std::cout <<"elem_hexB_selector: "<< elem_hexB_selector << std::endl;

  stk::mesh::put_field_on_mesh( f0 , elem_hexA_selector, 8u , nullptr);
  ASSERT_THROW(
        stk::mesh::put_field_on_mesh( f0 , elem_hexA_selector, 4u , nullptr),
        std::runtime_error
        );
  stk::mesh::put_field_on_mesh( f0 , elem_hexB_selector, 4u , nullptr);

  stk::mesh::print( oss , "  " , f0 );

  meta_data.commit();

  bulk_data.modification_begin();

  stk::mesh::PartVector parts;
  parts.push_back(&hex8s);
  ASSERT_THROW(
        bulk_data.declare_element(1, parts),
        std::runtime_error
        );

}

TEST(UnitTestField, writeFieldsWithSameName)
{
  std::string mesh_name = "mesh_fields_with_same_name.e";
  MPI_Comm communicator = MPI_COMM_WORLD;
  const std::string fieldName = "MyFieldForElementsAndNodes";
  const double nodeInitialValue = 1.25;
  const double elemInitialValue = 8.39;
  const double time = 1.0;

  // Create the mesh with fields with the same name
  {
    stk::io::StkMeshIoBroker stkIo(communicator);

    const std::string generatedFileName = "generated:4x4x16";
    size_t index = stkIo.add_mesh_database(generatedFileName, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();

    stk::mesh::Field<double> &nodeField = stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, fieldName, 1);
    stk::mesh::put_field_on_mesh(nodeField, stkIo.meta_data().universal_part(), &nodeInitialValue);

    stk::mesh::Field<double> &elemField = stkIo.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, fieldName, 1);
    stk::mesh::put_field_on_mesh(elemField, stkIo.meta_data().universal_part(), &elemInitialValue);

    stkIo.populate_bulk_data();

    size_t fh = stkIo.create_output_mesh(mesh_name, stk::io::WRITE_RESULTS);
    stkIo.write_output_mesh(fh);
    stkIo.add_field(fh, nodeField);
    stkIo.add_field(fh, elemField);
    stkIo.begin_output_step(fh, time);
    stkIo.write_defined_output_fields(fh);
    stkIo.end_output_step(fh);
  }

  // Verify that the fields were written out to disk properly
  {
    Ioss::DatabaseIO *resultsDb = Ioss::IOFactory::create("exodus", mesh_name, Ioss::READ_MODEL, communicator);
    Ioss::Region results(resultsDb);
    const int goldNumSteps = 1;
    EXPECT_EQ(goldNumSteps, results.get_property("state_count").get_int());
    // Should be 1 nodal field on database named "disp";
    Ioss::NodeBlock *nb = results.get_node_blocks()[0];
    const unsigned goldNumNodeFields = 1;
    EXPECT_EQ(goldNumNodeFields, nb->field_count(Ioss::Field::TRANSIENT));
    EXPECT_TRUE(nb->field_exists(fieldName));
    Ioss::ElementBlock *eb = results.get_element_blocks()[0];
    const unsigned goldNumElemFields = 1;
    EXPECT_EQ(goldNumElemFields, eb->field_count(Ioss::Field::TRANSIENT));
    EXPECT_TRUE(eb->field_exists(fieldName));

    const int step = 1;
    double db_time = results.begin_state(step);
    EXPECT_EQ(time, db_time);

    std::vector<double> nodeFieldData;
    nb->get_field_data(fieldName, nodeFieldData);
    for (size_t node = 0; node < nodeFieldData.size(); node++) {
      EXPECT_EQ(nodeInitialValue, nodeFieldData[node]);
    }

    std::vector<double> elemFieldData;
    eb->get_field_data(fieldName, elemFieldData);
    for (size_t elem = 0; elem < elemFieldData.size(); elem++) {
      EXPECT_EQ(elemInitialValue, elemFieldData[elem]);
    }

    results.end_state(step);
  }

  // Verify that we can read the mesh back into memory correctly
  {
    stk::io::StkMeshIoBroker stkIo(communicator);

    size_t index = stkIo.add_mesh_database(mesh_name, stk::io::READ_MESH);
    stkIo.set_active_mesh(index);
    stkIo.create_input_mesh();
    const double badInitialData = -1.2345;

    stk::mesh::Field<double> &nodeField = stkIo.meta_data().declare_field<double>(stk::topology::NODE_RANK, fieldName, 1);
    stk::mesh::put_field_on_mesh(nodeField, stkIo.meta_data().universal_part(), &badInitialData);

    stk::mesh::Field<double> &elemField = stkIo.meta_data().declare_field<double>(stk::topology::ELEMENT_RANK, fieldName, 1);
    stk::mesh::put_field_on_mesh(elemField, stkIo.meta_data().universal_part(), &badInitialData);

    stkIo.populate_bulk_data();
    stkIo.add_input_field(stk::io::MeshField(nodeField, fieldName));
    stkIo.add_input_field(stk::io::MeshField(elemField, fieldName));
    stkIo.read_defined_input_fields(time);
    stk::mesh::BulkData &mesh = stkIo.bulk_data();
    stk::mesh::MetaData &metaData = stkIo.meta_data();

    const stk::mesh::BucketVector &nodeBuckets = mesh.get_buckets(stk::topology::NODE_RANK, metaData.locally_owned_part());
    for (size_t bucket_i=0 ; bucket_i<nodeBuckets.size() ; ++bucket_i) {
      stk::mesh::Bucket &nodeBucket = *nodeBuckets[bucket_i];
      for (size_t node_i=0 ; node_i<nodeBucket.size() ; ++node_i) {
        double * nodeData = stk::mesh::field_data(nodeField,nodeBucket.bucket_id(),node_i);
        EXPECT_EQ(nodeInitialValue, *nodeData);
      }
    }

    const stk::mesh::BucketVector &elemBuckets = mesh.get_buckets(stk::topology::ELEM_RANK, metaData.locally_owned_part());
    for (size_t bucket_i=0 ; bucket_i<elemBuckets.size() ; ++bucket_i) {
      stk::mesh::Bucket &elemBucket = *elemBuckets[bucket_i];
      for (size_t elem_i=0 ; elem_i<elemBucket.size() ; ++elem_i) {
        double * elemData = stk::mesh::field_data(elemField,elemBucket.bucket_id(),elem_i);
        EXPECT_EQ(elemInitialValue, *elemData);
      }
    }

    // Test Field accessor functions:
    stk::mesh::Field<double> *myTemplatedField = stk::mesh::get_field_by_name<double>(fieldName, metaData);
    ASSERT_TRUE(myTemplatedField != NULL);
    EXPECT_TRUE( &nodeField == myTemplatedField);
    stk::mesh::FieldBase *myFieldBase = stk::mesh::get_field_by_name(fieldName, metaData);
    ASSERT_TRUE(myFieldBase != NULL);
    EXPECT_TRUE( &nodeField == myFieldBase);
  }

  stk::unit_test_util::delete_mesh(mesh_name);
}

//////////////////////

void write_mesh_with_fields(stk::io::StkMeshIoBroker &stkIo, const std::string& filename, stk::mesh::EntityRank rank, const std::vector<std::string>& field_names)
{
  size_t fh = stkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  for(size_t i=0;i<field_names.size();++i)
  {
    stk::mesh::FieldBase* field = stkIo.bulk_data().mesh_meta_data().get_field(rank, field_names[i]);
    ASSERT_TRUE(field!=nullptr);
    stkIo.add_field(fh, *field);
  }

  double timeValue = 1.0;
  stkIo.begin_output_step(fh, timeValue);
  stkIo.write_defined_output_fields(fh);
  stkIo.end_output_step(fh);
}

class SolutionCases
{
public:
  const int numSolnCases = 3;
  enum cases { STATIC, TRANSIENT, EIGEN, NUM_SOLUTION_CASES };
  std::vector<std::string> solnNames = { "static", "trans", "eigen" };

  SolutionCases() : fieldsPerSolutionCase(numSolnCases) {}
  ~SolutionCases() {}

  const std::vector<std::string> get_solution_case_names() const { return solnNames; }

  const std::string get_part_name_for_solution_case(size_t i) const { return solnNames[i]; }

  void add_static_field(const std::string& fieldname)
  {
    fieldsPerSolutionCase[STATIC].push_back(fieldname);
  }

  void add_eigen_field(const std::string& fieldname)
  {
    fieldsPerSolutionCase[EIGEN].push_back(fieldname);
  }

  void add_transient_field(const std::string& fieldname)
  {
    fieldsPerSolutionCase[TRANSIENT].push_back(fieldname);
  }

  size_t get_num_solution_cases() const { return NUM_SOLUTION_CASES; }

  const std::vector<std::string> get_fields_for_case(int i) const
  {
    return fieldsPerSolutionCase[i];
  }

  void set_up_additional_fields_on_mesh(stk::mesh::MetaData& meta,
                                        stk::mesh::Part &part,
                                        stk::mesh::EntityRank rank,
                                        const std::vector<std::string>& fieldNames)
  {
    for(size_t j = 0; j < fieldNames.size(); ++j)
    {
      stk::mesh::Field<double>& field = meta.declare_field<double>(rank, fieldNames[j], 1u);
      stk::mesh::Selector sel = part;
      stk::mesh::put_field_on_mesh(field, sel, &init_val);
    }
  }

private:
  std::vector<std::vector<std::string>> fieldsPerSolutionCase;
  double init_val = 1.0;
};

SolutionCases setup_solution_cases()
{
  SolutionCases solnCases;
  solnCases.add_static_field("displacement");
  solnCases.add_transient_field("displacement");
  solnCases.add_transient_field("acceleration");
  solnCases.add_eigen_field("displacement");
  return solnCases;
}

void setup_parts_associated_with_solution_cases(SolutionCases &solnCases, stk::mesh::MetaData& meta, stk::mesh::EntityRank rank, stk::mesh::PartVector &partVector)
{
  for(size_t i=0;i<solnCases.get_num_solution_cases();++i)
  {
    std::string partName = solnCases.get_part_name_for_solution_case(i);
    stk::mesh::Part* part = nullptr;
    if(rank!=stk::topology::ELEM_RANK)
      part = &meta.declare_part(partName, rank);
    else
      part = &meta.declare_part_with_topology(partName, stk::topology::HEX_8);

    stk::io::put_io_part_attribute(*part);
    partVector.push_back(part);
    solnCases.set_up_additional_fields_on_mesh(meta, *part, rank, solnCases.get_fields_for_case(i));
  }
}

void move_entities_into_solution_part(int soln_index, stk::mesh::BulkData& bulk, stk::mesh::PartVector& partVector, stk::mesh::EntityRank rank, const stk::mesh::EntityVector& nodes)
{
  std::vector<stk::mesh::PartVector> addPartsPerEntity(nodes.size(), {partVector[soln_index]});
  std::vector<stk::mesh::PartVector> removePartsPerEntity;

  if(soln_index==0)
  {
    removePartsPerEntity.assign(nodes.size(), stk::mesh::PartVector{});
  }
  else
  {
    removePartsPerEntity.assign(nodes.size(), stk::mesh::PartVector{partVector[soln_index-1]} );
  }

  bulk.batch_change_entity_parts(nodes, addPartsPerEntity, removePartsPerEntity);

  stk::mesh::Field<double> *dispField = bulk.mesh_meta_data().get_field<double>(rank, "displacement");
  for(size_t j=soln_index; j<nodes.size(); ++j)
  {
    double *data = stk::mesh::field_data(*dispField, nodes[j]);
    *data = static_cast<double>(soln_index);
  }
}

void verify_acceleration_is_not_on_entities(stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank)
{
  stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(rank, "acceleration");
  stk::mesh::Selector accelEntities = stk::mesh::selectField(*field);
  unsigned numAccelEntities = stk::mesh::count_selected_entities(accelEntities, bulk.buckets(rank));
  EXPECT_EQ(0u, numAccelEntities);
}

void verify_fields_are_on_entities(const std::string& filename, stk::mesh::EntityRank rank, const std::vector<std::string>& fieldnames, size_t goldNum)
{
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
  stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
  stk::mesh::BulkData& bulk = *bulkPtr;
  int stepNum = 1;
  double time = 1.0;
  stk::io::fill_mesh_save_step_info(filename, bulk, stepNum, time);

  for(const std::string& field_name : fieldnames)
  {
    stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(rank, field_name);
    STK_ThrowRequireWithSierraHelpMsg(field!=nullptr);
    stk::mesh::Selector selector = stk::mesh::selectField(*field) & meta.locally_owned_part();
    unsigned numAccelEntities = stk::mesh::count_selected_entities(selector, bulk.buckets(rank));
    EXPECT_EQ(goldNum, numAccelEntities);
  }
}

class FieldFixture : public stk::unit_test_util::MeshFixture
{
protected:
  void test_solution_case_with_rank(stk::mesh::EntityRank rank)
  {
    setup_empty_mesh(m_auraOption);
    SolutionCases solnCases = setup_solution_cases();

    stk::mesh::PartVector solutionCasePartVector;
    setup_parts_associated_with_solution_cases(solnCases, get_meta(), rank, solutionCasePartVector);

    stk::io::fill_mesh("generated:1x1x2", get_bulk());

    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(get_bulk());

    stk::mesh::EntityVector locallyOwnedEntities;
    stk::mesh::Selector selector = get_meta().locally_owned_part();
    stk::mesh::get_selected_entities(selector, get_bulk().buckets(rank), locallyOwnedEntities);

    for(size_t solnIndex=0;solnIndex<solnCases.get_num_solution_cases();++solnIndex)
    {
      move_entities_into_solution_part(solnIndex, get_bulk(), solutionCasePartVector, rank, locallyOwnedEntities);
      std::string filename = "junk-" + solnCases.get_solution_case_names()[solnIndex] + ".g";
      EXPECT_NO_THROW(write_mesh_with_fields(stkIo, filename, rank, solnCases.get_fields_for_case(solnIndex)));
      verify_fields_are_on_entities(filename, rank, solnCases.get_fields_for_case(solnIndex), locallyOwnedEntities.size());
      stk::unit_test_util::delete_mesh(filename);
    }

    verify_acceleration_is_not_on_entities(get_bulk(), rank);
  }

  stk::mesh::BulkData::AutomaticAuraOption m_auraOption = stk::mesh::BulkData::NO_AUTO_AURA;
};

TEST_F(FieldFixture, totalNgpFieldDataBytes)
{
  setup_empty_mesh(m_auraOption);
  stk::mesh::Part & partA = get_meta().declare_part("partA", stk::topology::ELEM_RANK);
  stk::mesh::Part & partB = get_meta().declare_part("partB", stk::topology::ELEM_RANK);
  stk::mesh::Field<double> &field = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "doubleField");
  const double scalarInitValue = 3.14;
  const double vectorInitValue[] = {1., 2., 3., 4., 5.};
  stk::mesh::put_field_on_mesh(field, partA, &scalarInitValue);
  stk::mesh::put_field_on_mesh(field, partB, 5, vectorInitValue);

  const int numElemsPerDim = 10;
  stk::io::fill_mesh(stk::unit_test_util::get_mesh_spec(numElemsPerDim), get_bulk());
  const int totalNumElements = numElemsPerDim * numElemsPerDim * numElemsPerDim;

  stk::mesh::EntityVector elements;
  stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK), elements);
  stk::mesh::EntityVector partAElements;
  stk::mesh::EntityVector partBElements;
  for (const stk::mesh::Entity & element : elements) {
    if (get_bulk().identifier(element) <= (totalNumElements/2)) {
      partAElements.push_back(element);
    }
    else {
      partBElements.push_back(element);
    }
  }
  get_bulk().modification_begin();
  get_bulk().change_entity_parts<stk::mesh::PartVector>(partAElements, {&partA});
  get_bulk().change_entity_parts<stk::mesh::PartVector>(partBElements, {&partB});
  get_bulk().modification_end();

  const size_t numBuckets = get_bulk().buckets(stk::topology::ELEM_RANK).size();
  const size_t bucketCapacity = stk::mesh::get_default_maximum_bucket_capacity();
  const size_t numPerEntityA = 1;
  const size_t numPerEntityB = 5;
  const size_t maxNumPerEntity = std::max(numPerEntityA, numPerEntityB);
  const size_t bytesPerScalar = sizeof(double);
  const size_t expectedTotalFieldDataBytes = (numBuckets * bucketCapacity * maxNumPerEntity * bytesPerScalar);
  EXPECT_LE(stk::mesh::get_total_ngp_field_allocation_bytes(field), expectedTotalFieldDataBytes);
}

TEST_F(FieldFixture, writingDifferentNodalFieldsPerSolutionCase)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)<3)
    test_solution_case_with_rank(stk::topology::NODE_RANK);
}

TEST_F(FieldFixture, DISABLED_writingDifferentElementFieldsPerSolutionCase)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD)<3)
    test_solution_case_with_rank(stk::topology::ELEM_RANK);
}

TEST_F(FieldFixture, fenceWithoutNgpField)
{
  setup_empty_mesh(m_auraOption);
  stk::mesh::Field<double> &field = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "doubleField");

  EXPECT_NO_THROW(field.fence());
}

class LateFieldFixtureNoTest : public stk::unit_test_util::MeshFixtureNoTest
{
protected:
  LateFieldFixtureNoTest() {}
  LateFieldFixtureNoTest(unsigned spatial_dim) : MeshFixtureNoTest(spatial_dim) {}
  LateFieldFixtureNoTest(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
    : MeshFixtureNoTest(spatial_dim,entityRankNames)  {}

  virtual ~LateFieldFixtureNoTest()
  {
    reset_mesh();
  }
};


class LateFieldFixture : public LateFieldFixtureNoTest, public ::ngp_testing::Test
{
protected:
  void custom_allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                            std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager)
  {
    stk::mesh::MeshBuilder builder(communicator);
    builder.set_spatial_dimension(m_spatialDim);
    builder.set_entity_rank_names(m_entityRankNames);
    builder.set_aura_option(auraOption);
    builder.set_field_data_manager(std::move(fieldDataManager));

    bulkData = builder.create();
    metaData = bulkData->mesh_meta_data_ptr();
  }

  void setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                                std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager)
  {
    custom_allocate_bulk(auraOption, std::move(fieldDataManager));
    metaData->enable_late_fields();
  }

  virtual void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                                std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager)
  {
    custom_allocate_bulk(auraOption, std::move(fieldDataManager));
    metaData->enable_late_fields();
  }

  void create_part(const std::string & partName, stk::mesh::EntityRank rank, bool addIoPartAttribute = true)
  {
    stk::mesh::Part & part = get_meta().declare_part(partName, rank);
    if (addIoPartAttribute) {
      stk::io::put_io_part_attribute(part);
    }
  }

  template <typename T>
  stk::mesh::Field<T> & declare_field(const std::string & fieldName, stk::mesh::EntityRank rank)
  {
    const int numStates = 1;
    stk::mesh::Field<T> & field = get_meta().declare_field<T>(rank, fieldName, numStates);
    return field;
  }

  template <typename T>
  void put_field(stk::mesh::Field<T> & field, const stk::mesh::Part & part)
  {
    const T initVal = 123;
    stk::mesh::put_field_on_mesh(field, part, &initVal);
  }

  template <typename T>
  void set_field_values_with_scale_factor(stk::mesh::Field<T> & field, T scaleFactor)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(field.entity_rank(), field);
    for (stk::mesh::Bucket * bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        T * data = stk::mesh::field_data(field, node);
        data[0] = get_bulk().identifier(node) * scaleFactor;
      }
    }
  }

  template <typename T>
  void expect_field_values_with_scale_factor(stk::mesh::Field<T> & field, T scaleFactor)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(field.entity_rank(), field);
    for (stk::mesh::Bucket * bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const T * data = stk::mesh::field_data(field, node);
        EXPECT_EQ(static_cast<T>(get_bulk().identifier(node) * scaleFactor), data[0]) << "For field: " << field.name();
      }
    }
  }

  template <typename T>
  void setup_add_late_first_field(stk::mesh::EntityRank rank,
                                  std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                      std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    // Note that we still have a nodal coordinates field

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(lateField, 1);
    expect_field_values_with_scale_factor(lateField, 1);
  }

  template <typename T>
  void setup_add_late_field(stk::mesh::EntityRank rank,
                            std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_late_field_multiple_buckets(stk::mesh::EntityRank rank,
                                             std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                 std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);

    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    get_bulk().modification_begin();
    if (get_bulk().is_valid(elem2)) {
      stk::mesh::PartVector addParts(1, get_meta().get_part("block_2"));
      stk::mesh::PartVector removeParts(1, get_meta().get_part("block_1"));
      get_bulk().change_entity_parts(elem2, addParts, removeParts);
    }
    get_bulk().modification_end();

    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_late_field_multiple_duplicate_put_field(stk::mesh::EntityRank rank,
                                                         std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                             std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField, 2);

    put_field(lateField, get_meta().universal_part());

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_late_field_multiple_different_put_field(stk::mesh::EntityRank rank,
                                                         std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                             std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);

    put_field(lateField, *get_meta().get_part("block_1"));
    put_field(lateField, *get_meta().get_part("block_2"));
    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_two_late_fields_sequential(stk::mesh::EntityRank rank,
                                            std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField1 = declare_field<T>("late_field1", rank);
    put_field(lateField1, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField1, 2);

    stk::mesh::Field<T> & lateField2 = declare_field<T>("late_field2", rank);
    put_field(lateField2, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField2, 3);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField1, 2);
    expect_field_values_with_scale_factor(lateField2, 3);
  }

  template <typename T>
  void setup_add_two_late_fields_interleaved(stk::mesh::EntityRank rank,
                                             std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                 std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField1 = declare_field<T>("late_field1", rank);
    stk::mesh::Field<T> & lateField2 = declare_field<T>("late_field2", rank);

    put_field(lateField1, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField1, 2);

    put_field(lateField2, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField2, 3);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField1, 2);
    expect_field_values_with_scale_factor(lateField2, 3);
  }

  template <typename T>
  void setup_add_two_late_fields_out_of_order(stk::mesh::EntityRank rank,
                                              std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                  std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField1 = declare_field<T>("late_field1", rank);
    stk::mesh::Field<T> & lateField2 = declare_field<T>("late_field2", rank);

    put_field(lateField2, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField2, 3);

    put_field(lateField1, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField1, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField1, 2);
    expect_field_values_with_scale_factor(lateField2, 3);
  }

  template <typename T1, typename T2>
  void setup_add_two_late_fields_different_type_out_of_order(stk::mesh::EntityRank rank,
                                                             std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                                 std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T1> & earlyField = declare_field<T1>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor<T1>(earlyField, 1);

    stk::mesh::Field<T1> & lateField1    = declare_field<T1>("late_field1", rank);
    stk::mesh::Field<T2> & lateField2 = declare_field<T2>("late_field2", rank);

    put_field<T2>(lateField2, get_meta().universal_part());
    set_field_values_with_scale_factor<T2>(lateField2, 3);

    put_field<T1>(lateField1, get_meta().universal_part());
    set_field_values_with_scale_factor<T1>(lateField1, 2);

    expect_field_values_with_scale_factor<T1>(earlyField, 1);
    expect_field_values_with_scale_factor<T1>(lateField1, 2);
    expect_field_values_with_scale_factor<T2>(lateField2, 3);
  }

  template <typename T>
  void setup_add_two_late_fields_different_rank_out_of_order(stk::mesh::EntityRank rank1, stk::mesh::EntityRank rank2,
                                                             std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                                                 std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", stk::topology::ELEM_RANK);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField1 = declare_field<T>("late_field1", rank1);
    stk::mesh::Field<T> & lateField2 = declare_field<T>("late_field2", rank2);

    put_field(lateField2, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField2, 3);

    put_field(lateField1, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField1, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField1, 2);
    expect_field_values_with_scale_factor(lateField2, 3);
  }

  template <typename T>
  void setup_add_early_field_to_late_part(stk::mesh::EntityRank rank,
                                          std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                              std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    const bool addIoPartAttribute = false;
    create_part("block_1", stk::topology::ELEM_RANK, addIoPartAttribute);
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, *get_meta().get_part("block_1"));
    stk::io::fill_mesh("generated:1x1x2", *bulkData);

    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    get_bulk().modification_begin();
    if (get_bulk().is_valid(elem2)) {
      stk::mesh::PartVector addParts(1, get_meta().get_part("block_2"));
      stk::mesh::PartVector removeParts(1, get_meta().get_part("block_1"));
      get_bulk().change_entity_parts(elem2, addParts, removeParts);
    }
    get_bulk().modification_end();

    set_field_values_with_scale_factor(earlyField, 1);

    get_meta().declare_part("block_3", stk::topology::ELEM_RANK);
    get_bulk().modification_begin();
    if (get_bulk().is_valid(elem2)) {
      stk::mesh::PartVector addParts(1, get_meta().get_part("block_3"));
      stk::mesh::PartVector removeParts;
      get_bulk().change_entity_parts(elem2, addParts, removeParts);
    }
    get_bulk().modification_end();

    put_field(earlyField, *get_meta().get_part("block_3"));

    set_field_values_with_scale_factor(earlyField, 1);

    expect_field_values_with_scale_factor(earlyField, 1);
  }

  template <typename T>
  void setup_add_late_field_to_late_part(stk::mesh::EntityRank rank,
                                         std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager =
                                             std::unique_ptr<stk::mesh::FieldDataManager>())
  {
    setup_empty_mesh_with_field_data_manager(stk::mesh::BulkData::NO_AUTO_AURA, std::move(fieldDataManager));
    const bool addIoPartAttribute = false;
    create_part("block_1", stk::topology::ELEM_RANK, addIoPartAttribute);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    stk::io::fill_mesh("generated:1x1x2", *bulkData);
    set_field_values_with_scale_factor(earlyField, 1);

    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Entity elem2 = get_bulk().get_entity(stk::topology::ELEM_RANK, 2);
    get_bulk().modification_begin();
    if (get_bulk().is_valid(elem2)) {
      stk::mesh::PartVector addParts(1, get_meta().get_part("block_2"));
      stk::mesh::PartVector removeParts(1, get_meta().get_part("block_1"));
      get_bulk().change_entity_parts(elem2, addParts, removeParts);
    }
    get_bulk().modification_end();

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(earlyField, *get_meta().get_part("block_2"));

    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

};

TEST_F(LateFieldFixture, addLateIntFirstElementField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_first_field<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addLateIntNodalField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addLateIntElementField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, disable_late_fields)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field<int>(stk::topology::ELEM_RANK);

  get_meta().disable_late_fields();
  stk::mesh::Field<int>& f = declare_field<int>("another_late_field", stk::topology::ELEM_RANK);
  EXPECT_ANY_THROW(stk::mesh::put_field_on_mesh(f, get_meta().universal_part(), nullptr));
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleBuckets)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_multiple_buckets<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDuplicatePutField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDuplicatePutField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDifferentPutField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDifferentPutField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_sequential)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_sequential<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_sequential)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_sequential<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_interleaved)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_interleaved<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_interleaved)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_interleaved<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_out_of_order<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_out_of_order<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleNodalFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleElementFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleNodalFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleElementFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addTwoLateNodalAndElementFields_outOfOrder)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_two_late_fields_different_rank_out_of_order<int>(stk::topology::NODE_RANK, stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addEarlyNodalFieldToLatePart)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_early_field_to_late_part<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addEarlyElementFieldToLatePart)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_early_field_to_late_part<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addLateNodalFieldToLatePart)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_to_late_part<int>(stk::topology::NODE_RANK);
}

TEST_F(LateFieldFixture, addLateElementFieldToLatePart)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  setup_add_late_field_to_late_part<int>(stk::topology::ELEM_RANK);
}

TEST_F(LateFieldFixture, addLateIntFirstElementFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_first_field<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntNodalFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntElementFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleBucketsContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_multiple_buckets<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDuplicatePutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDuplicatePutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDifferentPutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDifferentPutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_sequentialContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_sequential<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_sequentialContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_sequential<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_interleavedContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_interleaved<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_interleavedContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_interleaved<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_out_of_order<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_out_of_order<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addTwoLateNodalAndElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_two_late_fields_different_rank_out_of_order<int>(stk::topology::NODE_RANK, stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addEarlyNodalFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_early_field_to_late_part<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addEarlyElementFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_early_field_to_late_part<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateNodalFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_to_late_part<int>(stk::topology::NODE_RANK, std::move(fieldDataManager));
}

TEST_F(LateFieldFixture, addLateElementFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  auto fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>();
  setup_add_late_field_to_late_part<int>(stk::topology::ELEM_RANK, std::move(fieldDataManager));
}

TEST(SharedSidesetField, verifySidesetFieldAfterMeshRead) {
  std::string serialOutputMeshName = "ARB.e";

  if (stk::parallel_machine_size(MPI_COMM_WORLD) == 2)
  {
    const std::string fieldName = "surface_1_df";
    const double initValue = 123.0;

    // Build the target serial mesh and write it to a file
    if (stk::parallel_machine_rank(MPI_COMM_WORLD) == 0)
    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_SELF);
      stk::mesh::BulkData& bulk = *bulkPtr;

      stk::unit_test_util::create_AB_mesh_with_sideset_and_distribution_factors(bulk,
                                                                                               stk::unit_test_util::LEFT,
                                                                                               stk::unit_test_util::DECREASING,
                                                                                               fieldName,
                                                                                               initValue);

      stk::io::write_mesh_with_fields(serialOutputMeshName, bulk, 1, 1.0);
    }

    {
      std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(3, MPI_COMM_WORLD);
      stk::mesh::BulkData& bulk = *bulkPtr;
      stk::mesh::MetaData& meta = bulkPtr->mesh_meta_data();
      stk::io::StkMeshIoBroker stkIo;
      stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));

      stkIo.set_bulk_data(bulk);
      stkIo.add_mesh_database(serialOutputMeshName, stk::io::READ_MESH);
      stkIo.create_input_mesh();
      stkIo.add_all_mesh_fields_as_input_fields();

      stkIo.populate_bulk_data();

      const stk::mesh::BucketVector& buckets = bulk.get_buckets(stk::topology::FACE_RANK, meta.universal_part());
      ASSERT_EQ(1u, buckets.size());
      for (stk::mesh::Bucket* bucket : buckets)
      {
        ASSERT_EQ(1u, bucket->size());
        stk::mesh::Entity face = (*bucket)[0];
        stk::mesh::FieldBase* field = meta.get_field(stk::topology::FACE_RANK, fieldName);
        EXPECT_NE(nullptr, field);
        unsigned numEntries = stk::mesh::field_scalars_per_entity(*field, face);
        double* fieldData = reinterpret_cast<double*>(stk::mesh::field_data(*field, face));
        for (unsigned entry = 0; entry < numEntries; ++entry)
        {
          EXPECT_NEAR(initValue, fieldData[entry], 1e-12);
        }
      }
    }
  }
  unlink(serialOutputMeshName.c_str());
}

void create_node(stk::mesh::BulkData & bulk, stk::mesh::EntityId nodeId, stk::mesh::Part & part)
{
  bulk.modification_begin();
  bulk.declare_node(nodeId, stk::mesh::PartVector{&part});
  bulk.modification_end();
}

void change_node_parts(stk::mesh::BulkData & bulk, stk::mesh::EntityId nodeId,
                       stk::mesh::Part & addPart, stk::mesh::Part & removePart)
{
  const stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
  bulk.modification_begin();
  bulk.change_entity_parts(node, stk::mesh::PartVector{&addPart}, stk::mesh::PartVector{&removePart});
  bulk.modification_end();
}

class VariableCapacityBuckets : public ::testing::Test
{
public:
  VariableCapacityBuckets()
  {
  }

  void build_empty_mesh(unsigned initialBucketCapacity, unsigned maximumBucketCapacity)
  {
    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
    stk::mesh::get_updated_ngp_mesh(*m_bulk);
  }

protected:
  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData * m_meta;
};

void check_num_buckets(const stk::mesh::BulkData & bulk, unsigned expectedNumBuckets)
{
  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(buckets.size(), expectedNumBuckets);

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  ASSERT_EQ(ngpMesh.num_buckets(stk::topology::NODE_RANK), expectedNumBuckets);
}

void check_bucket_sizes(const stk::mesh::BulkData & bulk, const std::vector<unsigned> & expectedBucketSizes)
{
  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(buckets.size(), expectedBucketSizes.size());
  for (unsigned i = 0; i < buckets.size(); ++i) {
    EXPECT_EQ(buckets[i]->size(), expectedBucketSizes[i]);
  }

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  ASSERT_EQ(ngpMesh.num_buckets(stk::topology::NODE_RANK), expectedBucketSizes.size());
  for (unsigned i = 0; i < ngpMesh.num_buckets(stk::topology::NODE_RANK); ++i) {
    EXPECT_EQ(ngpMesh.get_bucket(stk::topology::NODE_RANK, i).size(), expectedBucketSizes[i]);
  }
}

void check_bucket_capacities(const stk::mesh::BulkData & bulk, const std::vector<unsigned> & expectedBucketCapacities)
{
  const stk::mesh::BucketVector & buckets = bulk.buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(buckets.size(), expectedBucketCapacities.size());
  for (unsigned i = 0; i < buckets.size(); ++i) {
    EXPECT_EQ(buckets[i]->capacity(), expectedBucketCapacities[i]);
  }

  stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
  ASSERT_EQ(ngpMesh.num_buckets(stk::topology::NODE_RANK), expectedBucketCapacities.size());
  for (unsigned i = 0; i < ngpMesh.num_buckets(stk::topology::NODE_RANK); ++i) {
    EXPECT_EQ(ngpMesh.get_bucket(stk::topology::NODE_RANK, i).capacity(), expectedBucketCapacities[i]);
  }
}

constexpr stk::topology::rank_t bucketRank = stk::topology::NODE_RANK;

TEST_F(VariableCapacityBuckets, createNodes_initialCapacity1_maxCapacity1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 1);

  stk::mesh::Part & block_1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  {
    SCOPED_TRACE("Create Node 1");
    create_node(*m_bulk, 1, block_1);

    check_num_buckets(*m_bulk, 1);
    check_bucket_sizes(*m_bulk, {1});
    check_bucket_capacities(*m_bulk, {1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 2");
    create_node(*m_bulk, 2, block_1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {1, 1});
    check_bucket_capacities(*m_bulk, {1, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1}}, {"block_1", {2}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 3");
    create_node(*m_bulk, 3, block_1);

    check_num_buckets(*m_bulk, 3);
    check_bucket_sizes(*m_bulk, {1, 1, 1});
    check_bucket_capacities(*m_bulk, {1, 1, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1}}, {"block_1", {2}}, {"block_1", {3}}}, bucketRank);
  }
}

TEST_F(VariableCapacityBuckets, createNodes_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  stk::mesh::Part & block_1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  {
    SCOPED_TRACE("Create Node 1");
    create_node(*m_bulk, 1, block_1);

    check_num_buckets(*m_bulk, 1);
    check_bucket_sizes(*m_bulk, {1});
    check_bucket_capacities(*m_bulk, {2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 2");
    create_node(*m_bulk, 2, block_1);

    check_num_buckets(*m_bulk, 1);
    check_bucket_sizes(*m_bulk, {2});
    check_bucket_capacities(*m_bulk, {2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 3");
    create_node(*m_bulk, 3, block_1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
  }
}

TEST_F(VariableCapacityBuckets, createNodes_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  stk::mesh::Part & block_1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  {
    SCOPED_TRACE("Create Node 1");
    create_node(*m_bulk, 1, block_1);

    check_num_buckets(*m_bulk, 1);
    check_bucket_sizes(*m_bulk, {1});
    check_bucket_capacities(*m_bulk, {1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 2");
    create_node(*m_bulk, 2, block_1);

    check_num_buckets(*m_bulk, 1);
    check_bucket_sizes(*m_bulk, {2});
    check_bucket_capacities(*m_bulk, {2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Create Node 3");
    create_node(*m_bulk, 3, block_1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
  }
}

TEST_F(VariableCapacityBuckets, changeNodeParts_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);
  stk::mesh::Part & block2 = m_meta->declare_part_with_topology("block_2", stk::topology::NODE);

  {
    SCOPED_TRACE("Create Nodes 1, 2, 3");
    create_node(*m_bulk, 1, block1);
    create_node(*m_bulk, 2, block1);
    create_node(*m_bulk, 3, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Change parts for Node 1");
    change_node_parts(*m_bulk, 1, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {2, 3}}, {"block_2", {1}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Change parts for Node 2");
    change_node_parts(*m_bulk, 2, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {1, 2});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {3}}, {"block_2", {1, 2}}}, bucketRank);
  }
  {
    SCOPED_TRACE("Change parts for Node 3");
    change_node_parts(*m_bulk, 3, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_2", {1, 2}}, {"block_2", {3}}}, bucketRank);
  }
}

TEST_F(VariableCapacityBuckets, changeNodeParts_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);
  stk::mesh::Part & block2 = m_meta->declare_part_with_topology("block_2", stk::topology::NODE);

  {
    SCOPED_TRACE("Create Nodes 1, 2, 3");
    create_node(*m_bulk, 1, block1);
    create_node(*m_bulk, 2, block1);
    create_node(*m_bulk, 3, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
  }

  {
    SCOPED_TRACE("Change parts for Node 1");
    change_node_parts(*m_bulk, 1, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {2, 3}}, {"block_2", {1}}}, bucketRank);
  }

  {
    SCOPED_TRACE("Change parts for Node 2");
    change_node_parts(*m_bulk, 2, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {1, 2});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {3}}, {"block_2", {1, 2}}}, bucketRank);
  }

  {
    SCOPED_TRACE("Change parts for Node 3");
    change_node_parts(*m_bulk, 3, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_2", {1, 2}}, {"block_2", {3}}}, bucketRank);
  }
}

TEST_F(VariableCapacityBuckets, initialMeshConstruction_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  {
    stk::mesh::EntityIdVector ids{1, 2, 3};
    stk::mesh::EntityVector newNodes;
    m_bulk->modification_begin();
    m_bulk->declare_entities(stk::topology::NODE_RANK, ids, stk::mesh::PartVector{&block1}, newNodes);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
    m_bulk->modification_end();
  }
}

TEST_F(VariableCapacityBuckets, initialMeshConstruction_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  {
    stk::mesh::EntityIdVector ids{1, 2, 3};
    stk::mesh::EntityVector newNodes;
    m_bulk->modification_begin();
    m_bulk->declare_entities(stk::topology::NODE_RANK, ids, stk::mesh::PartVector{&block1}, newNodes);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
    m_bulk->modification_end();
  }
}


void create_node_with_data(stk::mesh::BulkData & bulk, stk::mesh::EntityId nodeId,
                           const stk::mesh::Field<int> & field)
{
  bulk.modification_begin();
  const stk::mesh::Entity node = bulk.declare_node(nodeId);
  bulk.modification_end();

  *stk::mesh::field_data(field, node) = nodeId;
}

void create_node_with_multistate_data(stk::mesh::BulkData & bulk, stk::mesh::EntityId nodeId,
                                      const stk::mesh::Field<int> & field1, const stk::mesh::Field<int> & field2)
{
  bulk.modification_begin();
  const stk::mesh::Entity node = bulk.declare_node(nodeId);
  bulk.modification_end();

  *stk::mesh::field_data(field1, node) = 100 + nodeId;
  *stk::mesh::field_data(field2, node) = 200 + nodeId;
}

void create_node_with_data(stk::mesh::BulkData & bulk, stk::mesh::EntityId nodeId,
                           const stk::mesh::Field<int> & field, stk::mesh::Part & part)
{
  bulk.modification_begin();
  const stk::mesh::Entity node = bulk.declare_node(nodeId, stk::mesh::PartVector{&part});
  bulk.modification_end();

  *stk::mesh::field_data(field, node) = nodeId;
}

enum class FieldDataManagerType
{
  DEFAULT,
  CONTIGUOUS
};

using FieldDataType = int;

class VariableCapacityFieldData : public ::testing::TestWithParam<FieldDataManagerType>
{
public:
  VariableCapacityFieldData()
    : m_fieldDataManager(nullptr)
  {
  }

  void build_empty_mesh(unsigned initialBucketCapacity, unsigned maximumBucketCapacity)
  {
    std::unique_ptr<stk::mesh::FieldDataManager> fieldDataManager;
    if (GetParam() == FieldDataManagerType::DEFAULT) {
      const unsigned numRanks = 1;
      const unsigned alignment = 4;
      fieldDataManager = std::make_unique<stk::mesh::DefaultFieldDataManager>(numRanks, alignment);
    }
    else if (GetParam() == FieldDataManagerType::CONTIGUOUS) {
      const unsigned extraCapacity = 0;
      const unsigned alignment = 4;
      fieldDataManager = std::make_unique<stk::mesh::ContiguousFieldDataManager>(extraCapacity, alignment);
    }

    m_fieldDataManager = fieldDataManager.get();

    stk::mesh::MeshBuilder builder(MPI_COMM_WORLD);
    builder.set_spatial_dimension(3);
    builder.set_initial_bucket_capacity(initialBucketCapacity);
    builder.set_maximum_bucket_capacity(maximumBucketCapacity);
    builder.set_field_data_manager(std::move(fieldDataManager));
    m_bulk = builder.create();
    m_meta = &m_bulk->mesh_meta_data();
  }

  int expected_bytes_allocated_host(const stk::mesh::BucketVector & buckets, int dataSize)
  {
    int extraCapacity = 0;
    if (GetParam() == FieldDataManagerType::CONTIGUOUS) {
      extraCapacity = dynamic_cast<stk::mesh::ContiguousFieldDataManager&>(*m_fieldDataManager).get_extra_capacity();
    }

    return std::accumulate(buckets.begin(), buckets.end(), 0,
                           [&](int currentValue, const stk::mesh::Bucket * bucket) {
                              const int allocationCount = (GetParam() == FieldDataManagerType::DEFAULT) ? bucket->capacity()
                                                                                                        : bucket->size();
                              return currentValue + stk::adjust_up_to_alignment_boundary(dataSize*allocationCount,
                                                                                         m_fieldDataManager->get_alignment_bytes());
                           }) + extraCapacity;
  }

  int expected_bytes_allocated_device(const stk::mesh::BulkData & bulk,
                                      const stk::mesh::FieldBase & stkField,
                                      const stk::mesh::BucketVector & buckets, int dataSize)
  {
    return buckets.size() * bulk.get_maximum_bucket_capacity() * stkField.max_size() * dataSize;
  }

  void check_expected_bytes_allocated(const stk::mesh::BulkData & bulk,
                                      const stk::mesh::FieldBase & stkField)
  {
    const int dataSize = sizeof(FieldDataType);
    const unsigned fieldOrdinal = stkField.mesh_meta_data_ordinal();
    const stk::mesh::BucketVector & buckets = m_bulk->buckets(stk::topology::NODE_RANK);
    const int bytesAllocated = m_fieldDataManager->get_num_bytes_allocated_on_field(fieldOrdinal);
    ASSERT_EQ(bytesAllocated, expected_bytes_allocated_host(buckets, dataSize));

#ifdef STK_USE_DEVICE_MESH
    if (std::is_same_v<stk::mesh::NgpField<FieldDataType>, stk::mesh::DeviceField<FieldDataType>>) {
      stk::mesh::NgpField<FieldDataType> & ngpField = stk::mesh::get_updated_ngp_field<FieldDataType>(stkField);
      const auto fieldData = stk::mesh::impl::get_device_data(ngpField);
      const int bytesAllocatedDevice = fieldData.extent(0) * fieldData.extent(1) * fieldData.extent(2) *
                                       sizeof(FieldDataType);
      ASSERT_EQ(bytesAllocatedDevice, expected_bytes_allocated_device(bulk, stkField, buckets, dataSize));
    }
#endif
  }

  void check_field_values(const stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & stkField)
  {
    const stk::mesh::BucketVector & buckets = m_bulk->buckets(stk::topology::NODE_RANK);
    unsigned nodeIdx = 0;
    for (const stk::mesh::Bucket * bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const FieldDataType fieldValue = *static_cast<FieldDataType*>(stk::mesh::field_data(stkField, node));
        const FieldDataType expectedValue = bulk.identifier(node);
        EXPECT_EQ(fieldValue, expectedValue);
      }
    }

    const unsigned numNodes = stk::mesh::count_entities(bulk, stk::topology::NODE_RANK,
                                                        bulk.mesh_meta_data().universal_part());
    Kokkos::View<FieldDataType*> deviceValues("deviceValues", numNodes);
    Kokkos::View<FieldDataType*>::HostMirror hostValuesFromDevice = Kokkos::create_mirror_view(deviceValues);

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    stk::mesh::NgpField<FieldDataType> ngpField = stk::mesh::get_updated_ngp_field<FieldDataType>(stkField);
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(size_t /*index*/) {
        unsigned nodeGlobalIndex = 0;
        const unsigned numBuckets = ngpMesh.num_buckets(stk::topology::NODE_RANK);
        for (unsigned bucketId = 0; bucketId < numBuckets; ++bucketId) {
          const stk::mesh::NgpMesh::BucketType & deviceBucket = ngpMesh.get_bucket(stk::topology::NODE_RANK, bucketId);
          const unsigned numNodesInBucket = deviceBucket.size();
          for (unsigned nodeOrdinal = 0; nodeOrdinal < numNodesInBucket; ++nodeOrdinal) {
            const stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(deviceBucket[nodeOrdinal]);
            deviceValues[nodeGlobalIndex++] = ngpField.get(nodeIndex, 0);
          }
        }
      });

    Kokkos::deep_copy(hostValuesFromDevice, deviceValues);

    nodeIdx = 0;
    for (const stk::mesh::Bucket * bucket : buckets) {
      for (stk::mesh::Entity node : *bucket) {
        const FieldDataType expectedValue = bulk.identifier(node);
        EXPECT_EQ(hostValuesFromDevice[nodeIdx++], expectedValue);
      }
    }
  }

  void check_field_value(const stk::mesh::BulkData & bulk, const stk::mesh::FieldBase & stkField,
                         const stk::mesh::EntityId nodeId, const int expectedValue)
  {
    stk::mesh::Entity node = bulk.get_entity(stk::topology::NODE_RANK, nodeId);
    const FieldDataType fieldValue = *static_cast<FieldDataType*>(stk::mesh::field_data(stkField, node));
    EXPECT_EQ(fieldValue, expectedValue);

    Kokkos::View<FieldDataType*> deviceValue("deviceValues", 1);
    Kokkos::View<FieldDataType*>::HostMirror hostValueFromDevice = Kokkos::create_mirror_view(deviceValue);

    stk::mesh::NgpMesh & ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);
    stk::mesh::NgpField<FieldDataType> ngpField = stk::mesh::get_updated_ngp_field<FieldDataType>(stkField);
    Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(size_t /*index*/) {
        deviceValue[0] = ngpField.get(ngpMesh, node, 0);
      });

    Kokkos::deep_copy(hostValueFromDevice, deviceValue);

    EXPECT_EQ(hostValueFromDevice[0], expectedValue);
  }

protected:
  stk::mesh::FieldDataManager * m_fieldDataManager;
  std::unique_ptr<stk::mesh::BulkData> m_bulk;
  stk::mesh::MetaData * m_meta;
};

INSTANTIATE_TEST_SUITE_P(VariableCapacityBuckets, VariableCapacityFieldData,
                         ::testing::Values(FieldDataManagerType::DEFAULT, FieldDataManagerType::CONTIGUOUS));

TEST_P(VariableCapacityFieldData, createNodes_initialCapacity1_maxCapacity1)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 1);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Node 1");
    create_node_with_data(*m_bulk, 1, field);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 2");
    create_node_with_data(*m_bulk, 2, field);

    check_num_buckets(*m_bulk, 2);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 3");
    create_node_with_data(*m_bulk, 3, field);

    check_num_buckets(*m_bulk, 3);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
}

TEST_P(VariableCapacityFieldData, createNodes_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Node 1");
    create_node_with_data(*m_bulk, 1, field);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 2");
    create_node_with_data(*m_bulk, 2, field);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 3");
    create_node_with_data(*m_bulk, 3, field);

    check_num_buckets(*m_bulk, 2);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
}

TEST_P(VariableCapacityFieldData, createNodes_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Node 1");
    create_node_with_data(*m_bulk, 1, field);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 2");
    create_node_with_data(*m_bulk, 2, field);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Create Node 3");
    create_node_with_data(*m_bulk, 3, field);

    check_num_buckets(*m_bulk, 2);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
}

TEST_P(VariableCapacityFieldData, createNodes_initialCapacity2_maxCapacity2_withMultistateFieldData)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  const unsigned numberOfStates = 2;
  auto & field1 = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1", numberOfStates);
  auto & field2 = dynamic_cast<stk::mesh::Field<int>&>(*field1.field_state(stk::mesh::FieldState::StateOld));
  stk::mesh::put_field_on_mesh(field1, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Node 1");
    create_node_with_multistate_data(*m_bulk, 1, field1, field2);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field1);
    check_expected_bytes_allocated(*m_bulk, field2);
    check_field_value(*m_bulk, field1, 1, 101);
    check_field_value(*m_bulk, field2, 1, 201);
  }

  m_bulk->update_field_data_states();  // Default to not rotate on device; will rotate anyway during mesh mod below

  {
    SCOPED_TRACE("Create Node 2");
    create_node_with_multistate_data(*m_bulk, 2, field1, field2);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field1);
    check_expected_bytes_allocated(*m_bulk, field2);
    check_field_value(*m_bulk, field1, 1, 201);  // Flipped values due to state rotation
    check_field_value(*m_bulk, field2, 1, 101);
    check_field_value(*m_bulk, field1, 2, 102);  // Unflipped values written into new state layout
    check_field_value(*m_bulk, field2, 2, 202);
  }
}

TEST_P(VariableCapacityFieldData, createNodes_initialCapacity1_maxCapacity2_withMultistateFieldData)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  const unsigned numberOfStates = 2;
  auto & field1 = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1", numberOfStates);
  auto & field2 = dynamic_cast<stk::mesh::Field<int>&>(*field1.field_state(stk::mesh::FieldState::StateOld));
  stk::mesh::put_field_on_mesh(field1, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Node 1");
    create_node_with_multistate_data(*m_bulk, 1, field1, field2);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field1);
    check_expected_bytes_allocated(*m_bulk, field2);
    check_field_value(*m_bulk, field1, 1, 101);
    check_field_value(*m_bulk, field2, 1, 201);
    create_node_with_multistate_data(*m_bulk, 1, field1, field2);
  }

  m_bulk->update_field_data_states();  // Default to not rotate on device; will rotate anyway during mesh mod below

  {
    SCOPED_TRACE("Create Node 2");
    create_node_with_multistate_data(*m_bulk, 2, field1, field2);

    check_num_buckets(*m_bulk, 1);
    check_expected_bytes_allocated(*m_bulk, field1);
    check_expected_bytes_allocated(*m_bulk, field2);
    check_field_value(*m_bulk, field1, 1, 201);  // Flipped values due to state rotation
    check_field_value(*m_bulk, field2, 1, 101);
    check_field_value(*m_bulk, field1, 2, 102);  // Unflipped values written into new state layout
    check_field_value(*m_bulk, field2, 2, 202);
  }
}

TEST_P(VariableCapacityFieldData, changeNodeParts_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);
  stk::mesh::Part & block2 = m_meta->declare_part_with_topology("block_2", stk::topology::NODE);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Nodes 1, 2, 3");
    create_node_with_data(*m_bulk, 1, field, block1);
    create_node_with_data(*m_bulk, 2, field, block1);
    create_node_with_data(*m_bulk, 3, field, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Move Node 1 from block_1 to block_2");
    change_node_parts(*m_bulk, 1, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {2, 3}}, {"block_2", {1}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Move Node 2 from block_1 to block_2");
    change_node_parts(*m_bulk, 2, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {1, 2});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {3}}, {"block_2", {1, 2}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
  {
    SCOPED_TRACE("Move Node 3 from block_1 to block_2");
    change_node_parts(*m_bulk, 3, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_2", {1, 2}}, {"block_2", {3}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
}

TEST_P(VariableCapacityFieldData, changeNodeParts_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);
  stk::mesh::Part & block2 = m_meta->declare_part_with_topology("block_2", stk::topology::NODE);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, m_meta->universal_part(), nullptr);

  {
    SCOPED_TRACE("Create Nodes 1, 2, 3");
    create_node_with_data(*m_bulk, 1, field, block1);
    create_node_with_data(*m_bulk, 2, field, block1);
    create_node_with_data(*m_bulk, 3, field, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {1, 2}}, {"block_1", {3}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Move Node 1 from block_1 to block_2");
    change_node_parts(*m_bulk, 1, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {2, 3}}, {"block_2", {1}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }

  {
    SCOPED_TRACE("Move Node 2 from block_1 to block_2");
    change_node_parts(*m_bulk, 2, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {1, 2});
    check_bucket_capacities(*m_bulk, {2, 2});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_1", {3}}, {"block_2", {1, 2}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
  {
    SCOPED_TRACE("Move Node 3 from block_1 to block_2");
    change_node_parts(*m_bulk, 3, block2, block1);

    check_num_buckets(*m_bulk, 2);
    check_bucket_sizes(*m_bulk, {2, 1});
    check_bucket_capacities(*m_bulk, {2, 1});
    ngp_unit_test_utils::check_bucket_layout(*m_bulk, {{"block_2", {1, 2}}, {"block_2", {3}}}, bucketRank);
    check_expected_bytes_allocated(*m_bulk, field);
    check_field_values(*m_bulk, field);
  }
}

TEST_P(VariableCapacityFieldData, initialMeshConstruction_initialCapacity2_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(2, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, block1, nullptr);

  {
    m_bulk->deactivate_field_updating();
    stk::mesh::EntityIdVector ids{1, 2, 3};
    stk::mesh::EntityVector newNodes;

    m_bulk->modification_begin();
    m_bulk->declare_entities(stk::topology::NODE_RANK, ids, stk::mesh::PartVector{&block1}, newNodes);
    m_bulk->allocate_field_data();

    check_num_buckets(*m_bulk, 2);
    check_expected_bytes_allocated(*m_bulk, field);
    m_bulk->modification_end();
  }
}

TEST_P(VariableCapacityFieldData, initialMeshConstruction_initialCapacity1_maxCapacity2)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) return;

  build_empty_mesh(1, 2);

  stk::mesh::Part & block1 = m_meta->declare_part_with_topology("block_1", stk::topology::NODE);

  auto & field = m_meta->declare_field<int>(stk::topology::NODE_RANK, "field1");
  stk::mesh::put_field_on_mesh(field, block1, nullptr);

  {
    m_bulk->deactivate_field_updating();
    stk::mesh::EntityIdVector ids{1, 2, 3};
    stk::mesh::EntityVector newNodes;

    m_bulk->modification_begin();
    m_bulk->declare_entities(stk::topology::NODE_RANK, ids, stk::mesh::PartVector{&block1}, newNodes);
    m_bulk->allocate_field_data();

    check_num_buckets(*m_bulk, 2);
    check_expected_bytes_allocated(*m_bulk, field);
    m_bulk->modification_end();
  }
}

} //namespace <anonymous>

