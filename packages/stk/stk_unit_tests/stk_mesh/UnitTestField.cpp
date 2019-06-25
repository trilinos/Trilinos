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
#include <stddef.h>                     // for size_t, NULL
#include <unistd.h>                     // for unlink
#include <iostream>                     // for ostream, operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <string>                       // for string, operator==, etc
#include <vector>                       // for vector
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include "stk_util/environment/WallTime.hpp"
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL
#include "Ioss_ElementBlock.h"          // for ElementBlock
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Property.h"              // for Property
#include "Ioss_Region.h"                // for Region, etc
#include "Shards_Array.hpp"
#include "mpi.h"                        // for ompi_communicator_t, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/MeshField.hpp"         // for MeshField
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_bytes_per_entity, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for operator<<, Selector, etc
#include "stk_mesh/base/Types.hpp"      // for BucketVector, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace Ioss { class DatabaseIO; }
namespace stk { namespace mesh { class Bucket; } }

namespace {

const stk::topology::rank_t NODE_RANK = stk::topology::NODE_RANK;

typedef shards::ArrayDimTag::size_type size_type;

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( CTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( DTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( ETAG )


TEST(UnitTestField, testFieldMaxSize)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for some test fields
  typedef stk::mesh::Field<double>                          rank_zero_field;
  typedef stk::mesh::Field<double,ATAG>                     rank_one_field;
  typedef stk::mesh::Field<double,ATAG,BTAG>                rank_two_field;
  typedef stk::mesh::Field<double,ATAG,BTAG,CTAG>           rank_three_field;
  typedef stk::mesh::Field<double,ATAG,BTAG,CTAG,DTAG>      rank_four_field;
  typedef stk::mesh::Field<double,ATAG,BTAG,CTAG,DTAG,ETAG> rank_five_field;

  const std::string name0("test_field_0");
  const std::string name1("test_field_1");
  const std::string name2("test_field_2");
  const std::string name3("test_field_3");
  const std::string name4("test_field_4");
  const std::string name5("test_field_5");

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field  >( NODE_RANK, name0 );
  rank_one_field   & f1 = meta_data.declare_field< rank_one_field   >( NODE_RANK, name1 );
  rank_two_field   & f2 = meta_data.declare_field< rank_two_field   >( NODE_RANK, name2 );
  rank_three_field & f3 = meta_data.declare_field< rank_three_field >( NODE_RANK, name3 );
  rank_four_field  & f4 = meta_data.declare_field< rank_four_field  >( NODE_RANK, name4 );
  rank_five_field  & f5 = meta_data.declare_field< rank_five_field  >( NODE_RANK, name5 );

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK );
  stk::mesh::Part & p2 = meta_data.declare_part("P2", NODE_RANK );
  stk::mesh::Part & p3 = meta_data.declare_part("P3", NODE_RANK );
  stk::mesh::Part & p4 = meta_data.declare_part("P4", NODE_RANK );
  stk::mesh::Part & p5 = meta_data.declare_part("P5", NODE_RANK );

  stk::mesh::put_field_on_mesh( f0, p0, (stk::mesh::FieldTraits<rank_zero_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f1, p1, 10, (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f2, p2, 10, 20, (stk::mesh::FieldTraits<rank_two_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f3, p3, 10, 20, 30, (stk::mesh::FieldTraits<rank_three_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f4, p4, 10, 20, 30, 40, (stk::mesh::FieldTraits<rank_four_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f5, p5, 10, 20, 30, 40, 50, (stk::mesh::FieldTraits<rank_five_field>::data_type*) nullptr);

  meta_data.commit();

  EXPECT_EQ( f0.max_size(stk::topology::NODE_RANK), 1u );
  EXPECT_EQ( f0.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f0.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f0.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f1.max_size(stk::topology::NODE_RANK), 10u );
  EXPECT_EQ( f1.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f1.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f1.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f2.max_size(stk::topology::NODE_RANK), 200u );
  EXPECT_EQ( f2.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f2.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f2.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f3.max_size(stk::topology::NODE_RANK), 6000u );
  EXPECT_EQ( f3.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f3.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f3.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f4.max_size(stk::topology::NODE_RANK), 240000u );
  EXPECT_EQ( f4.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f4.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f4.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f5.max_size(stk::topology::NODE_RANK), 12000000u );
  EXPECT_EQ( f5.max_size(stk::topology::EDGE_RANK), 0u );
  EXPECT_EQ( f5.max_size(stk::topology::FACE_RANK), 0u );
  EXPECT_EQ( f5.max_size(stk::topology::ELEMENT_RANK), 0u );

  EXPECT_EQ( f0.field_array_rank(), 0u ); // Field Rank NOT entity rank
  EXPECT_EQ( f1.field_array_rank(), 1u );
  EXPECT_EQ( f2.field_array_rank(), 2u );
  EXPECT_EQ( f3.field_array_rank(), 3u );
  EXPECT_EQ( f4.field_array_rank(), 4u );
  EXPECT_EQ( f5.field_array_rank(), 5u );

}

TEST(UnitTestField, testFieldWithSelector)
{
  stk::ParallelMachine pm = MPI_COMM_SELF ;
  std::ostringstream oss; // to test printing of things w/out spamming cout

  // specifications for test field
  typedef stk::mesh::Field<double>    rank_zero_field ;

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  rank_zero_field  & f0 = meta_data.declare_field< rank_zero_field >( NODE_RANK, name0 );

  stk::mesh::Part & p0 = meta_data.declare_part("P0", NODE_RANK );
  stk::mesh::Part & p1 = meta_data.declare_part("P1", NODE_RANK );

  stk::mesh::Selector select_p0 = p0;
  std::cout <<"select_p0: "<< select_p0 << std::endl;

  stk::mesh::put_field_on_mesh( f0 , select_p0 , (stk::mesh::FieldTraits<rank_zero_field>::data_type*) nullptr);

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

  typedef stk::mesh::Field<double,shards::ArrayDimension>           rank_one_field ;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  rank_one_field  & f0 = meta_data.declare_field< rank_one_field >( stk::topology::ELEMENT_RANK, name0 );

  stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  stk::mesh::Part & elements = meta_data.declare_part("Elements", elem_rank);
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );
  stk::mesh::Part & tet4s = meta_data.declare_part("Tet4", elem_rank );

  stk::mesh::Selector elem_hex_selector = elements & hex8s;
  stk::mesh::Selector elem_tet_selector = elements & tet4s;
  std::cout <<"elem_hex_selector: "<< elem_hex_selector << std::endl;
  std::cout <<"elem_tet_selector: "<< elem_tet_selector << std::endl;

  stk::mesh::put_field_on_mesh( f0 , elem_hex_selector, 8u , (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr);
  stk::mesh::put_field_on_mesh( f0 , elem_tet_selector, 4u , (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr);

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

  typedef stk::mesh::Field<double,shards::ArrayDimension>           rank_one_field ;
  // specifications for test field

  const std::string name0("test_field_0");

  const int spatial_dimension = 3;
  stk::mesh::MetaData meta_data( spatial_dimension );
  stk::mesh::BulkData bulk_data( meta_data , pm );

  rank_one_field  & f0 = meta_data.declare_field< rank_one_field >( stk::topology::ELEMENT_RANK, name0 );

  stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  stk::mesh::Part & hex8s = meta_data.declare_part("Hex8", elem_rank );

  stk::mesh::Part & universal_part = meta_data.universal_part();
  stk::mesh::Selector elem_hexA_selector = hex8s;
  stk::mesh::Selector elem_hexB_selector = universal_part & hex8s;

  std::cout <<"elem_hexA_selector: "<< elem_hexA_selector << std::endl;
  std::cout <<"elem_hexB_selector: "<< elem_hexB_selector << std::endl;

  stk::mesh::put_field_on_mesh( f0 , elem_hexA_selector, 8u , (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr);
  ASSERT_THROW(
    stk::mesh::put_field_on_mesh( f0 , elem_hexA_selector, 4u , (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr),
    std::runtime_error
  );
  stk::mesh::put_field_on_mesh( f0 , elem_hexB_selector, 4u , (stk::mesh::FieldTraits<rank_one_field>::data_type*) nullptr);

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

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ATAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( BTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( CTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( DTAG )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( ETAG )


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

        stk::mesh::Field<double> &nodeField = stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, fieldName, 1);
        stk::mesh::put_field_on_mesh(nodeField, stkIo.meta_data().universal_part(), &nodeInitialValue);

        stk::mesh::Field<double> &elemField = stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName, 1);
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

        stk::mesh::Field<double> &nodeField = stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK, fieldName, 1);
        stk::mesh::put_field_on_mesh(nodeField, stkIo.meta_data().universal_part(), &badInitialData);

        stk::mesh::Field<double> &elemField = stkIo.meta_data().declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK, fieldName, 1);
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
        stk::mesh::Field<double> *myTemplatedField = stk::mesh::get_field_by_name<stk::mesh::Field<double> >(fieldName, metaData);
        ASSERT_TRUE(myTemplatedField != NULL);
        EXPECT_TRUE( &nodeField == myTemplatedField);
        stk::mesh::FieldBase *myFieldBase = stk::mesh::get_field_by_name(fieldName, metaData);
        ASSERT_TRUE(myFieldBase != NULL);
        EXPECT_TRUE( &nodeField == myFieldBase);
    }

    unlink(mesh_name.c_str());
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

    SolutionCases() : fieldsPerSolutionCase(numSolnCases) {};
    ~SolutionCases() {};

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
            stk::mesh::Field<double>& field =
                    meta.declare_field<stk::mesh::Field<double>>(rank,
                            fieldNames[j],  1u);
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

    stk::mesh::Field<double> *dispField = bulk.mesh_meta_data().get_field<stk::mesh::Field<double>>(rank, "displacement");
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
    stk::mesh::MetaData meta;
    stk::mesh::BulkData bulk(meta, MPI_COMM_WORLD, stk::mesh::BulkData::NO_AUTO_AURA);
    int stepNum = 1;
    double time = 1.0;
    stk::io::fill_mesh_save_step_info(filename, bulk, stepNum, time);

    for(const std::string& field_name : fieldnames)
    {
        stk::mesh::FieldBase* field = bulk.mesh_meta_data().get_field(rank, field_name);
        ThrowRequireWithSierraHelpMsg(field!=nullptr);
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
        SolutionCases solnCases = setup_solution_cases();

        stk::mesh::PartVector solutionCasePartVector;
        setup_parts_associated_with_solution_cases(solnCases, get_meta(), rank, solutionCasePartVector);

        setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA);

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
        }

        verify_acceleration_is_not_on_entities(get_bulk(), rank);
    }
};

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

class LateFieldFixtureNoTest : public stk::unit_test_util::MeshFixtureNoTest
{
protected:
  LateFieldFixtureNoTest() {}
  LateFieldFixtureNoTest(unsigned spatial_dim) : MeshFixtureNoTest(spatial_dim) {}
  LateFieldFixtureNoTest(unsigned spatial_dim, const std::vector<std::string>& entityRankNames)
    : MeshFixtureNoTest(spatial_dim,entityRankNames)  {}

  stk::mesh::FieldDataManager * fieldDataManager = nullptr;

  virtual ~LateFieldFixtureNoTest()
  {
    delete bulkData;
    delete metaData;
    delete fieldDataManager;
    bulkData = nullptr;
    metaData = nullptr;
    fieldDataManager = nullptr;
  }
};


class LateFieldFixture : public LateFieldFixtureNoTest, public ::ngp_testing::Test
{
protected:
  virtual void allocate_bulk(stk::mesh::BulkData::AutomaticAuraOption auraOption, stk::mesh::FieldDataManager * fieldDataManager)
  {
    if (nullptr == metaData) {
      allocate_meta();
    }

    bulkData = new stk::mesh::BulkData(get_meta(), communicator, auraOption,
#ifdef SIERRA_MIGRATION
                                       false, // add_fmwk_data
#endif
                                       fieldDataManager);
  }

  virtual void setup_mesh(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption, stk::mesh::FieldDataManager * fieldDataManager)
  {
    allocate_bulk(auraOption, fieldDataManager);
    metaData->enable_late_fields();
    stk::io::fill_mesh(meshSpecification, *bulkData);
  }

  virtual void setup_empty_mesh(stk::mesh::BulkData::AutomaticAuraOption auraOption, stk::mesh::FieldDataManager * fieldDataManager)
  {
    allocate_bulk(auraOption, fieldDataManager);
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
    stk::mesh::Field<T> & field = get_meta().declare_field<stk::mesh::Field<T>>(rank, fieldName, numStates);
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
  void setup_add_late_first_field(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    // Note that we still have a nodal coordinates field

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    set_field_values_with_scale_factor(lateField, 1);
    expect_field_values_with_scale_factor(lateField, 1);
  }

  template <typename T>
  void setup_add_late_field(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_performance_of_early_field(stk::mesh::EntityRank rank, int numberOfFields, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    for (int i = 0; i < numberOfFields; ++i) {
      stk::mesh::Field<T> & earlyField = declare_field<T>("early_field" + std::to_string(i), rank);
      put_field(earlyField, get_meta().universal_part());
    }
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    get_bulk().deactivate_field_updating();  // Shut down extra bucket churn, to mimic real file I/O
    stk::io::fill_mesh("generated:100x100x100", get_bulk());
  }

  template <typename T>
  void setup_performance_of_late_field(stk::mesh::EntityRank rank, int numberOfFields, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    setup_mesh("generated:100x100x100", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    for (int i = 0; i < numberOfFields; ++i) {
      stk::mesh::Field<T> & lateField = declare_field<T>("late_field" + std::to_string(i), rank);
      put_field(lateField, get_meta().universal_part());
    }
  }

  template <typename T>
  void setup_performance_of_single_late_field(stk::mesh::EntityRank rank, int numberOfEarlyFields, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    for (int i = 0; i < numberOfEarlyFields; ++i) {
      stk::mesh::Field<T> & earlyField = declare_field<T>("early_field" + std::to_string(i), rank);
      put_field(earlyField, get_meta().universal_part());
    }
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    get_bulk().deactivate_field_updating();  // Shut down extra bucket churn, to mimic real file I/O
    stk::io::fill_mesh("generated:100x100x100", get_bulk());

    const double startTime = stk::wall_time();

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());

    const double stopTime = stk::wall_time();
    std::cout << "Time for adding single late field: " << stopTime-startTime << " s" << std::endl;
  }


  template <typename T>
  void setup_add_late_field_multiple_buckets(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);

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
  void setup_add_late_field_multiple_duplicate_put_field(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);
    put_field(lateField, get_meta().universal_part());
    set_field_values_with_scale_factor(lateField, 2);

    put_field(lateField, get_meta().universal_part());

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_late_field_multiple_different_put_field(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
    set_field_values_with_scale_factor(earlyField, 1);

    stk::mesh::Field<T> & lateField = declare_field<T>("late_field", rank);

    put_field(lateField, *get_meta().get_part("block_1"));
    put_field(lateField, *get_meta().get_part("block_2"));
    set_field_values_with_scale_factor(lateField, 2);

    expect_field_values_with_scale_factor(earlyField, 1);
    expect_field_values_with_scale_factor(lateField, 2);
  }

  template <typename T>
  void setup_add_two_late_fields_sequential(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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
  void setup_add_two_late_fields_interleaved(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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
  void setup_add_two_late_fields_out_of_order(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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
  void setup_add_two_late_fields_different_type_out_of_order(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T1> & earlyField = declare_field<T1>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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
  void setup_add_two_late_fields_different_rank_out_of_order(stk::mesh::EntityRank rank1, stk::mesh::EntityRank rank2, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", stk::topology::ELEM_RANK);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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
  void setup_add_early_field_to_late_part(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    const bool addIoPartAttribute = false;
    create_part("block_1", stk::topology::ELEM_RANK, addIoPartAttribute);
    create_part("block_2", stk::topology::ELEM_RANK);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, *get_meta().get_part("block_1"));
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);

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
  void setup_add_late_field_to_late_part(stk::mesh::EntityRank rank, stk::mesh::FieldDataManager * fieldDataManager = nullptr) {
    const bool addIoPartAttribute = false;
    create_part("block_1", stk::topology::ELEM_RANK, addIoPartAttribute);
    stk::mesh::Field<T> & earlyField = declare_field<T>("early_field", rank);
    put_field(earlyField, get_meta().universal_part());
    setup_mesh("generated:1x1x2", stk::mesh::BulkData::NO_AUTO_AURA, fieldDataManager);
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

TEST_F(LateFieldFixture, DISABLED_performanceOfEarlyField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  const double startTime = stk::wall_time();
  setup_performance_of_early_field<int>(stk::topology::NODE_RANK, 64);
  const double stopTime = stk::wall_time();
  std::cout << "Time for early field registration: " << stopTime-startTime << " s" << std::endl;
}

TEST_F(LateFieldFixture, DISABLED_performanceOfLateField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  const double startTime = stk::wall_time();
  setup_performance_of_late_field<int>(stk::topology::NODE_RANK, 64);
  const double stopTime = stk::wall_time();
  std::cout << "Time for late field registration: " << stopTime-startTime << " s" << std::endl;
}

TEST_F(LateFieldFixture, DISABLED_performanceOfSingleLateField)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  setup_performance_of_single_late_field<int>(stk::topology::NODE_RANK, 64);
}

TEST_F(LateFieldFixture, addLateIntFirstElementFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_first_field<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntNodalFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntElementFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleBucketsContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_multiple_buckets<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDuplicatePutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDuplicatePutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_multiple_duplicate_put_field<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntNodalField_multipleDifferentPutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateIntElementField_multipleDifferentPutFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_multiple_different_put_field<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_sequentialContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_sequential<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_sequentialContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_sequential<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_interleavedContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_interleaved<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_interleavedContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_interleaved<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_out_of_order<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_out_of_order<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateIntAndDoubleElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_different_type_out_of_order<int, double>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleNodalFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateShortAndDoubleElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_different_type_out_of_order<short, double>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addTwoLateNodalAndElementFields_outOfOrderContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_two_late_fields_different_rank_out_of_order<int>(stk::topology::NODE_RANK, stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addEarlyNodalFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_early_field_to_late_part<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addEarlyElementFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_early_field_to_late_part<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateNodalFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_to_late_part<int>(stk::topology::NODE_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, addLateElementFieldToLatePartContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 2) return;
  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_add_late_field_to_late_part<int>(stk::topology::ELEM_RANK, fieldDataManager);
}

TEST_F(LateFieldFixture, DISABLED_performanceOfEarlyFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  const double startTime = stk::wall_time();
  setup_performance_of_early_field<int>(stk::topology::NODE_RANK, 64, fieldDataManager);
  const double stopTime = stk::wall_time();
  std::cout << "Time for early field registration: " << stopTime-startTime << " s" << std::endl;
}

TEST_F(LateFieldFixture, DISABLED_performanceOfLateFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  const double startTime = stk::wall_time();
  setup_performance_of_late_field<int>(stk::topology::NODE_RANK, 64, fieldDataManager);
  const double stopTime = stk::wall_time();
  std::cout << "Time for late field registration: " << stopTime-startTime << " s" << std::endl;
}

TEST_F(LateFieldFixture, DISABLED_performanceOfSingleLateFieldContiguous)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 16) return;

  fieldDataManager = new stk::mesh::ContiguousFieldDataManager;
  setup_performance_of_single_late_field<int>(stk::topology::NODE_RANK, 64, fieldDataManager);
}

} //namespace <anonymous>

