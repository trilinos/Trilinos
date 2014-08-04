/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <stk_io/StkMeshIoBroker.hpp>

#include <stk_mesh/base/BulkData.hpp>

TEST(UnitTestNodeBucketsHaveValidTopology, testUnit)
{
  std::string generated_mesh("generated:2x3x4");

  stk::ParallelMachine comm = MPI_COMM_WORLD;
  stk::io::StkMeshIoBroker meshReader(comm);
  meshReader.add_mesh_database(generated_mesh, stk::io::READ_MESH);

  meshReader.create_input_mesh();
  meshReader.populate_bulk_data();

  const stk::mesh::BulkData& stkmesh = meshReader.bulk_data();

  const stk::mesh::BucketVector& nodeBuckets = stkmesh.buckets(stk::topology::NODE_RANK);

  for(size_t i=0; i<nodeBuckets.size(); ++i) {
    const stk::mesh::Bucket& bucket = *nodeBuckets[i];

    EXPECT_EQ(stk::topology::NODE, bucket.topology());
  }
}

