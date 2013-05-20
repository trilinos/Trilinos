/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/fixtures/GearsFixture.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <sstream>

namespace {

// Test that our debug printing works

STKUNIT_UNIT_TEST( UnitTestDebugDump, MetaData )
{
  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_WORLD, 1,
                                            stk::mesh::fixtures::GearParams(0.01, 0.4, 1.5, -0.4, 0.4));
  fixture.meta_data.commit();

  // Doesn't check anything, but at least makes sure it passes
  std::ostringstream out;
  fixture.meta_data.dump_all_meta_info(out);
}

STKUNIT_UNIT_TEST( UnitTestDebugDump, BulkData )
{
  const unsigned NX = 3;
  const unsigned NY = 1;
  const unsigned NZ = 1;
  stk::mesh::fixtures::HexFixture hf(MPI_COMM_WORLD, NX, NY, NZ);
  hf.m_fem_meta.commit();
  hf.generate_mesh();

  // Doesn't check anything, but at least makes sure it passes
  std::ostringstream out;
  hf.m_bulk_data.dump_all_mesh_info(out);

  // Uncomment to see output
  //std::cout << out.str() << std::endl;
}

} // namespace
