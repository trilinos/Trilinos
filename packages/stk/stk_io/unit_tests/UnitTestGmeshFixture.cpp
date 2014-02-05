/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <Ioss_Utils.h>

#include <assert.h>

enum { SpaceDim = 3 };

STKUNIT_UNIT_TEST(UnitTestGmeshFixture, testUnit)
{
  const size_t num_x = 1;
  const size_t num_y = 2;
  const size_t num_z = 3;
  const size_t num_surf = 6;
  std::string config_mesh = Ioss::Utils::to_string(num_x) + "x" +
                            Ioss::Utils::to_string(num_y) + "x" +
                            Ioss::Utils::to_string(num_z) + "|sideset:xXyYzZ";
  stk::io::util::Gmesh_STKmesh_Fixture fixture(MPI_COMM_WORLD, config_mesh);

  fixture.commit();

  const std::vector<std::string> & sideset_names = fixture.getSidesetNames();
  STKUNIT_ASSERT_EQUAL( num_surf, sideset_names.size() );

  for( size_t i = 0; i < num_surf; ++i ) {
    std::string surf_name =  (std::string)"surface_" + Ioss::Utils::to_string(i+1);
    STKUNIT_ASSERT(surf_name == sideset_names[i]);
  }

  // Needed to test field data
  stk::mesh::Field<double,stk::mesh::Cartesian> * coord_field =
    fixture.getMetaData().get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");
  STKUNIT_ASSERT( coord_field );

  const stk::mesh::PartVector & side_parts = fixture.getSideParts();
  STKUNIT_ASSERT_EQUAL( sideset_names.size(), side_parts.size() );
}

