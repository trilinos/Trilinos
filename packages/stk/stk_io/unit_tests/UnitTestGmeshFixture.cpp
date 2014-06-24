/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <Ioss_Utils.h>                 // for Utils
#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <gtest/gtest.h>
#include <string>                       // for allocator, operator+, etc
#include <vector>                       // for vector
#include "gtest/gtest.h"                // for AssertHelper
#include "stk_mesh/base/Types.hpp"      // for PartVector




enum { SpaceDim = 3 };

TEST(UnitTestGmeshFixture, testUnit)
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
  ASSERT_EQ( num_surf, sideset_names.size() );

  for( size_t i = 0; i < num_surf; ++i ) {
    std::string surf_name =  (std::string)"surface_" + Ioss::Utils::to_string(i+1);
    ASSERT_TRUE(surf_name == sideset_names[i]);
  }

  // Needed to test field data
  stk::mesh::Field<double,stk::mesh::Cartesian> * coord_field =
    fixture.getMetaData().get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >(stk::topology::NODE_RANK, "coordinates");
  ASSERT_TRUE( coord_field );

  const stk::mesh::PartVector & side_parts = fixture.getSideParts();
  ASSERT_EQ( sideset_names.size(), side_parts.size() );
}

