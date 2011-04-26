/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <stk_io/util/Gmesh_STKmesh_Fixture.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#  include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>

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
  stk::mesh::fem::FEMMetaData & fem_meta = stk::mesh::fem::FEMMetaData::get( fixture.getMetaData() );

  fixture.commit();

  const std::vector<std::string> & sideset_names = fixture.getSidesetNames();
  STKUNIT_ASSERT_EQUAL( num_surf, sideset_names.size() );

  for( size_t i = 0; i < num_surf; ++i ) {
    std::string surf_name =  (std::string)"surface_" + Ioss::Utils::to_string(i+1);
    STKUNIT_ASSERT(surf_name == sideset_names[i]);
  }

  std::vector<size_t> num_surf_elem(3);
  num_surf_elem[0] = num_y * num_z;
  num_surf_elem[1] = num_x * num_z;
  num_surf_elem[2] = num_x * num_y;

  for( size_t i = 0; i < num_surf/2; ++i )
  {
    STKUNIT_ASSERT_EQUAL( num_surf_elem[i], fixture.getSurfElemCount(2*i) );
    STKUNIT_ASSERT_EQUAL( num_surf_elem[i], fixture.getSurfElemCount(2*i+1) );
  }

  const size_t total_side_count = 2 * (num_surf_elem[0]+num_surf_elem[1]+num_surf_elem[2]);
  STKUNIT_ASSERT_EQUAL( total_side_count, fixture.getSideCount() );

  const size_t total_elem_count = num_x * num_y * num_z;
  STKUNIT_ASSERT_EQUAL( total_elem_count, fixture.getElemCount() );

  const size_t total_node_count = (num_x+1) * (num_y+1) * (num_z+1);
  STKUNIT_ASSERT_EQUAL( total_node_count, fixture.getNodeCount() );
  
  // Needed to test field data
  stk::mesh::Field<double,stk::mesh::Cartesian> * coord_field =
    fixture.getMetaData().get_field<stk::mesh::Field<double,stk::mesh::Cartesian> >("coordinates");
  STKUNIT_ASSERT( coord_field );

  // All side buckets
  const std::vector<stk::mesh::Bucket*> & all_side_buckets = fixture.getBulkData().buckets( fem_meta.side_rank() );

  std::vector<stk::mesh::Entity *> entities;

  const stk::mesh::PartVector & side_parts = fixture.getSideParts();
  STKUNIT_ASSERT_EQUAL( sideset_names.size(), side_parts.size() );

  for( size_t ifset = 0; ifset < side_parts.size(); ++ifset )
  {
    std::pair<int, double> expected = fixture.getSurfCoordInfo(ifset);

    stk::mesh::Selector selector = *side_parts[ifset];
    entities.clear();
    stk::mesh::get_selected_entities(selector, all_side_buckets, entities);
    STKUNIT_ASSERT_EQUAL( fixture.getSurfElemCount(ifset), entities.size() );

    for ( size_t i = 0 ; i < entities.size() ; ++i ) {
      stk::mesh::Entity & side = *entities[i] ;

      const CellTopologyData * cell_topology = stk::mesh::fem::get_cell_topology_new(side).getCellTopologyData();

      STKUNIT_ASSERT( cell_topology );

      stk::mesh::PairIterRelation rel = side.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );

      STKUNIT_ASSERT_EQUAL( cell_topology->node_count, rel.size() );

      for ( unsigned j = 0 ; j < cell_topology->node_count ; ++j ) 
      {
        stk::mesh::Entity & rel_node = *rel[j].entity();
        double * coords = stk::mesh::field_data(*coord_field, rel_node);
	STKUNIT_ASSERT( coords );
        //std::cerr << "( " << coords[0] << ", " << coords[1] << ", " << coords[2] << ")" << std::endl;

        STKUNIT_ASSERT_DOUBLE_EQ(expected.second, coords[expected.first]);
      }
      //std::cerr << std::endl;
    }
  }
}

