/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <iostream>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <stk_mesh/fixtures/GearsFixture.hpp>
#include <stk_mesh/fixtures/Gear.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <Shards_BasicTopologies.hpp>

#define PI     3.14159265358979
#define TWO_PI 6.28210184121061

namespace stk {
namespace mesh {
namespace fixtures {

typedef shards::Hexahedron<8> Hex8 ;
typedef shards::Wedge<6>      Wedge6 ;




GearsFixture::GearsFixture( ParallelMachine pm, size_t num_gears)
  : NUM_GEARS(num_gears)
    , meta_data( TopologicalMetaData::entity_rank_names(SpatialDimension) )
    , bulk_data( meta_data , pm )
    , top_data( meta_data, SpatialDimension )
    , cylindrical_coord_part( meta_data.declare_part("cylindrical_coord_part", SpatialDimension))
    , hex_part( top_data.declare_part<Hex8>("hex8_part"))
    , wedge_part( top_data.declare_part<Wedge6>("wedge6_part"))
    , cartesian_coord_field( meta_data.declare_field<CartesianField>("coordinates", 1))
    , displacement_field( meta_data.declare_field<CartesianField>("displacement", 2))
    , translation_field( meta_data.declare_field<CartesianField>("translation", 2))
    , cylindrical_coord_field( meta_data.declare_field<CylindricalField>("cylindrical_coordinates", 2))
    , m_gears()
  {

    put_field(
        cartesian_coord_field,
        top_data.node_rank,
        meta_data.universal_part(),
        SpatialDimension
        );

    put_field(
        displacement_field,
        top_data.node_rank,
        meta_data.universal_part(),
        SpatialDimension
        );

    put_field(
        translation_field,
        top_data.node_rank,
        cylindrical_coord_part,
        SpatialDimension
        );

    put_field(
        cylindrical_coord_field,
        NodeRank,
        cylindrical_coord_part,
        SpatialDimension
        );

    m_gears.resize(NUM_GEARS);

    for ( size_t i = 0; i < NUM_GEARS; ++i) {
      std::ostringstream oss;
      oss << "Gear_" << i ;
      m_gears[i] = new Gear (
          meta_data,
          bulk_data,
          meta_data.declare_part(oss.str(),SpatialDimension),
          cylindrical_coord_part,
          hex_part,
          wedge_part,
          cartesian_coord_field,
          displacement_field,
          translation_field,
          cylindrical_coord_field
          );
    }
  }

GearsFixture::~GearsFixture() {

  for( std::vector<Gear *>::iterator i = m_gears.begin();
       i != m_gears.end();
       ++i )
  {
    delete *i;
    *i = NULL;
  }
}

void GearsFixture::generate_mesh() {

  bulk_data.modification_begin();

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  //create the gears on a line
  for( size_t i = 0; i < m_gears.size(); ++i) {
    if (( (i*p_size)/m_gears.size()) == p_rank) {
      Gear & gear = get_gear(i);
      gear.generate_gear();

      GearData data;
    }
  }

  bulk_data.modification_end();

  communicate_model_fields();

}


void GearsFixture::communicate_model_fields() {

  //copy the field data to the aura nodes

  std::vector< const FieldBase *> fields;

  fields.push_back(& cartesian_coord_field);
  fields.push_back(& translation_field);
  fields.push_back(& cylindrical_coord_field);

  communicate_field_data(bulk_data.shared_aura(), fields);
}

} // fixtures
} // mesh
} // stk



