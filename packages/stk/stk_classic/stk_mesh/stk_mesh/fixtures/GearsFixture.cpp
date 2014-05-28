/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
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
#include <set>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fixtures/GearsFixture.hpp>
#include <stk_mesh/fixtures/Gear.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <Shards_BasicTopologies.hpp>

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

const unsigned ONE_STATE = 1;
const unsigned TWO_STATE = 2;

typedef shards::Hexahedron<8> Hex8 ;
typedef shards::Wedge<6>      Wedge6 ;

} // namespace

namespace stk {
namespace mesh {
namespace fixtures {

GearsFixture::GearsFixture( ParallelMachine pm, size_t num_gears, GearParams gear_params)
  : NUM_GEARS(num_gears)
  , meta_data( SpatialDimension )
  , bulk_data( fem::FEMMetaData::get_meta_data(meta_data) , pm )
  , element_rank( meta_data.element_rank() )
  , cylindrical_coord_part( meta_data.declare_part("cylindrical_coord_part", element_rank))
  , hex_part( fem::declare_part<Hex8>(meta_data, "hex8_part"))
  , wedge_part( fem::declare_part<Wedge6>(meta_data, "wedge6_part"))
  , cartesian_coord_field( meta_data.declare_field<CartesianField>("coordinates", ONE_STATE))
  , displacement_field( meta_data.declare_field<CartesianField>("displacement", TWO_STATE))
  , translation_field( meta_data.declare_field<CartesianField>("translation", ONE_STATE))
  , cylindrical_coord_field( meta_data.declare_field<CylindricalField>("cylindrical_coordinates", ONE_STATE))
  , m_gears()
  {

    put_field(
        cartesian_coord_field,
        NODE_RANK,
        meta_data.universal_part(),
        SpatialDimension
        );

    put_field(
        displacement_field,
        NODE_RANK,
        meta_data.universal_part(),
        SpatialDimension
        );

    put_field(
        translation_field,
        NODE_RANK,
        cylindrical_coord_part,
        SpatialDimension
        );

    put_field(
        cylindrical_coord_field,
        NODE_RANK,
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
          cylindrical_coord_field,
          gear_params.element_size,
          gear_params.radius_min,
          gear_params.radius_max,
          gear_params.height_min,
          gear_params.height_max
          );
  }
}

GearsFixture::~GearsFixture()
{
  for( std::vector<Gear *>::iterator i = m_gears.begin();
       i != m_gears.end();
       ++i )
  {
    delete *i;
    *i = NULL;
  }
}

void GearsFixture::generate_mesh() {

  // Parallel collective call:
  bulk_data.modification_begin();

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  //create the gears on a line
  for( size_t i = 0; i < m_gears.size(); ++i) {
    if (( (i*p_size)/m_gears.size()) == p_rank) {
      Gear & gear = get_gear(i);
      gear.generate_gear();
    } else {
      // Parallel synchronization:
      std::vector<size_t> empty_requests(meta_data.entity_rank_count(), 0);
      EntityVector empty_entities;
      bulk_data.generate_new_entities(empty_requests, empty_entities);
    }
  }

  // Parallel collective call:
  bulk_data.modification_end();

  // Parallel collective call:
  communicate_model_fields();

  for ( size_t i = 0 ; i < m_gears.size() ; ++i ) {
    // Parallel collective call:
    distribute_gear_across_processors(get_gear(i),cylindrical_coord_field);
  }

}

void GearsFixture::communicate_model_fields()
{
  //copy the field data to the aura nodes

  std::vector< const FieldBase *> fields;

  fields.push_back(& cartesian_coord_field);
  fields.push_back(& translation_field);
  fields.push_back(& cylindrical_coord_field);
  fields.push_back(& displacement_field.field_of_state(stk::mesh::StateNew));
  fields.push_back(& displacement_field.field_of_state(stk::mesh::StateOld));

  // Parallel collective call:
  communicate_field_data(bulk_data.shared_aura(), fields);
}


double scale_angle_2pi(double angle) {
  while ( angle < 0.0 ) {
    angle += TWO_PI;
  }
  while ( angle >= TWO_PI) {
    angle -= TWO_PI;
  }
  return angle;
}


void select_nodal_data(
    GearsFixture::CylindricalField & cylindrical_coord_field, 
    Entity & element, 
    double & radius, 
    double & angle, 
    double & height
    )
{
  radius = 0.0;
  angle = TWO_PI;
  height = 0.0;
  PairIterRelation node_relations = element.relations(NODE_RANK);
  int numNodes = node_relations.second - node_relations.first;
  for ( ; node_relations.first != node_relations.second ; ++(node_relations.first) ) {
    Entity * node = node_relations.first->entity();
    EntityArray<GearsFixture::CylindricalField> cylindrical_data( cylindrical_coord_field, *node);
    radius += cylindrical_data(0);
    angle  = std::min(angle,cylindrical_data(1));
    height += cylindrical_data(2);
  }
  radius /= numNodes;
  height /= numNodes;
}

// Parallel collective call:
void distribute_gear_across_processors(Gear & gear, GearsFixture::CylindricalField & cylindrical_coord_field)
{
  BulkData & bulk_data = gear.bulk_data;

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();
  
  EntityProcVec elements_to_change_owner;

  Selector locally_owned = gear.meta_data.locally_owned_part();
  if (p_rank == 0) {
    BucketVector all_elements;
    stk::mesh::get_buckets(locally_owned,bulk_data.buckets(gear.meta_data.element_rank()),all_elements);
    std::set<Entity *> node_set; // First come first serve nodal movement.
    for (BucketVector::iterator it = all_elements.begin() ; it != all_elements.end() ; ++it) {
      Bucket & b = **it;
      for (size_t i=0 ; i<b.size() ; ++i) {
        Entity * element = &b[i];
        double radius = 0.0;
        double angle  = 0.0;
        double height = 0.0;
        select_nodal_data(cylindrical_coord_field, *element,radius,angle,height);
        unsigned destination_processor_rank = destination_processor(gear,radius,angle,height,p_rank,p_size);
        elements_to_change_owner.push_back(EntityProc(element,destination_processor_rank));
        // Now add all related nodes to list to move to this processor:
        PairIterRelation node_relations = element->relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
        for ( ; node_relations.first != node_relations.second ; ++(node_relations.first) ) {
          Entity * node = node_relations.first->entity();
          if (node_set.count(node)==0) {
            elements_to_change_owner.push_back(EntityProc(node,destination_processor_rank));
            node_set.insert(node);
          }
        }
      }
    }
  }

  // Parallel collective call:
  bulk_data.modification_begin();
  bulk_data.change_entity_owner(elements_to_change_owner);
  // Parallel collective call:
  bulk_data.modification_end();

  // Print out how many ended up on each processor:
  //{
  //  BucketVector local_elements;
  //  stk::mesh::get_buckets(locally_owned,bulk_data.buckets(Element),local_elements);
  //  BucketVector local_nodes;
  //  stk::mesh::get_buckets(locally_owned,bulk_data.buckets(Node),local_nodes);
  //  size_t element_count = 0;
  //  for (BucketVector::iterator it = local_elements.begin() ; it != local_elements.end() ; ++it) {
  //    Bucket & b = **it;
  //    element_count += b.size();
  //  }
  //  size_t node_count = 0;
  //  for (BucketVector::iterator it = local_nodes.begin() ; it != local_nodes.end() ; ++it) {
  //    Bucket & b = **it;
  //    node_count += b.size();
  //  }
  //  std::cout << "Proc " << p_rank << ": element count = " << element_count << " node count = " << node_count << std::endl;
  //}
  
}


double floor0(double value)
{
  double result = std::floor( std::fabs( value ) );
  return (value < 0.0) ? -result : result;
}


void scale_p_rank(unsigned & p_rank, unsigned p_size)
{
  if (p_rank >= p_size) {
    p_rank = p_size-1;
  }
}

unsigned destination_processor(const Gear & gear, double rad, double angle, double height, unsigned p_rank, unsigned p_size)
{
  unsigned result = 0;
  // Distribute elements across angles: (not working perfectly yet)
  angle = scale_angle_2pi(angle);
  result = static_cast<unsigned>(floor0((angle/TWO_PI)*p_size));
  
  // Distribute elements across radius:
  //result = static_cast<unsigned>(floor0((rad-gear.rad_min)/(gear.rad_max-gear.rad_min)*p_size));

  // Distribute elements across height:
  //result = static_cast<unsigned>(floor0((height-gear.height_min)/(gear.height_max-gear.height_min)*p_size));
  
  // Distribute elements randomly:
  //result = std::rand() % p_size;

  // Distribute in round-robin fashion:
  //static unsigned dest_proc = 0;
  //result = dest_proc++ % p_size;

  // Distribute all to processor 2:
  //result = 1;

  scale_p_rank(result,p_size);
  
  //std::cout << "P"<<p_rank<<"("<<p_size<<"):  (r,t,z) = ("<<rad<<", "<<angle<<", "<<height<<"), dest = " << result << std::endl;

  return result;
}

} // fixtures
} // mesh
} // stk



