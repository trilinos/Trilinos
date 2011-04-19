/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <cmath>
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
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fixtures/Gear.hpp>

namespace {

const stk::mesh::EntityRank NODE_RANK = stk::mesh::fem::FEMMetaData::NODE_RANK;

}

namespace stk {
namespace mesh {
namespace fixtures {

Gear::Gear(
    fem::FEMMetaData & meta,
    BulkData & bulk,
    Part & gear,
    Part & arg_cylindrical_coord_part,
    Part & hex,
    Part & wedge,
    CartesianField   & arg_cartesian_coord_field,
    CartesianField   & arg_displacement_field,
    CartesianField   & arg_translation_field,
    CylindricalField & arg_cylindrical_coord_field,
    double arg_element_size,
    double arg_radius_min,
    double arg_radius_max,
    double arg_height_min,
    double arg_height_max
    ) :
      element_size(arg_element_size)
    , rad_min(arg_radius_min)
    , rad_max(arg_radius_max)
    , height_min(arg_height_min)
    , height_max(arg_height_max)

    , angle_num(static_cast<size_t>(TWO_PI/element_size))
    , rad_num(static_cast<size_t>(2 +(rad_max - rad_min)/element_size))
    , height_num(static_cast<size_t>(2 +(height_max - height_min)/element_size))

    , angle_increment( TWO_PI / static_cast<double>(angle_num))
    , rad_increment( (rad_max-rad_min) / static_cast<double>(rad_num-1))
    , height_increment( (height_max-height_min) / static_cast<double>(height_num-1))

    , num_elements( angle_num * (rad_num-1) * (height_num-1) )
    , num_nodes( angle_num * (rad_num) * (height_num) )

    , meta_data(meta)
    , bulk_data(bulk)

    , gear_part(gear)
    , cylindrical_coord_part(arg_cylindrical_coord_part)
    , hex_part(hex)
    , wedge_part(wedge)

    , cartesian_coord_field(arg_cartesian_coord_field)
    , displacement_field(arg_displacement_field)
    , translation_field(arg_translation_field)
    , cylindrical_coord_field(arg_cylindrical_coord_field)
{}

//
//-----------------------------------------------------------------------------
//

void Gear::populate_fields(stk::mesh::FieldState state) {

  //setup the cylindrical_coord_field on the hex nodes
  for ( size_t ir = 0 ; ir < rad_num-1; ++ir ) {
    const double rad = rad_min + rad_increment * ir ;

    for ( size_t ia = 0 ; ia < angle_num; ++ia ) {
      const double angle = angle_increment * ia ;

      for ( size_t iz = 0 ; iz < height_num; ++iz ) {
        const double height = height_min + height_increment * iz ;

        Entity & node = get_node(iz,ir,ia);

        double * const cylindrical_data = field_data( cylindrical_coord_field , node );
        double * const translation_data = field_data( translation_field , node );
        double * const cartesian_data = field_data( cartesian_coord_field , node );
        double * const displacement_data = field_data( displacement_field.field_of_state(state) , node );

        cylindrical_data[0] = rad ;
        cylindrical_data[1] = angle ;
        cylindrical_data[2] = height;

        translation_data[0] = 0;
        translation_data[1] = 0;
        translation_data[2] = 0;

        displacement_data[0] = 0;
        displacement_data[1] = 0;
        displacement_data[2] = 0;

        cartesian_data[0] = rad * std::cos(angle);
        cartesian_data[1] = rad * std::sin(angle);
        cartesian_data[2] = height;
      }
    }
  }

  //setup the cylindrical_coord_field on the wedge nodes
  {
    size_t ir = rad_num-1;
    //const double rad = rad_min + rad_increment * ir ;
    const double rad = 1.1*rad_max;

    for ( size_t ia = 0 ; ia < angle_num; ++ia ) {
      const double angle = angle_increment * (ia + ia +1.0)/2.0;

      for ( size_t iz = 0 ; iz < height_num; ++iz ) {
        const double height = height_min + height_increment * iz ;

        Entity & node = get_node(iz,ir,ia);

        double * const cylindrical_data = field_data( cylindrical_coord_field , node );
        double * const translation_data = field_data( translation_field , node );
        double * const cartesian_data = field_data( cartesian_coord_field , node );
        double * const displacement_data = field_data( displacement_field.field_of_state(state) , node );

        cylindrical_data[0] = rad ;
        cylindrical_data[1] = angle ;
        cylindrical_data[2] = height;

        translation_data[0] = 0;
        translation_data[1] = 0;
        translation_data[2] = 0;

        displacement_data[0] = 0;
        displacement_data[1] = 0;
        displacement_data[2] = 0;

        cartesian_data[0] = rad * std::cos(angle);
        cartesian_data[1] = rad * std::sin(angle);
        cartesian_data[2] = height;
      }
    }
  }
}


//
//-----------------------------------------------------------------------------
//

void Gear::generate_gear()
{
  const stk::mesh::EntityRank element_rank = meta_data.element_rank();

  std::vector<size_t> requests(meta_data.entity_rank_count(), 0);
  requests[NODE_RANK]     = num_nodes;
  requests[element_rank] = num_elements;

  // Parallel collective call:
  bulk_data.generate_new_entities(requests, gear_entities);

  //setup hex elements
  {
    std::vector<Part*> add_parts, remove_parts ;
    add_parts.push_back( & cylindrical_coord_part );
    add_parts.push_back( & gear_part );
    add_parts.push_back( & hex_part );

    for ( size_t ir = 0 ; ir < rad_num -2 ; ++ir ) {
      for ( size_t ia = 0 ; ia < angle_num; ++ia ) {
        for ( size_t iz = 0 ; iz < height_num -1 ; ++iz ) {

          Entity & elem = get_element(iz, ir, ia);
          bulk_data.change_entity_parts(elem, add_parts, remove_parts);

          const size_t ia_1 = ( ia + 1 ) % angle_num ;
          const size_t ir_1 = ir + 1 ;
          const size_t iz_1 = iz + 1 ;

          Entity * node[8] ;

          node[0] = & get_node(iz  , ir  , ia_1 );
          node[1] = & get_node(iz  , ir  , ia   );
          node[2] = & get_node(iz_1, ir  , ia   );
          node[3] = & get_node(iz_1, ir  , ia_1 );
          node[4] = & get_node(iz  , ir_1, ia_1 );
          node[5] = & get_node(iz  , ir_1, ia   );
          node[6] = & get_node(iz_1, ir_1, ia   );
          node[7] = & get_node(iz_1, ir_1, ia_1 );

          for ( size_t j = 0 ; j < 8 ; ++j ) {
            bulk_data.declare_relation( elem , * node[j] , j );
          }
        }
      }
    }
  }

  //setup wedges elements
  {
    std::vector<Part*> add_parts, remove_parts ;
    add_parts.push_back( & cylindrical_coord_part );
    add_parts.push_back( & gear_part );
    add_parts.push_back( & wedge_part );

    size_t ir = rad_num-2 ;
    for ( size_t ia = 0 ; ia < angle_num; ++ia ) {
      for ( size_t iz = 0 ; iz < height_num -1 ; ++iz ) {

        Entity & elem = get_element(iz, ir, ia);
        bulk_data.change_entity_parts(elem, add_parts, remove_parts);

        const size_t ia_1 = ( ia + 1 ) % angle_num ;
        const size_t ir_1 = ir + 1 ;
        const size_t iz_1 = iz + 1 ;

        Entity * node[6] ;

        node[0] = & get_node(iz  , ir  , ia_1 );
        node[1] = & get_node(iz  , ir  , ia   );
        node[2] = & get_node(iz  , ir_1, ia   );
        node[3] = & get_node(iz_1, ir  , ia_1 );
        node[4] = & get_node(iz_1, ir  , ia   );
        node[5] = & get_node(iz_1, ir_1, ia   );

        for ( size_t j = 0 ; j < 6 ; ++j ) {
          bulk_data.declare_relation( elem , * node[j] , j );
        }
      }
    }
  }

  //cylindrical and cartesian coordinates
  //are 2 state fields.  Need to update
  //both states at construction
  populate_fields(stk::mesh::StateOld);
  populate_fields(stk::mesh::StateNew);

}

//
//-----------------------------------------------------------------------------
//

void Gear::move( const GearMovement & data) {

  enum { Node = 0 };

  Selector select = gear_part & cylindrical_coord_part & (meta_data.locally_owned_part() | meta_data.globally_shared_part());

  BucketVector all_node_buckets = bulk_data.buckets(NODE_RANK);

  BucketVector node_buckets;

  get_buckets(
      select,
      all_node_buckets,
      node_buckets
      );

  for (BucketVector::iterator b_itr = node_buckets.begin();
      b_itr != node_buckets.end();
      ++b_itr)
  {
    Bucket & b = **b_itr;
    BucketArray<CylindricalField> cylindrical_data( cylindrical_coord_field, b);  // ONE STATE
    BucketArray<CartesianField>   translation_data( translation_field, b); // ONE STATE
    const BucketArray<CartesianField>   old_coordinate_data( cartesian_coord_field, b); // ONE STATE
    BucketArray<CartesianField>   new_displacement_data( displacement_field.field_of_state(stk::mesh::StateNew), b); // TWO STATE

    double new_coordinate_data[3] = {0,0,0};
    for (size_t i = 0; i < b.size(); ++i) {
      int index = i;

      const double   radius = cylindrical_data(0,index);
            double & angle  = cylindrical_data(1,index);
      const double   height = cylindrical_data(2,index);


      angle += data.rotation;

      if ( angle < 0.0) {
        angle += TWO_PI;
      }
      else if ( angle > TWO_PI) {
        angle -= TWO_PI;
      }

      translation_data(0,index) += data.x;
      translation_data(1,index) += data.y;
      translation_data(2,index) += data.z;

      new_coordinate_data[0] = translation_data(0,index) + radius * std::cos(angle);
      new_coordinate_data[1] = translation_data(1,index) + radius * std::sin(angle);
      new_coordinate_data[2] = translation_data(2,index) + height;

      new_displacement_data(0,index) = new_coordinate_data[0] - old_coordinate_data(0,index);
      new_displacement_data(1,index) = new_coordinate_data[1] - old_coordinate_data(1,index);
      new_displacement_data(2,index) = new_coordinate_data[2] - old_coordinate_data(2,index);

    }
  }
}

} // fixtures
} // mesh
} // stk
