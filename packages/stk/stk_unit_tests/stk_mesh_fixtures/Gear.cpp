// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_mesh/fixtures/Gear.hpp>
#include <cmath>                        // for cos, sin
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Selector.hpp>   // for operator&, Selector, etc
#include <stk_mesh/base/Types.hpp>      // for BucketVector, EntityRank
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_topology/topology.hpp"    // for topology, etc



namespace {

const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

}

namespace stk {
namespace mesh {
namespace fixtures {

Gear::Gear(
    MetaData & meta,
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

void Gear::populate_fields(FieldState state) {

  //setup the cylindrical_coord_field on the hex nodes
  for ( size_t ir = 0 ; ir < rad_num-1; ++ir ) {
    const double rad = rad_min + rad_increment * ir ;

    for ( size_t ia = 0 ; ia < angle_num; ++ia ) {
      const double angle = angle_increment * ia ;

      for ( size_t iz = 0 ; iz < height_num; ++iz ) {
        const double height = height_min + height_increment * iz ;

        Entity node = get_node(iz,ir,ia);

        double * const cylindrical_data = stk::mesh::field_data( cylindrical_coord_field , node );
        double * const translation_data = stk::mesh::field_data( translation_field , node );
        double * const cartesian_data = stk::mesh::field_data( cartesian_coord_field , node );
        double * const displacement_data = stk::mesh::field_data( displacement_field.field_of_state(state) , node );

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

        Entity node = get_node(iz,ir,ia);

        double * const cylindrical_data = stk::mesh::field_data( cylindrical_coord_field , node );
        double * const translation_data = stk::mesh::field_data( translation_field , node );
        double * const cartesian_data = stk::mesh::field_data( cartesian_coord_field , node );
        double * const displacement_data = stk::mesh::field_data( displacement_field.field_of_state(state) , node );

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
  const EntityRank element_rank = stk::topology::ELEMENT_RANK;

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

          Entity elem = get_element(iz, ir, ia);
          bulk_data.change_entity_parts(elem, add_parts, remove_parts);

          const size_t ia_1 = ( ia + 1 ) % angle_num ;
          const size_t ir_1 = ir + 1 ;
          const size_t iz_1 = iz + 1 ;

          Entity node[8] ;

          node[0] = get_node(iz  , ir  , ia   );
          node[1] = get_node(iz_1, ir  , ia   );
          node[2] = get_node(iz_1, ir_1, ia   );
          node[3] = get_node(iz  , ir_1, ia   );
          node[4] = get_node(iz  , ir  , ia_1 );
          node[5] = get_node(iz_1, ir  , ia_1 );
          node[6] = get_node(iz_1, ir_1, ia_1 );
          node[7] = get_node(iz  , ir_1, ia_1 );

          for ( size_t j = 0 ; j < 8 ; ++j ) {
            bulk_data.declare_relation( elem , node[j] , j );
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

        Entity elem = get_element(iz, ir, ia);
        bulk_data.change_entity_parts(elem, add_parts, remove_parts);

        const size_t ia_1 = ( ia + 1 ) % angle_num ;
        const size_t ir_1 = ir + 1 ;
        const size_t iz_1 = iz + 1 ;

        Entity node[6] ;

        node[0] = get_node(iz  , ir  , ia_1 );
        node[1] = get_node(iz  , ir  , ia   );
        node[2] = get_node(iz  , ir_1, ia   );
        node[3] = get_node(iz_1, ir  , ia_1 );
        node[4] = get_node(iz_1, ir  , ia   );
        node[5] = get_node(iz_1, ir_1, ia   );

        for ( size_t j = 0 ; j < 6 ; ++j ) {
          bulk_data.declare_relation( elem , node[j] , j );
        }
      }
    }
  }

  //cylindrical and cartesian coordinates
  //are 2 state fields.  Need to update
  //both states at construction
  populate_fields(StateOld);
  populate_fields(StateNew);

}

//
//-----------------------------------------------------------------------------
//

void Gear::move( const GearMovement & data) {

  enum { Node = 0 };

  Selector select = gear_part & cylindrical_coord_part & (meta_data.locally_owned_part() | meta_data.globally_shared_part());

  BucketVector const& node_buckets = bulk_data.get_buckets(NODE_RANK, select);

  for (BucketVector::const_iterator b_itr = node_buckets.begin();
      b_itr != node_buckets.end();
      ++b_itr)
  {
    Bucket & b = **b_itr;
    double* cylindrical_data = stk::mesh::field_data( cylindrical_coord_field, b);  // ONE STATE
    double*   translation_data = stk::mesh::field_data( translation_field, b); // ONE STATE
    const double*   old_coordinate_data = stk::mesh::field_data( cartesian_coord_field, b); // ONE STATE
    double*   new_displacement_data = stk::mesh::field_data( displacement_field.field_of_state(StateNew), b); // TWO STATE

    double new_coordinate_data[3] = {0,0,0};
    for (size_t i = 0; i < b.size(); ++i) {
      int index = i;

      const double   radius = cylindrical_data[0+index*3];
            double & angle  = cylindrical_data[1+index*3];
      const double   height = cylindrical_data[2+index*3];


      angle += data.rotation;

      if ( angle < 0.0) {
        angle += TWO_PI;
      }
      else if ( angle > TWO_PI) {
        angle -= TWO_PI;
      }

      translation_data[0+index*3] += data.x;
      translation_data[1+index*3] += data.y;
      translation_data[2+index*3] += data.z;

      new_coordinate_data[0] = translation_data[0+index*3] + radius * std::cos(angle);
      new_coordinate_data[1] = translation_data[1+index*3] + radius * std::sin(angle);
      new_coordinate_data[2] = translation_data[2+index*3] + height;

      new_displacement_data[0+index*3] = new_coordinate_data[0] - old_coordinate_data[0+index*3];
      new_displacement_data[1+index*3] = new_coordinate_data[1] - old_coordinate_data[1+index*3];
      new_displacement_data[2+index*3] = new_coordinate_data[2] - old_coordinate_data[2+index*3];
    }
  }
}

} // fixtures
} // mesh
} // stk
