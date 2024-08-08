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

#include <stk_util/stk_config.h>        // for STK_HAS_MPI
#include <algorithm>                    // for min
#include <cmath>                        // for fabs, floor
#include <iostream>                     // for ostringstream, etc
#include <set>                          // for set
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/Types.hpp>      // for BucketVector, EntityProcVec, etc
#include <stk_unit_test_utils/stk_mesh_fixtures/Gear.hpp>   // for Gear, TWO_PI
#include <stk_unit_test_utils/stk_mesh_fixtures/GearsFixture.hpp>
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/FieldState.hpp"  // for FieldState::StateNew, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine
#include "stk_io/IossBridge.hpp"

namespace stk { namespace mesh { class FieldBase; } }

namespace {

const unsigned ONE_STATE = 1;
const unsigned TWO_STATE = 2;

} // namespace

namespace stk {
namespace mesh {
namespace fixtures {

GearsFixture::GearsFixture(ParallelMachine pm, size_t num_gears, GearParams gear_params)
  : NUM_GEARS(num_gears),
    bulk_data_ptr(stk::unit_test_util::build_mesh(SpatialDimension, pm)),
    bulk_data(*bulk_data_ptr),
    meta_data(bulk_data.mesh_meta_data()),
    cylindrical_coord_part( meta_data.declare_part("cylindrical_coord_part", stk::topology::ELEMENT_RANK)),
    hex_part( meta_data.declare_part_with_topology("hex8_part", stk::topology::HEX_8)),
    wedge_part( meta_data.declare_part_with_topology("wedge6_part", stk::topology::WEDGE_6)),
    m_gears()
{
  cartesian_coord_field   = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates", ONE_STATE);
  displacement_field      = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "displacement", TWO_STATE);
  translation_field       = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "translation", ONE_STATE);
  cylindrical_coord_field = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "cylindrical_coordinates", ONE_STATE);

  put_field_on_mesh(*cartesian_coord_field, meta_data.universal_part(), SpatialDimension, nullptr);
  stk::io::set_field_output_type(*cartesian_coord_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*displacement_field, meta_data.universal_part(), SpatialDimension, nullptr);
  stk::io::set_field_output_type(*displacement_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*translation_field, cylindrical_coord_part, SpatialDimension, nullptr);
  stk::io::set_field_output_type(*translation_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*cylindrical_coord_field, cylindrical_coord_part, SpatialDimension, nullptr);

  m_gears.resize(NUM_GEARS);

  for ( size_t i = 0; i < NUM_GEARS; ++i) {
    std::ostringstream oss;
    oss << "Gear_" << i;
    m_gears[i] = new Gear(meta_data,
                          bulk_data,
                          meta_data.declare_part(oss.str(),static_cast<EntityRank>(SpatialDimension)),
                          cylindrical_coord_part,
                          hex_part,
                          wedge_part,
                          *cartesian_coord_field,
                          *displacement_field,
                          *translation_field,
                          *cylindrical_coord_field,
                          gear_params.element_size,
                          gear_params.radius_min,
                          gear_params.radius_max,
                          gear_params.height_min,
                          gear_params.height_max);
  }
}

GearsFixture::~GearsFixture()
{
  for (std::vector<Gear *>::iterator i = m_gears.begin(); i != m_gears.end(); ++i) {
    delete *i;
    *i = nullptr;
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
    distribute_gear_across_processors(get_gear(i), *cylindrical_coord_field);
  }

}

void GearsFixture::communicate_model_fields()
{
  //copy the field data to the aura nodes

  std::vector< const FieldBase *> fields;

  fields.push_back(cartesian_coord_field);
  fields.push_back(translation_field);
  fields.push_back(cylindrical_coord_field);
  fields.push_back(&displacement_field->field_of_state(stk::mesh::StateNew));
  fields.push_back(&displacement_field->field_of_state(stk::mesh::StateOld));

  // Parallel collective call:
#if defined( STK_HAS_MPI)
  communicate_field_data(bulk_data.aura_ghosting(), fields);
#endif
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


void select_nodal_data(const BulkData& mesh,
                       GearsFixture::CylindricalField & cylindrical_coord_field,
                       Entity element,
                       double & radius,
                       double & angle,
                       double & height)
{
  radius = 0.0;
  angle = TWO_PI;
  height = 0.0;

  int numNodes = mesh.num_nodes(element);
  Entity const *elem_nodes = mesh.begin_nodes(element);
  for (int i = 0; i < numNodes; ++i)
  {
    Entity node = elem_nodes[i];
    const MeshIndex& mi = mesh.mesh_index(node);
    const Bucket& bucket = *mi.bucket;
    double *cylindrical_data = stk::mesh::field_data(cylindrical_coord_field, bucket, mi.bucket_ordinal);
    radius += cylindrical_data[0];
    angle  = std::min(angle,cylindrical_data[1]);
    height += cylindrical_data[2];
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
    BucketVector const& all_elements = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_owned);
    std::set<Entity> node_set; // First come first serve nodal movement.
    for (BucketVector::const_iterator it = all_elements.begin() ; it != all_elements.end() ; ++it) {
      Bucket & b = **it;
      for (size_t i=0 ; i<b.size() ; ++i) {
        Entity element = b[i];
        double radius = 0.0;
        double angle  = 0.0;
        double height = 0.0;
        select_nodal_data(bulk_data, cylindrical_coord_field, element, radius, angle, height);
        unsigned destination_processor_rank = destination_processor(gear,radius,angle,height,p_rank,p_size);
        elements_to_change_owner.push_back(EntityProc(element,destination_processor_rank));

        // Now add all related nodes to list to move to this processor:
        Entity const *elem_nodes_j = bulk_data.begin_nodes(element);
        Entity const *elem_nodes_e = bulk_data.end_nodes(element);
        for ( ; elem_nodes_j != elem_nodes_e; ++elem_nodes_j)
        {
          Entity node = *elem_nodes_j;
          if (node_set.count(node)==0) {
            elements_to_change_owner.push_back(EntityProc(node,destination_processor_rank));
            node_set.insert(node);
          }
        }
      }
    }
  }

  // Parallel collective call:
  bulk_data.change_entity_owner(elements_to_change_owner);
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

  scale_p_rank(result,p_size);

  return result;
}

namespace simple_fields {

GearsFixture::GearsFixture(ParallelMachine pm, size_t num_gears, stk::mesh::fixtures::GearParams gear_params)
  : NUM_GEARS(num_gears),
    bulk_data_ptr(stk::unit_test_util::build_mesh(SpatialDimension, pm)),
    bulk_data(*bulk_data_ptr),
    meta_data(bulk_data.mesh_meta_data()),
    cylindrical_coord_part( meta_data.declare_part("cylindrical_coord_part", stk::topology::ELEMENT_RANK)),
    hex_part( meta_data.declare_part_with_topology("hex8_part", stk::topology::HEX_8)),
    wedge_part( meta_data.declare_part_with_topology("wedge6_part", stk::topology::WEDGE_6)),
    m_gears()
{
  cartesian_coord_field   = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "coordinates", ONE_STATE);
  displacement_field      = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "displacement", TWO_STATE);
  translation_field       = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "translation", ONE_STATE);
  cylindrical_coord_field = &meta_data.declare_field<double>(stk::topology::NODE_RANK, "cylindrical_coordinates", ONE_STATE);

  put_field_on_mesh(*cartesian_coord_field, meta_data.universal_part(), SpatialDimension, nullptr);
  stk::io::set_field_output_type(*cartesian_coord_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*displacement_field, meta_data.universal_part(), SpatialDimension, nullptr);
  stk::io::set_field_output_type(*displacement_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*translation_field, cylindrical_coord_part, SpatialDimension, nullptr);
  stk::io::set_field_output_type(*translation_field, stk::io::FieldOutputType::VECTOR_3D);

  put_field_on_mesh(*cylindrical_coord_field, cylindrical_coord_part, SpatialDimension, nullptr);

  m_gears.resize(NUM_GEARS);

  for ( size_t i = 0; i < NUM_GEARS; ++i) {
    std::ostringstream oss;
    oss << "Gear_" << i;
    m_gears[i] = new stk::mesh::fixtures::Gear(meta_data,
                          bulk_data,
                          meta_data.declare_part(oss.str(),static_cast<EntityRank>(SpatialDimension)),
                          cylindrical_coord_part,
                          hex_part,
                          wedge_part,
                          *cartesian_coord_field,
                          *displacement_field,
                          *translation_field,
                          *cylindrical_coord_field,
                          gear_params.element_size,
                          gear_params.radius_min,
                          gear_params.radius_max,
                          gear_params.height_min,
                          gear_params.height_max);
  }
}

GearsFixture::~GearsFixture()
{
  for (std::vector<stk::mesh::fixtures::Gear *>::iterator i = m_gears.begin(); i != m_gears.end(); ++i) {
    delete *i;
    *i = nullptr;
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
      stk::mesh::fixtures::Gear & gear = get_gear(i);
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
    stk::mesh::fixtures::distribute_gear_across_processors(get_gear(i), *cylindrical_coord_field);
  }

}

void GearsFixture::communicate_model_fields()
{
  //copy the field data to the aura nodes

  std::vector< const FieldBase *> fields;

  fields.push_back(cartesian_coord_field);
  fields.push_back(translation_field);
  fields.push_back(cylindrical_coord_field);
  fields.push_back(&displacement_field->field_of_state(stk::mesh::StateNew));
  fields.push_back(&displacement_field->field_of_state(stk::mesh::StateOld));

  // Parallel collective call:
#if defined( STK_HAS_MPI)
  communicate_field_data(bulk_data.aura_ghosting(), fields);
#endif
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


void select_nodal_data(const BulkData& mesh,
                       stk::mesh::fixtures::GearsFixture::CylindricalField & cylindrical_coord_field,
                       Entity element,
                       double & radius,
                       double & angle,
                       double & height)
{
  radius = 0.0;
  angle = TWO_PI;
  height = 0.0;

  int numNodes = mesh.num_nodes(element);
  Entity const *elem_nodes = mesh.begin_nodes(element);
  for (int i = 0; i < numNodes; ++i)
  {
    Entity node = elem_nodes[i];
    const MeshIndex& mi = mesh.mesh_index(node);
    const Bucket& bucket = *mi.bucket;
    double *cylindrical_data = stk::mesh::field_data(cylindrical_coord_field, bucket, mi.bucket_ordinal);
    radius += cylindrical_data[0];
    angle  = std::min(angle,cylindrical_data[1]);
    height += cylindrical_data[2];
  }
  radius /= numNodes;
  height /= numNodes;
}

// Parallel collective call:
void distribute_gear_across_processors(stk::mesh::fixtures::Gear & gear,
                                       stk::mesh::fixtures::GearsFixture::CylindricalField & cylindrical_coord_field)
{
  BulkData & bulk_data = gear.bulk_data;

  const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  EntityProcVec elements_to_change_owner;

  Selector locally_owned = gear.meta_data.locally_owned_part();
  if (p_rank == 0) {
    BucketVector const& all_elements = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_owned);
    std::set<Entity> node_set; // First come first serve nodal movement.
    for (BucketVector::const_iterator it = all_elements.begin() ; it != all_elements.end() ; ++it) {
      Bucket & b = **it;
      for (size_t i=0 ; i<b.size() ; ++i) {
        Entity element = b[i];
        double radius = 0.0;
        double angle  = 0.0;
        double height = 0.0;
        select_nodal_data(bulk_data, cylindrical_coord_field, element, radius, angle, height);
        unsigned destination_processor_rank = stk::mesh::fixtures::destination_processor(gear,radius,angle,height,p_rank,p_size);
        elements_to_change_owner.push_back(EntityProc(element,destination_processor_rank));

        // Now add all related nodes to list to move to this processor:
        Entity const *elem_nodes_j = bulk_data.begin_nodes(element);
        Entity const *elem_nodes_e = bulk_data.end_nodes(element);
        for ( ; elem_nodes_j != elem_nodes_e; ++elem_nodes_j)
        {
          Entity node = *elem_nodes_j;
          if (node_set.count(node)==0) {
            elements_to_change_owner.push_back(EntityProc(node,destination_processor_rank));
            node_set.insert(node);
          }
        }
      }
    }
  }

  // Parallel collective call:
  bulk_data.change_entity_owner(elements_to_change_owner);
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

unsigned destination_processor(const stk::mesh::fixtures::Gear & gear, double rad, double angle, double height, unsigned p_rank, unsigned p_size)
{
  unsigned result = 0;
  // Distribute elements across angles: (not working perfectly yet)
  angle = scale_angle_2pi(angle);
  result = static_cast<unsigned>(floor0((angle/TWO_PI)*p_size));

  scale_p_rank(result,p_size);

  return result;
}

} // namespace simple_fields

} // fixtures
} // mesh
} // stk
