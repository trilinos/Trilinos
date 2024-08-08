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

#include <sstream>                      // for ostringstream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase::Restriction, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/stk_mesh_fixtures/heterogeneous_mesh.hpp>
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

//--------------------------------------------------------------------
//

void heterogeneous_mesh_meta_data(stk::mesh::MetaData & meta_data, const VectorFieldType & node_coord)
{
  stk::mesh::Part & universal = meta_data.universal_part();
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("hexes", stk::topology::HEX_8));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("wedges", stk::topology::WEDGE_6));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("tets", stk::topology::TET_4));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("pyramids", stk::topology::PYRAMID_5));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("quad_shells", stk::topology::SHELL_QUAD_4));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("tri_shells", stk::topology::SHELL_TRI_3));

  const stk::mesh::FieldBase::Restriction & res = stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

  if ( res.num_scalars_per_entity() != 3 ) {
    std::ostringstream msg ;
    msg << "stk_mesh/unit_tests/heterogenous_mesh_meta_data FAILED, coordinate dimension must be 3 != " << res.num_scalars_per_entity() ;
    throw std::runtime_error( msg.str() );
  }
}

//--------------------------------------------------------------------
/*----------------------------------------------------------------------
 * Internal use-case #5 mesh generation.
 *
 * Three hexes, three wedges, three tets, two pyramids,
 * three quad shells, and three triangle shells.
 *
 *  Z = 0 plane:
 *
 *    Y
 *    ^   9      10
 *    !   *-------*
 *    !  / \     / \
 *    ! /   \   /   \
 *     /     \ /     \
 *    *-------*-------*-------*
 *   5|      6|      7|      8|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *    1       2       3       4
 *
 *  Z = -1 plane:
 *
 *    Y
 *    ^  19      20
 *    !   *-------*
 *    !  / \     / \
 *    ! /   \   /   \
 *     /     \ /     \
 *    *-------*-------*-------*
 *  15|     16|     17|     18|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *   11      12      13      14
 *
 *
 *  Last node (#21) at Z = -2, translated from node #16
 *----------------------------------------------------------------------*/

enum { node_count = 21 };

enum { number_hex = 3 };
enum { number_wedge = 3 };
enum { number_tetra = 3 };
enum { number_pyramid = 2 };
enum { number_shell_quad = 3 };
enum { number_shell_tri = 3 };

namespace {

typedef stk::topology::topology_type<stk::topology::HEX_8> Hex8;
typedef stk::topology::topology_type<stk::topology::WEDGE_6> Wedge6;
typedef stk::topology::topology_type<stk::topology::TET_4> Tet4;
typedef stk::topology::topology_type<stk::topology::PYRAMID_5> Pyramid5;
typedef stk::topology::topology_type<stk::topology::SHELL_QUAD_4> ShellQuad4;
typedef stk::topology::topology_type<stk::topology::SHELL_TRI_3> ShellTri3;

static const double node_coord_data[ node_count ][ SpatialDim ] = {
  { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
  { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
  { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
  { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
  { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
  { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
  { 1 , 1 , -2 }
};

static const stk::mesh::EntityIdVector hex_node_ids[number_hex] {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 }
};

static const stk::mesh::EntityIdVector wedge_node_ids[number_wedge] {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 20,  19,  16,  10 ,  9 ,  6} ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 }
};

static const stk::mesh::EntityIdVector tetra_node_ids[number_tetra] {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 }
};

static const stk::mesh::EntityIdVector pyramid_node_ids[number_pyramid] {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 }
};

static const stk::mesh::EntityIdVector shell_quad_node_ids[number_shell_quad] {
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 }
};

static const stk::mesh::EntityIdVector shell_tri_node_ids[number_shell_tri] {
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 }
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void heterogeneous_mesh_bulk_data(stk::mesh::BulkData & bulk_data, const VectorFieldType & node_coord )
{
  static const char method[] = "stk_mesh::fixtures::heterogenous_mesh_bulk_data" ;

  bulk_data.modification_begin();

  const stk::mesh::MetaData & meta_data = bulk_data.mesh_meta_data();

  stk::mesh::Part & hex_block        = * meta_data.get_part("hexes",method);
  stk::mesh::Part & wedge_block      = * meta_data.get_part("wedges",method);
  stk::mesh::Part & tetra_block      = * meta_data.get_part("tets",method);
  stk::mesh::Part & pyramid_block    = * meta_data.get_part("pyramids",method);
  stk::mesh::Part & quad_shell_block = * meta_data.get_part("quad_shells",method);
  stk::mesh::Part & tri_shell_block  = * meta_data.get_part("tri_shells",method);

  unsigned elem_id = 1 ;

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, hex_block, elem_id, hex_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, wedge_block, elem_id, wedge_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_tetra ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, tetra_block, elem_id, tetra_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, pyramid_block, elem_id, pyramid_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, quad_shell_block, elem_id, shell_quad_node_ids[i]);
  }

  for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, tri_shell_block, elem_id, shell_tri_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {

    stk::mesh::Entity const node = bulk_data.get_entity( stk::topology::NODE_RANK , i + 1 );

    double * const coord = stk::mesh::field_data( node_coord , node );

    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  bulk_data.modification_end();
}

namespace simple_fields {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

//--------------------------------------------------------------------
//

void heterogeneous_mesh_meta_data(stk::mesh::MetaData & meta_data, const VectorFieldType & node_coord)
{
  stk::mesh::Part & universal = meta_data.universal_part();
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("hexes", stk::topology::HEX_8));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("wedges", stk::topology::WEDGE_6));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("tets", stk::topology::TET_4));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("pyramids", stk::topology::PYRAMID_5));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("quad_shells", stk::topology::SHELL_QUAD_4));
  stk::io::put_io_part_attribute(meta_data.declare_part_with_topology("tri_shells", stk::topology::SHELL_TRI_3));

  const stk::mesh::FieldBase::Restriction & res = stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

  if ( res.num_scalars_per_entity() != 3 ) {
    std::ostringstream msg ;
    msg << "stk_mesh/unit_tests/heterogenous_mesh_meta_data FAILED, coordinate dimension must be 3 != " << res.num_scalars_per_entity() ;
    throw std::runtime_error( msg.str() );
  }
}

//--------------------------------------------------------------------
/*----------------------------------------------------------------------
 * Internal use-case #5 mesh generation.
 *
 * Three hexes, three wedges, three tets, two pyramids,
 * three quad shells, and three triangle shells.
 *
 *  Z = 0 plane:
 *
 *    Y
 *    ^   9      10
 *    !   *-------*
 *    !  / \     / \
 *    ! /   \   /   \
 *     /     \ /     \
 *    *-------*-------*-------*
 *   5|      6|      7|      8|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *    1       2       3       4
 *
 *  Z = -1 plane:
 *
 *    Y
 *    ^  19      20
 *    !   *-------*
 *    !  / \     / \
 *    ! /   \   /   \
 *     /     \ /     \
 *    *-------*-------*-------*
 *  15|     16|     17|     18|
 *    |       |       |       |
 *    |       |       |       |
 *    *-------*-------*-------*    ----> X
 *   11      12      13      14
 *
 *
 *  Last node (#21) at Z = -2, translated from node #16
 *----------------------------------------------------------------------*/

enum { node_count = 21 };

enum { number_hex = 3 };
enum { number_wedge = 3 };
enum { number_tetra = 3 };
enum { number_pyramid = 2 };
enum { number_shell_quad = 3 };
enum { number_shell_tri = 3 };

namespace {

typedef stk::topology::topology_type<stk::topology::HEX_8> Hex8;
typedef stk::topology::topology_type<stk::topology::WEDGE_6> Wedge6;
typedef stk::topology::topology_type<stk::topology::TET_4> Tet4;
typedef stk::topology::topology_type<stk::topology::PYRAMID_5> Pyramid5;
typedef stk::topology::topology_type<stk::topology::SHELL_QUAD_4> ShellQuad4;
typedef stk::topology::topology_type<stk::topology::SHELL_TRI_3> ShellTri3;

static const double node_coord_data[ node_count ][ SpatialDim ] = {
  { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
  { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
  { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
  { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
  { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
  { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
  { 1 , 1 , -2 }
};

static const stk::mesh::EntityIdVector hex_node_ids[number_hex] {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 }
};

static const stk::mesh::EntityIdVector wedge_node_ids[number_wedge] {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 20,  19,  16,  10 ,  9 ,  6} ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 }
};

static const stk::mesh::EntityIdVector tetra_node_ids[number_tetra] {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 }
};

static const stk::mesh::EntityIdVector pyramid_node_ids[number_pyramid] {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 }
};

static const stk::mesh::EntityIdVector shell_quad_node_ids[number_shell_quad] {
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 }
};

static const stk::mesh::EntityIdVector shell_tri_node_ids[number_shell_tri] {
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 }
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void heterogeneous_mesh_bulk_data(stk::mesh::BulkData & bulk_data, const VectorFieldType & node_coord )
{
  static const char method[] = "stk_mesh::fixtures::heterogenous_mesh_bulk_data" ;

  bulk_data.modification_begin();

  const stk::mesh::MetaData & meta_data = bulk_data.mesh_meta_data();

  stk::mesh::Part & hex_block        = * meta_data.get_part("hexes",method);
  stk::mesh::Part & wedge_block      = * meta_data.get_part("wedges",method);
  stk::mesh::Part & tetra_block      = * meta_data.get_part("tets",method);
  stk::mesh::Part & pyramid_block    = * meta_data.get_part("pyramids",method);
  stk::mesh::Part & quad_shell_block = * meta_data.get_part("quad_shells",method);
  stk::mesh::Part & tri_shell_block  = * meta_data.get_part("tri_shells",method);

  unsigned elem_id = 1 ;

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, hex_block, elem_id, hex_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, wedge_block, elem_id, wedge_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_tetra ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, tetra_block, elem_id, tetra_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, pyramid_block, elem_id, pyramid_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, quad_shell_block, elem_id, shell_quad_node_ids[i]);
  }

  for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++elem_id ) {
    stk::mesh::declare_element( bulk_data, tri_shell_block, elem_id, shell_tri_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {

    stk::mesh::Entity const node = bulk_data.get_entity( stk::topology::NODE_RANK , i + 1 );

    double * const coord = stk::mesh::field_data( node_coord , node );

    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  bulk_data.modification_end();
}

} // namespace simple_fields

}}}

//----------------------------------------------------------------------

