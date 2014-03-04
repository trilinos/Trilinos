/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <vector>
#include <stdexcept>
#include <sstream>

#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

//--------------------------------------------------------------------
//

void heterogeneous_mesh_meta_data(
  stk::mesh::MetaData & meta_data ,
  VectorFieldType & node_coord )
{
  stk::mesh::Part & universal        = meta_data.universal_part();
  declare_part<shards::Hexahedron<8> >(meta_data, "hexes");
  declare_part<shards::Wedge<6> >(meta_data, "wedges");
  declare_part<shards::Tetrahedron<4> >(meta_data, "tets");
  declare_part<shards::Pyramid<5> >(meta_data, "pyramids");
  declare_part<shards::ShellQuadrilateral<4> >(meta_data, "quad_shells");
  declare_part<shards::ShellTriangle<3> >(meta_data, "tri_shells");
  
  const stk::mesh::FieldBase::Restriction & res =
    stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

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

static const double node_coord_data[ node_count ][ SpatialDim ] = {
  { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
  { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
  { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
  { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
  { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
  { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
  { 1 , 1 , -2 } };

static const stk::mesh::EntityId hex_node_ids[number_hex][ shards::Hexahedron<8> ::node_count ] = {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId wedge_node_ids[number_wedge][ shards::Wedge<6> ::node_count ] = {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 } };

static const stk::mesh::EntityId tetra_node_ids[number_tetra][ shards::Tetrahedron<4> ::node_count ] = {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 } };

static const stk::mesh::EntityId pyramid_node_ids[number_pyramid][ shards::Pyramid<5> ::node_count ] = {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 } };

static const stk::mesh::EntityId shell_quad_node_ids[number_shell_quad][ shards::ShellQuadrilateral<4> ::node_count ]={
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId shell_tri_node_ids[number_shell_tri][ shards::ShellTriangle<3> ::node_count ] ={
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 } };

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void heterogeneous_mesh_bulk_data(
  stk::mesh::BulkData & bulk_data ,
  const VectorFieldType & node_coord )
{
  static const char method[] =
    "stk_mesh::fixtures::heterogenous_mesh_bulk_data" ;

  bulk_data.modification_begin();

  const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(bulk_data);

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

}}}

//----------------------------------------------------------------------

