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

#include <vector>
#include <stdexcept>
#include <sstream>

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

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

using namespace stk ;

namespace stk_examples {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef mesh::Field<double,mesh::Cartesian> VectorFieldType ;

//--------------------------------------------------------------------
//
// main driver for use-case 5: heterogeneous element mesh.
//

void use_case_5_generate_mesh_meta_data(
  stk::mesh::MetaData & meta_data ,
  VectorFieldType & node_coord )
{
  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;

  mesh::Part & universal        = meta_data.universal_part();
  mesh::Part & block_hex        = meta_data.declare_part("hexes",element_rank);
  mesh::Part & block_wedge      = meta_data.declare_part("wedges",element_rank);
  mesh::Part & block_tet        = meta_data.declare_part("tets",element_rank);
  mesh::Part & block_pyramid    = meta_data.declare_part("pyramids",element_rank);
  mesh::Part & block_quad_shell = meta_data.declare_part("quad_shells",element_rank);
  mesh::Part & block_tri_shell  = meta_data.declare_part("tri_shells",element_rank);

  stk::io::put_io_part_attribute(block_hex);
  stk::io::put_io_part_attribute(block_wedge);
  stk::io::put_io_part_attribute(block_tet);
  stk::io::put_io_part_attribute(block_pyramid);
  stk::io::put_io_part_attribute(block_quad_shell);
  stk::io::put_io_part_attribute(block_tri_shell);

  stk::mesh::set_cell_topology(block_hex            , stk::mesh::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8>         >()));
  stk::mesh::set_cell_topology(block_wedge          , stk::mesh::CellTopology(shards::getCellTopologyData<shards::Wedge<6>              >()));
  stk::mesh::set_cell_topology(block_tet            , stk::mesh::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4>        >()));
  stk::mesh::set_cell_topology(block_pyramid        , stk::mesh::CellTopology(shards::getCellTopologyData<shards::Pyramid<5>            >()));
  stk::mesh::set_cell_topology(block_quad_shell     , stk::mesh::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<4> >()));
  stk::mesh::set_cell_topology(block_tri_shell      , stk::mesh::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<3>      >()));

  const mesh::FieldBase::Restriction & res =
    stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

  if ( res.num_scalars_per_entity() != 3 ) {
    std::ostringstream msg ;
    msg << "stk_examples::use_case_5_generate_mesh_meta_data FAILED, coordinate dimension must be 3 != " << res.num_scalars_per_entity() ;
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

static const stk::mesh::EntityId hex_node_ids[3][ shards::Hexahedron<8> ::node_count ] = {
  { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
  { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
  { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId wedge_node_ids[3][ shards::Wedge<6> ::node_count ] = {
  { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
  { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
  { 16 , 17 , 20 ,  6 ,  7 , 10 } };

static const stk::mesh::EntityId tetra_node_ids[3][ shards::Tetrahedron<4> ::node_count ] = {
  { 15 , 19 , 16 , 21 } ,
  { 19 , 20 , 16 , 21 } ,
  { 16 , 20 , 17 , 21 } };

static const stk::mesh::EntityId pyramid_node_ids[2][ shards::Pyramid<5> ::node_count ] = {
  { 11 , 15 , 16 , 12 , 21 } ,
  { 12 , 16 , 17 , 13 , 21 } };

static const stk::mesh::EntityId shell_quad_node_ids[3][ shards::ShellQuadrilateral<4> ::node_count ]={
  { 9 , 6 , 16 , 19 } ,
  { 6 , 7 , 17 , 16 } ,
  { 7 , 8 , 18 , 17 } };

static const stk::mesh::EntityId shell_tri_node_ids[3][ shards::ShellTriangle<3> ::node_count ] ={
  { 19 , 16 , 21 } ,
  { 16 , 17 , 21 } ,
  { 17 , 13 , 21 } };

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void use_case_5_generate_mesh_bulk_data(
  mesh::BulkData & bulk_data ,
  const VectorFieldType & node_coord )
{
  static const char method[] =
    "stk_examples::use_case_5_generate_mesh_bulk_data" ;

  bulk_data.modification_begin();

  const mesh::MetaData & meta_data = stk::mesh::MetaData::get(bulk_data);

  mesh::Part & hex_block        = * meta_data.get_part("hexes",method);
  mesh::Part & wedge_block      = * meta_data.get_part("wedges",method);
  mesh::Part & tetra_block      = * meta_data.get_part("tets",method);
  mesh::Part & pyramid_block    = * meta_data.get_part("pyramids",method);
  mesh::Part & quad_shell_block = * meta_data.get_part("quad_shells",method);
  mesh::Part & tri_shell_block  = * meta_data.get_part("tri_shells",method);

  unsigned elem_id = 1 ;

  for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, hex_block, elem_id, hex_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_wedge ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, wedge_block, elem_id, wedge_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_tetra ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, tetra_block, elem_id, tetra_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, pyramid_block, elem_id, pyramid_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, quad_shell_block, elem_id, shell_quad_node_ids[i]);
  }

  for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++elem_id ) {
    mesh::declare_element( bulk_data, tri_shell_block, elem_id, shell_tri_node_ids[i] );
  }

  for ( unsigned i = 0 ; i < node_count ; ++i ) {

    mesh::Entity const node = bulk_data.get_entity( stk::topology::NODE_RANK , i + 1 );

    double * const coord = stk::mesh::field_data( node_coord , node );

    coord[0] = node_coord_data[i][0] ;
    coord[1] = node_coord_data[i][1] ;
    coord[2] = node_coord_data[i][2] ;
  }

  bulk_data.modification_end();
}


void use_case_5_write_mesh( stk::ParallelMachine comm ,
                            const std::string & filename )
{
  Ioss::Init::Initializer init_db;
  stk::mesh::MetaData meta_data(SpatialDim);

  VectorFieldType & node_coord =
    meta_data.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

  mesh::put_field( node_coord , meta_data.universal_part() , SpatialDim );

  use_case_5_generate_mesh_meta_data( meta_data , node_coord );

  meta_data.commit();

  mesh::BulkData bulk_data( meta_data, comm );

  use_case_5_generate_mesh_bulk_data( bulk_data , node_coord );

  stk::io::StkMeshIoBroker stkIoMeshBroker;
  stkIoMeshBroker.set_bulk_data(bulk_data);
  size_t resultFileIndex = stkIoMeshBroker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  stkIoMeshBroker.write_output_mesh(resultFileIndex);
}

} // namespace stk_examples

//----------------------------------------------------------------------

