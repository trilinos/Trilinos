/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_linsys/DofMapper.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 1: DOF mapping use case.
// The function 'use_case_2_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk_linsys_usecases {

enum { SpatialDim = 3 };
enum { MaximumEntityRank = 6 };

//----------------------------------------------------------------------

typedef stk::mesh::Field<double,stk::mesh::Cartesian>    VectorFieldType ;
typedef stk::mesh::Field<double>              ScalarFieldType ;

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_2_generate_mesh(
  const std::string& mesh_options,
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  stk::mesh::Part & hex_block ,
  stk::mesh::Part & quad_shell_block );

void use_case_2_initialize_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const VectorFieldType & node_displ ,
  const VectorFieldType & node_rotat );

//--------------------------------------------------------------------
//
// main driver for use-case 1: DOF mapping.
//

bool use_case_2_driver( MPI_Comm comm ,
                        const std::string& mesh_options )
{
  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stk_linsys use case 1" << std::endl
              << "  Number Processes = " << stk::parallel_machine_size( comm )
              << std::endl ;
  }

  //--------------------------------------------------------------------

  {
    //------------------------------------------------------------------
    // Declare the mesh meta data: element blocks and associated fields

    stk::mesh::fem::FEMMetaData fem_meta;
    fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim) ) ;

    stk::mesh::MetaData & mesh_meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
    const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

    //--------------------------------
    // Element-block declarations typically occur when reading the
    // mesh-file meta-data, and thus won't usually appear in application code.
    // Declaring the element blocks and associating an element traits
    // with each element block.

    stk::mesh::Part & universal        = fem_meta.universal_part();
    stk::mesh::Part & block_hex        = fem_meta.declare_part("block_1", element_rank);
    stk::mesh::Part & block_quad_shell = fem_meta.declare_part("block_2", element_rank);

    stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
    stk::mesh::fem::CellTopology qshell_top(shards::getCellTopologyData<shards::ShellQuadrilateral<> >());
    stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    stk::mesh::fem::set_cell_topology( block_quad_shell, qshell_top );

    //--------------------------------
    // Declaring fields of specified types on all nodes:

    VectorFieldType & coordinates_field =
      stk::mesh::put_field(
        fem_meta.declare_field< VectorFieldType >( "coordinates" ) ,
        stk::mesh::fem::FEMMetaData::NODE_RANK , universal , SpatialDim );

    VectorFieldType & displacements_field =
      stk::mesh::put_field(
        fem_meta.declare_field< VectorFieldType >( "displacements" ) ,
        stk::mesh::fem::FEMMetaData::NODE_RANK , universal , SpatialDim );

    //--------------------------------
    // Put a scalar "pressure" field on all elements, just to use in demonstrating
    // DOF mappings below:

    ScalarFieldType & pressure_field =
      stk::mesh::put_field(
        fem_meta.declare_field< ScalarFieldType >("pressure"),
        element_rank, universal);

    //--------------------------------
    // rotation_field only exists on the shell-nodes:

    VectorFieldType & rotation_field =
      stk::mesh::put_field(
        fem_meta.declare_field< VectorFieldType >( "rotation" ),
        stk::mesh::fem::FEMMetaData::NODE_RANK , block_quad_shell , SpatialDim );

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    fem_meta.commit();

    //------------------------------------------------------------------
    // stk::mesh::BulkData bulk data conforming to the meta data.

    stk::mesh::BulkData mesh_bulk_data( mesh_meta_data , comm );

    // In a typical app, the mesh would be read from file at this point.
    // But in this use-case, we generate the mesh and initialize
    // field data to use-case defined values.

    use_case_2_generate_mesh(
      mesh_options ,
      mesh_bulk_data ,
      coordinates_field ,
      block_hex ,
      block_quad_shell );

    use_case_2_initialize_data(
      mesh_bulk_data ,
      coordinates_field ,
      displacements_field ,
      rotation_field );

    mesh_bulk_data.modification_end();

    //------------------------------------------------------------------

    const unsigned myProc = mesh_bulk_data.parallel_rank();

    stk::linsys::DofMapper dof_mapper(comm);

    if (myProc == 0) {
      std::cout << "Adding DOF mappings for displacements field for all locally-used "
        << "(owned and shared) nodes..." << std::endl;
    }

    stk::mesh::Selector owned_and_shared =
      fem_meta.locally_owned_part() |
      fem_meta.globally_shared_part();

    dof_mapper.add_dof_mappings(mesh_bulk_data, owned_and_shared,
                                stk::mesh::fem::FEMMetaData::NODE_RANK, displacements_field);

    if (myProc == 0) {
      std::cout << "Adding DOF mappings for pressure field for all locally-owned "
        << " elements..." << std::endl;
    }

    stk::mesh::Selector select_owned = fem_meta.locally_owned_part();
    dof_mapper.add_dof_mappings(mesh_bulk_data, select_owned,
                                element_rank, pressure_field);

    dof_mapper.finalize();

    if (myProc == 0) {
      std::cout << "Global Number of Indices: "
          << dof_mapper.get_fei_VectorSpace()->getGlobalNumIndices() << std::endl;
    }

    std::vector<stk::mesh::Entity*> nodes;
    stk::mesh::get_entities(mesh_bulk_data, stk::mesh::fem::FEMMetaData::NODE_RANK, nodes);

    std::vector<stk::mesh::Entity*> elems;
    stk::mesh::get_entities(mesh_bulk_data, element_rank, elems);

    //The object of this use-case is to perform reverse-lookups. i.e., given a
    //global_index, query to obtain entity type/id info.
    //To verify correctness we will first query to obtain the global_index, then do
    //the reverse-lookup and make sure we get back the right values.

    int global_index = 0;

    for(size_t i=0; i<nodes.size(); i+=1000) {
      //is the i-th node in the locally-used part? If not, continue.
      if (! owned_and_shared( nodes[i]->bucket())) continue;

      global_index = dof_mapper.get_global_index(stk::mesh::fem::FEMMetaData::NODE_RANK, nodes[i]->identifier(), displacements_field);
      std::cout << "Proc " << myProc << ", global index for node " << nodes[i]->identifier()
        << ", field '"<<displacements_field.name()<<"' is: " << global_index << std::endl;
      stk::mesh::EntityRank ent_type;
      stk::mesh::EntityId ent_id;
      const stk::mesh::FieldBase* field = NULL;
      int offset_into_field;
      dof_mapper.get_dof(global_index, ent_type, ent_id, field, offset_into_field);
      if (ent_type != stk::mesh::fem::FEMMetaData::NODE_RANK || ent_id != nodes[i]->identifier() ||
          field->name() != displacements_field.name()) {
        std::cout << "Reverse-lookup Test Failed" << std::endl;
      }
    }

    for(size_t i=0; i<elems.size(); i+=1000) {
      //is the i-th elem in the locally-owned part? If not, continue.
      if (!elems[i]->bucket().member(fem_meta.locally_owned_part())) continue;

      global_index = dof_mapper.get_global_index(element_rank, elems[i]->identifier(), pressure_field);
      std::cout << "Proc " << myProc << ", global index for element " << elems[i]->identifier()
        << ", field '"<<pressure_field.name()<<"' is: " << global_index << std::endl;
      stk::mesh::EntityRank ent_type;
      stk::mesh::EntityId ent_id;
      const stk::mesh::FieldBase* field = NULL;
      int offset_into_field;
      dof_mapper.get_dof(global_index, ent_type, ent_id, field, offset_into_field);
      if (ent_type != element_rank || ent_id != elems[i]->identifier() ||
          field->name() != pressure_field.name()) {
        std::cout << "Reverse-lookup Test Failed" << std::endl;
      }
    }

    if (mesh_options == "10x10x10+shell:y") {
      return dof_mapper.get_fei_VectorSpace()->getGlobalNumIndices() == 5093;
    }
    return true;
  }
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

void use_case_2_initialize_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const VectorFieldType & node_displ ,
  const VectorFieldType & node_rotat )
{
  const std::vector<stk::mesh::Bucket*> & buckets = mesh.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
    stk::mesh::Bucket & bucket = **k ;
    const unsigned length = bucket.size();
    const unsigned length_3 = length * 3 ;

    double * const coord = stk::mesh::field_data( node_coord , bucket.begin() );
    double * const displ = stk::mesh::field_data( node_displ , bucket.begin() );
    double * const rotat = stk::mesh::field_data( node_rotat , bucket.begin() );

    for ( unsigned i = 0 ; i < length_3 ; ++i ) {
      displ[i] = 0.1 * coord[i] ;
    }

    if ( rotat ) {
      for ( unsigned i = 0 ; i < length_3 ; ++i ) {
        rotat[i] = 0.01 * coord[i] ;
      }
    }
  }
}

} // namespace stk_linsys_usecases

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <generated/Iogn_GeneratedMesh.h>

namespace stk_linsys_usecases {

void use_case_2_generate_mesh(
  const std::string& mesh_options ,
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  stk::mesh::Part & hex_block ,
  stk::mesh::Part & quad_shell_block )
{
  mesh.modification_begin();

  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();

  double t = 0 ;
  size_t num_hex = 0 ;
  size_t num_shell = 0 ;
  size_t num_nodes = 0 ;
  size_t num_block = 0 ;
  int error_flag = 0 ;

  try {

    Iogn::GeneratedMesh gmesh( mesh_options, parallel_size, parallel_rank );

    num_nodes = gmesh.node_count_proc();
    num_block = gmesh.block_count();

    t = stk::wall_time();

    std::vector<int> node_map( num_nodes , 0 );

    gmesh.node_map( node_map );

    for ( size_t i = 1 ; i <= num_block ; ++i ) {
      const size_t                     num_elem = gmesh.element_count_proc(i);
      const std::pair<std::string,int> top_info = gmesh.topology_type(i);

      std::vector<int> elem_map(  num_elem, 0 );
      std::vector<int> elem_conn( num_elem * top_info.second, 0 );

      gmesh.element_map( i, elem_map );
      gmesh.connectivity( i , elem_conn );

      if ( top_info.second == 8 ) {

        for ( size_t j = 0 ; j < num_elem ; ++j ) {

          const int * const local_node_id = & elem_conn[ j * 8 ] ;

          const stk::mesh::EntityId node_id[8] = {
            local_node_id[0] ,
            local_node_id[1] ,
            local_node_id[2] ,
            local_node_id[3] ,
            local_node_id[4] ,
            local_node_id[5] ,
            local_node_id[6] ,
            local_node_id[7]
          };

          const stk::mesh::EntityId elem_id = elem_map[ j ];

          stk::mesh::fem::declare_element( mesh , hex_block , elem_id , node_id );

          ++num_hex ;
        }
      }
      else if ( top_info.second == 4 ) {

        for ( size_t j = 0 ; j < num_elem ; ++j ) {

          const int * const local_node_id = & elem_conn[ j * 4 ] ;

          const stk::mesh::EntityId node_id[4] = {
            local_node_id[0] ,
            local_node_id[1] ,
            local_node_id[2] ,
            local_node_id[3]
          };

          const stk::mesh::EntityId elem_id = elem_map[ j ];

          stk::mesh::fem::declare_element( mesh , quad_shell_block , elem_id , node_id );

          ++num_shell ;
        }
      }
    }

    std::vector<double> node_coordinates( 3 * node_map.size() );

    gmesh.coordinates( node_coordinates );

    if ( 3 * node_map.size() != node_coordinates.size() ) {
      std::ostringstream msg ;
      msg << "  P" << mesh.parallel_rank()
          << ": ERROR, node_map.size() = "
          << node_map.size()
          << " , node_coordinates.size() / 3 = "
          << ( node_coordinates.size() / 3 );
      throw std::runtime_error( msg.str() );
    }

    for ( unsigned i = 0 ; i < node_map.size() ; ++i ) {
      const unsigned i3 = i * 3 ;

      stk::mesh::Entity * const node = mesh.get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK , node_map[i] );

      if ( NULL == node ) {
        std::ostringstream msg ;
        msg << "  P:" << mesh.parallel_rank()
            << " ERROR, Node not found: "
            << node_map[i] << " = node_map[" << i << "]" ;
        throw std::runtime_error( msg.str() );
      }

      double * const data = field_data( node_coord , *node );
      data[0] = node_coordinates[ i3 + 0 ];
      data[1] = node_coordinates[ i3 + 1 ];
      data[2] = node_coordinates[ i3 + 2 ];
    }
  }
  catch ( const std::exception & X ) {
    std::cout << "  P:" << mesh.parallel_rank() << ": " << X.what()
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }
  catch( ... ) {
    std::cout << "  P:" << mesh.parallel_rank()
              << " Caught unknown exception"
              << std::endl ;
    std::cout.flush();
    error_flag = 1 ;
  }

  stk::all_reduce( mesh.parallel() , stk::ReduceMax<1>( & error_flag ) );

  if ( error_flag ) {
    std::string msg( "Failed mesh generation" );
    throw std::runtime_error( msg );
  }

  mesh.modification_end();

  double dt = stk::wall_dtime( t );

  stk::all_reduce( mesh.parallel() , stk::ReduceMax<1>( & dt ) );

  std::cout << "  P" << mesh.parallel_rank()
            << ": Meshed Hex = " << num_hex
            << " , Shell = " << num_shell
            << " , Node = " << num_nodes
            << " in " << dt << " sec"
            << std::endl ;
  std::cout.flush();
}

} // namespace stk_linsys_usecases

