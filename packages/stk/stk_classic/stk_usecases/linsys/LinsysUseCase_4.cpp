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
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SetOwners.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/LinsysFunctions.hpp>

#include <fei_Factory_Trilinos.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 1: DOF mapping use case.
// The function 'use_case_4_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk_linsys_usecases {

enum { SpatialDim = 3 };
enum { MaximumEntityRank = 6 };

//----------------------------------------------------------------------

typedef stk::mesh::Field<double,stk::mesh::Cartesian>    VectorFieldType ;
typedef stk::mesh::Field<double>              ScalarFieldType ;

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_4_generate_mesh(
  const std::string& mesh_options,
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  stk::mesh::Part & hex_block ,
  stk::mesh::Part & quad_shell_block );

void use_case_4_initialize_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const VectorFieldType & node_displ ,
  const VectorFieldType & node_rotat );

//--------------------------------------------------------------------
//
// main driver for use-case 1: DOF mapping.
//

bool use_case_4_driver( MPI_Comm comm ,
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

//    ScalarFieldType & pressure_field =
      stk::mesh::put_field(
        fem_meta.declare_field< ScalarFieldType >("pressure"),
        element_rank, universal);

    //--------------------------------
    // rotation_field only exists on the shell-nodes:

    VectorFieldType & rotation_field =
      stk::mesh::put_field(
        fem_meta.declare_field< VectorFieldType >( "rotation" ),
        stk::mesh::fem::FEMMetaData::NODE_RANK , block_quad_shell , SpatialDim );

    stk::mesh::Part& bcpart = fem_meta.declare_part("bcpart");

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

    use_case_4_generate_mesh(
      mesh_options ,
      mesh_bulk_data ,
      coordinates_field ,
      block_hex ,
      block_quad_shell );

    use_case_4_initialize_data(
      mesh_bulk_data ,
      coordinates_field ,
      displacements_field ,
      rotation_field );

    //Add a node to our boundary-condition part 'bcpart'.
    //let's choose the first locally-owned node. (This will produce a
    //different boundary-condition for different numbers of processors...
    //A more realistic case would simply pick a specific set of nodes
    //regardless of which processors they are on.)

    mesh_bulk_data.modification_begin();

    std::vector<stk::mesh::Entity*> local_nodes;
    stk::mesh::Selector select_owned(fem_meta.locally_owned_part());
    stk::mesh::get_selected_entities(select_owned,
                                     mesh_bulk_data.buckets(stk::mesh::fem::FEMMetaData::NODE_RANK),
                                     local_nodes);

    if (local_nodes.size() > 0) {
      stk::mesh::PartVector partvector;
      partvector.push_back(&bcpart);
      mesh_bulk_data.change_entity_parts(*local_nodes[0], partvector);
    }

    mesh_bulk_data.modification_end();

    //set owner-processors to lowest-sharing (stk::mesh defaults to
    //highest-sharing) If highest-sharing owns, then it isn't correct for the
    //way the fei library sets ownership of shared nodes for vectors etc.
    stk::mesh::set_owners<stk::mesh::LowestRankSharingProcOwns>( mesh_bulk_data );

    //------------------------------------------------------------------

    const unsigned myProc = mesh_bulk_data.parallel_rank();

    fei::SharedPtr<fei::Factory> feifactory(new Factory_Trilinos(comm));
    stk::linsys::LinearSystem ls(comm, feifactory);

    if (myProc == 0) {
      std::cout << "Adding element-node connectivities for displacements field for all locally-owned "
        << "elements..." << std::endl;
    }

    stk::linsys::add_connectivities(ls, element_rank, stk::mesh::fem::FEMMetaData::NODE_RANK,
                                    displacements_field, select_owned, mesh_bulk_data);

    ls.synchronize_mappings_and_structure();

    ls.create_fei_LinearSystem();

    fei::SharedPtr<fei::MatrixGraph> matgraph = ls.get_fei_MatrixGraph();
    fei::SharedPtr<fei::Matrix> matrix = ls.get_fei_LinearSystem()->getMatrix();
    fei::SharedPtr<fei::Vector> rhs = ls.get_fei_LinearSystem()->getRHS();

    {
      const std::vector<stk::mesh::Bucket*>& mesh_buckets = mesh_bulk_data.buckets(element_rank);
      std::vector<stk::mesh::Bucket*> part_buckets;
      stk::mesh::get_buckets(select_owned, mesh_buckets, part_buckets);

      stk::linsys::DofMapper& dof_mapper = ls.get_DofMapper();

      int field_id = dof_mapper.get_field_id(displacements_field);

      stk::mesh::Entity& first_entity = *(part_buckets[0]->begin());
      stk::mesh::PairIterRelation rel = first_entity.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
      int num_nodes_per_elem = rel.second - rel.first;

      int pattern_id = matgraph->definePattern(num_nodes_per_elem, stk::mesh::fem::FEMMetaData::NODE_RANK, field_id);

      std::vector<int> node_ids(num_nodes_per_elem);

      const int field_size = dof_mapper.get_fei_VectorSpace()->getFieldSize(field_id);
      const int matsize = num_nodes_per_elem*field_size*num_nodes_per_elem*field_size;
      const int vecsize = num_nodes_per_elem*field_size;

      std::vector<double> elem_matrix_1d(matsize, 0);
      std::vector<double*> elem_matrix_2d(vecsize);

      std::vector<double> elem_vector(vecsize, 0);

      for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
        elem_matrix_2d[i] = &elem_matrix_1d[i*vecsize];
      }

      //fill our dummy elem-matrix:
      for(size_t i=0; i<elem_matrix_2d.size(); ++i) {
        double* row = elem_matrix_2d[i];
        if (i>=1) row[i-1] = -1;
        row[i] = 2;
        if (i<elem_matrix_2d.size()-1) row[i+1] = -1;

        elem_vector[i] = 1;
      }

      std::vector<int> eqn_indices(vecsize);

      for(size_t i=0; i<part_buckets.size(); ++i) {
        stk::mesh::Bucket::iterator
          b_iter = part_buckets[i]->begin(),
                 b_end  = part_buckets[i]->end();
        for(; b_iter != b_end; ++b_iter) {
          stk::mesh::Entity& elem = *b_iter;
          rel = elem.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
          for(int j=0; rel.first != rel.second; ++rel.first, ++j) {
            node_ids[j] = rel.first->entity()->identifier();
          }

          matgraph->getPatternIndices(pattern_id, &node_ids[0], eqn_indices);

          matrix->sumIn(vecsize, &eqn_indices[0], vecsize, &eqn_indices[0],
                        &elem_matrix_2d[0]);
          rhs->sumIn(vecsize, &eqn_indices[0], &elem_vector[0]);
        }
      }

    }

    stk::linsys::dirichlet_bc(ls, mesh_bulk_data, bcpart, stk::mesh::fem::FEMMetaData::NODE_RANK,
                              displacements_field, 0, 3.14159265);

    ls.finalize_assembly();

    matrix->writeToFile("A.mtx");
    rhs->writeToFile("rhs.vec");
  }
  return true;
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

void use_case_4_initialize_data(
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

void use_case_4_generate_mesh(
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

    {

      for ( size_t i = 1 ; i <= num_block ; ++i ) {
        const size_t                        num_elem = gmesh.element_count_proc(i);
        const std::pair<std::string,int> top_info = gmesh.topology_type(i);

	std::vector<int> elem_map( num_elem , 0 );
        std::vector<int> elem_conn( num_elem * top_info.second );

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

