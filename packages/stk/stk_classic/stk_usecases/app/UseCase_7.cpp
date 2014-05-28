/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

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

#include <stk_algsup/AlgorithmRunner.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>

#include <common/gnu_malloc_hooks.hpp>

#include <app/performance_algorithms.hpp>
#include <app/UseCase_blas_algs.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 7: performance case.
// The function 'use_case_7_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk {
namespace app {

enum { SpatialDim = 3 };

//----------------------------------------------------------------------

typedef stk::mesh::Field<double,stk::mesh::Cartesian>    VectorFieldType ;
typedef stk::mesh::Field<double>              ScalarFieldType ;

// Specification for the aggressive gather pointer-field for elements.

typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerFieldType ;

// Centroid algorithm generic programming functions:

#include <mesh/centroid_algorithm.hpp>

//--------------------------------
// prototype for the function that will generate the use-case mesh.

void use_case_7_generate_mesh(
  const std::string& mesh_options,
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & elem_node_coord ,
  stk::mesh::Part & hex_block ,
  stk::mesh::Part & quad_shell_block );

void use_case_7_initialize_data(
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const VectorFieldType & node_displ ,
  const VectorFieldType & node_rotat );

bool use_case_7_driver(
  MPI_Comm comm ,
  int num_threads ,
  const std::string& thread_runner ,
  const std::string & mesh_options ,
  const unsigned num_trials );

//--------------------------------------------------------------------
//
// main driver for use-case 7: heterogeneous element mesh.
//

bool use_case_7_driver( MPI_Comm comm ,
                        bool performance_test ,
                        const std::string& mesh_options ,
                        int num_threads ,
                        const std::string& thread_runner )
{
  if ( ! stk::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh Use Case #7, begin" << std::endl ;
  }

  int num_trials = 1;
  if ( performance_test ) num_trials = 100;

  return use_case_7_driver( comm , num_threads, thread_runner, mesh_options, num_trials );

//The following 3 cases are run by 3 different xml files in the
//stk_app_rtest/performance/use_case_7 regression test.
//
//    1000x1x1000 plate: 1M hex, 1M shells
//        shells on min-y face
//
//    50x50x500 beam: 1.25M hex, 2.5K shells
//        shells on min-z face
//
//    100x100x100 cube: 1M hex, 60K shells
//        shells on all 6 faces
//
}

//--------------------------------------------------------------------

bool use_case_7_driver(
  MPI_Comm comm ,
  int num_threads ,
  const std::string& thread_runner ,
  const std::string & mesh_options ,
  const unsigned num_trials )
{
  // stk::mesh::BulkData bulk data CHUNK_SIZE = max entities per field data chunk.
  // (CHUNK_SIZE is the application's preferred "workset" size.)

  // Timing:
  //   [0] = BulkData creation
  //   [1] = BulkData synchronization
  //   [2] = Element computation #1
  //   [3] = Element computation #1 gathered
  //   [4] = Element computation #2
  //   [5] = Element computation #2 gathered
  //   [6] = Node computation #1
  //   [7] = Node computation #2
  //   [8] = stk::mesh::BulkData output
  //   [9] = stk::mesh::BulkData destruction

  double time_min[10] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double time_max[10] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double wtime = 0 ;

  double result_inner_product[2] = { 0 , 0 };
  std::vector<unsigned> hex_block_count ;
  std::vector<unsigned> shell_block_count ;


  const AlgorithmRunnerInterface* alg_runner = NULL ;
  if (thread_runner.empty() || thread_runner == "NonThreaded") {
    alg_runner = stk::algorithm_runner_non_thread();
  }
  else if (thread_runner == "TPI") {
    alg_runner = stk::algorithm_runner_tpi(num_threads);
  }
  else if ( thread_runner == std::string("TBB") ) {
    alg_runner = stk::algorithm_runner_tbb(num_threads);
  }

  if (alg_runner != NULL) {
    if (stk::parallel_machine_rank(comm) == 0)
      std::cout << "Using " << thread_runner
                << " algorithm runner, num_threads = " << num_threads
                << std::endl;
  } else {
    std::cout << "ERROR, failed to obtain requested AlgorithmRunner '"
              << thread_runner << "'." << std::endl;
    return false;
  }


  //--------------------------------------------------------------------

  reset_malloc_stats();

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh use case 7" << std::endl
              << "  Number Processes = " << stk::parallel_machine_size( comm )
              << std::endl ;
    if ( num_threads ) {
      std::cout << "  Threads/process  = " << num_threads << std::endl ;
    }
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  wtime = stk::wall_time();

  //------------------------------------------------------------------
  // Declare the mesh meta data: element blocks and associated fields
  stk::mesh::fem::FEMMetaData fem_meta;
  fem_meta.FEM_initialize(SpatialDim, stk::mesh::fem::entity_rank_names(SpatialDim) ) ;
  stk::mesh::MetaData & mesh_meta_data = stk::mesh::fem::FEMMetaData::get_meta_data(fem_meta);
  const stk::mesh::EntityRank element_rank = fem_meta.element_rank();

  {

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
      stk::mesh::put_field(fem_meta.declare_field< VectorFieldType >( "coordinates" ) , fem_meta.node_rank() , universal , SpatialDim );

    VectorFieldType & displacements_field =
      stk::mesh::put_field(fem_meta.declare_field< VectorFieldType >( "displacements" ) , fem_meta.node_rank() , universal , SpatialDim );

    VectorFieldType & residual_field =
      stk::mesh::put_field(fem_meta.declare_field< VectorFieldType >( "residual" ) ,fem_meta.node_rank() , universal , SpatialDim );

    //--------------------------------
    // rotation_field only exists on the shell-nodes:

    VectorFieldType & rotation_field =
      stk::mesh::put_field(fem_meta.declare_field< VectorFieldType >( "rotation" ), stk::mesh::fem::FEMMetaData::NODE_RANK , block_quad_shell , SpatialDim );

    //--------------------------------
    // Declare the coordinates field for all elements - for the centroid

    stk::mesh::put_field(coordinates_field , element_rank , universal , SpatialDim );

    // Declare rotation field on elements also - for the mean rotation.

    stk::mesh::put_field(rotation_field , element_rank , block_quad_shell , SpatialDim );

    //--------------------------------
    // Declare an aggressive "gather" field which is an
    // array of pointers to the element's nodes' coordinate field data.
    // The declaration specifies:
    //
    //     double * elem_node_coord[number_of_nodes]

    ElementNodePointerFieldType & elem_node_coord =
      fem_meta.declare_field< ElementNodePointerFieldType >( "elem_node_coord" );

    ElementNodePointerFieldType & elem_node_rot =
      fem_meta.declare_field< ElementNodePointerFieldType >( "elem_node_rot" );

    // Declare that the 'elem_node_coord' pointer field data
    // points to the 'coordinates_field' data on the nodes.

    fem_meta.declare_field_relation(
      elem_node_coord ,
      stk::mesh::fem::get_element_node_stencil(SpatialDim) ,
      coordinates_field );

    fem_meta.declare_field_relation(
      elem_node_rot ,
      stk::mesh::fem::get_element_node_stencil(SpatialDim) ,
      rotation_field );

    // Declare the size of the aggressive "gather" field
    //     double * elem_node_coord[ size = number_of_nodes ]
    // is the number of nodes per element.
    // This size is different for each element block.

    stk::mesh::put_field(elem_node_coord , element_rank , block_hex , shards::Hexahedron<> ::node_count );

    stk::mesh::put_field(elem_node_coord, element_rank, block_quad_shell,shards::ShellQuadrilateral<> ::node_count );

    stk::mesh::put_field(elem_node_rot, element_rank, block_quad_shell,shards::ShellQuadrilateral<> ::node_count );

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

    use_case_7_generate_mesh(
      mesh_options ,
      mesh_bulk_data ,
      coordinates_field ,
      elem_node_coord ,
      block_hex ,
      block_quad_shell );

    use_case_7_initialize_data(
      mesh_bulk_data ,
      coordinates_field ,
      displacements_field ,
      rotation_field );

    time_max[0] = stk::wall_dtime( wtime );

    mesh_bulk_data.modification_end();

    time_max[1] = stk::wall_dtime( wtime );

    //------------------------------------------------------------------

#ifdef USE_GNU_MALLOC_HOOKS
    if (stk::parallel_machine_rank(comm) == 0) {
      double net_alloc = alloc_MB() - freed_MB();
      std::cout << "Mesh creation:" << "\n   Total allocated: "
        << alloc_MB()<<"MB in "<<alloc_blks() << " blocks."
        << "\n   Total freed: " << freed_MB() << "MB in "
        << freed_blks() << " blocks."
        << "\n   Net allocated: "<<net_alloc << "MB."<<std::endl;
    }
#endif

    //------------------------------------------------------------------

    stk::mesh::Selector selector_all(universal);
    stk::mesh::Selector selector_hex(block_hex);
    stk::mesh::count_entities( selector_hex, mesh_bulk_data, hex_block_count );
    stk::mesh::Selector selector_quad_shell(block_quad_shell);
    stk::mesh::count_entities( selector_quad_shell, mesh_bulk_data, shell_block_count);

    verify_elem_node_coord( mesh_bulk_data, elem_node_coord, coordinates_field);

    const std::vector< stk::mesh::Bucket * >& node_buckets = mesh_bulk_data.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );
    const std::vector< stk::mesh::Bucket * >& elem_buckets = mesh_bulk_data.buckets( element_rank );

    ElementMeanValue elem_mean_hex(coordinates_field, elem_node_coord);
    ElementMeanValue_Gather elem_mean_gather_hex(coordinates_field, elem_node_coord);

    ElementMeanValue elem_mean_quad_shell(rotation_field, elem_node_rot);
    ElementMeanValue_Gather elem_mean_gather_quad_shell(rotation_field, elem_node_rot);

    double a_scale = 0.5 ;
    double b_scale = 0.5 ;
    double c_scale = 0.5 ;

    NodeScaleSum node_scale_sum(a_scale, residual_field,
                                b_scale, displacements_field,
                                c_scale, rotation_field);

    //------------------------------------------------------------------
    // Timing tests:

    for ( unsigned i = 0 ; i < num_trials ; ++i ) {

      wtime = stk::wall_time();

      // None of these Selector usages have a union part vector so we
      // pass in an empty union vector.
      mesh::PartVector empty_union_vector;

      alg_runner->run( selector_hex , empty_union_vector, elem_buckets , elem_mean_hex );

      time_max[2] += stk::wall_dtime( wtime );

      alg_runner->run( selector_hex , empty_union_vector, elem_buckets , elem_mean_gather_hex );

      time_max[3] += stk::wall_dtime( wtime );

      // Mean rotation of shell elements:

      alg_runner->run( selector_quad_shell , empty_union_vector, elem_buckets , elem_mean_quad_shell );

      time_max[4] += stk::wall_dtime( wtime );

      alg_runner->run( selector_quad_shell , empty_union_vector, elem_buckets , elem_mean_gather_quad_shell );

      time_max[5] += stk::wall_dtime( wtime );

      alg_runner->run( selector_all , empty_union_vector, node_buckets , node_scale_sum );

      time_max[6] += stk::wall_dtime( wtime );

      double result = dot(*alg_runner , mesh_bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK , residual_field , displacements_field);

      time_max[7] += stk::wall_dtime( wtime );

      if ( ! i ) { result_inner_product[0] = result ; }
      if ( i + 1 == num_trials ) { result_inner_product[1] = result ; }
    }

    //------------------------------

    wtime = stk::wall_time();
  }

  time_max[9] = stk::wall_dtime( wtime );

  bool success = true;

  if (mesh_options == "1000x1x1000+shell:y" && num_trials == 100) {
    double expected_final_inner_prod = 1.4035038e+10;
    double err = std::abs(expected_final_inner_prod - result_inner_product[1])/expected_final_inner_prod;
    success = err < 1.e-6;
    if ( !success ) {
      std::cout << "TEST FAILED: Inner-product Err="<<err << ", inner-prod: " << result_inner_product[1] << std::endl;
    }
  }
  else if (mesh_options == "50x50x500+shell:z" && num_trials == 100) {
    double expected_final_inner_prod = 1108943329.35;
    double err = std::abs(expected_final_inner_prod - result_inner_product[1])/expected_final_inner_prod;
    success = err < 1.e-6;
    if ( !success ) {
      std::cout << "TEST FAILED: Inner-product Err="<<err << ", inner-prod: " << result_inner_product[1] << std::endl;
    }
  }
  else if (mesh_options == "100x100x100+shell:xXyYzZ" && num_trials == 100) {
    double expected_final_inner_prod = 104245300.5;
    double err = std::abs(expected_final_inner_prod - result_inner_product[1])/expected_final_inner_prod;
    success = err < 1.e-6;
    if ( !success ) {
      std::cout << "TEST FAILED: Inner-product Err="<<err << ", inner-prod: " << result_inner_product[1] << std::endl;
    }
  }

  time_min[0] = time_max[0] ;
  time_min[1] = time_max[1] ;
  time_min[2] = time_max[2] ;
  time_min[3] = time_max[3] ;
  time_min[4] = time_max[4] ;
  time_min[5] = time_max[5] ;
  time_min[6] = time_max[6] ;
  time_min[7] = time_max[7] ;
  time_min[8] = time_max[8] ;
  time_min[9] = time_max[9] ;

  stk::all_reduce( comm , stk::ReduceMax<10>( time_max ) & stk::ReduceMin<10>( time_min ) );

  time_max[2] /= num_trials ;
  time_max[3] /= num_trials ;
  time_max[4] /= num_trials ;
  time_max[5] /= num_trials ;
  time_max[6] /= num_trials ;
  time_max[7] /= num_trials ;

  time_min[2] /= num_trials ;
  time_min[3] /= num_trials ;
  time_min[4] /= num_trials ;
  time_min[5] /= num_trials ;
  time_min[6] /= num_trials ;
  time_min[7] /= num_trials ;

  double centroid_grind_time = 0 ;
  double mean_rot_grind_time = 0 ;

  if ( hex_block_count[element_rank] + shell_block_count[element_rank] ) {
    centroid_grind_time = 1e6 *
      time_max[2] / ( hex_block_count[element_rank] + shell_block_count[element_rank] );
  }

  if ( shell_block_count[element_rank] ) {
    mean_rot_grind_time = 1e6 *
      time_max[4] / ( shell_block_count[element_rank] );
  }

  if ( ! stk::parallel_machine_rank( comm ) ) {
    std::cout
      << "stk_mesh performance use case results:" << std::endl
      << "  Number trials       = " << num_trials << std::endl
      << "  Mesh Options        = '" << mesh_options << "', "  << std::endl
                                    << std::endl
      << "  Inner-prod trial #1 = " << result_inner_product[0]
                                    << std::endl
      << "  Inner-prod trial #"     << num_trials
                                    << " = " << result_inner_product[1]
                                    << std::endl
      << "  Mesh generation     = " << time_min[0] << " : "
                                    << time_max[0] << " sec, min : max"
                                    << std::endl
      << "  Mesh synchronization = " << time_min[1] << " : "
                                    << time_max[1] << " sec, min : max"
                                    << std::endl
      << "  Elem centroid       = " << time_min[2] << " : "
                                    << time_max[2] << " sec, min : max"
                                    << std::endl
      << "  Elem centroid gather= " << time_min[3] << " : "
                                    << time_max[3] << " sec, min : max"
                                    << std::endl
      << "  Elem centroid grind = " << centroid_grind_time
                                    << " millisec / element"
                                    << std::endl
      << "  Elem mean-rot       = " << time_min[4] << " : "
                                    << time_max[4] << " sec, min : max"
                                    << std::endl
      << "  Elem mean-rot gather= " << time_min[5] << " : "
                                    << time_max[5] << " sec, min : max"
                                    << std::endl
      << "  Elem mean-rot grind = " << mean_rot_grind_time
                                    << " millisec / shell"
                                    << std::endl
      << "  Node scaled sum     = " << time_min[6] << " : "
                                    << time_max[6] << " sec, min : max"
                                    << std::endl
      << "  Node inner-prod     = " << time_min[7] << " : "
                                    << time_max[7] << " sec, min : max"
                                    << std::endl
      << "  Mesh output         = " << time_min[8] << " : "
                                    << time_max[8] << " sec, min : max"
                                    << std::endl
      << "  Mesh destruction    = " << time_min[9] << " : "
                                    << time_max[9] << " sec, min : max"
                                    << std::endl
      << std::endl ;
    std::cout << "Output for wiki table:"<<std::endl;
    std::cout << "||||||stk_mesh||"
           << time_min[0] << "||"
           << time_min[2]*num_trials << "||"
           << time_min[4]*num_trials << "||"
           << time_min[6]*num_trials << "||"
           << time_min[7]*num_trials << "||" << std::endl;
  }
  return success;
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

void use_case_7_initialize_data(
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

} // namespace app
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <generated/Iogn_GeneratedMesh.h>

namespace stk {
namespace app {

void use_case_7_generate_mesh(
  const std::string& mesh_options ,
  stk::mesh::BulkData & mesh ,
  const VectorFieldType & node_coord ,
  const ElementNodePointerFieldType & /*elem_node_coord */,
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
              node_map[ local_node_id[0] - 1 ] ,
              node_map[ local_node_id[1] - 1 ] ,
              node_map[ local_node_id[2] - 1 ] ,
              node_map[ local_node_id[3] - 1 ] ,
              node_map[ local_node_id[4] - 1 ] ,
              node_map[ local_node_id[5] - 1 ] ,
              node_map[ local_node_id[6] - 1 ] ,
              node_map[ local_node_id[7] - 1 ]
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
              node_map[ local_node_id[0] - 1 ] ,
              node_map[ local_node_id[1] - 1 ] ,
              node_map[ local_node_id[2] - 1 ] ,
              node_map[ local_node_id[3] - 1 ]
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

} // namespace app
} // namespace stk

