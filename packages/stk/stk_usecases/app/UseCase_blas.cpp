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

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/util/string_case_compare.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_algsup/AlgorithmRunner.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <common/gnu_malloc_hooks.hpp>
#include <app/AppUseCase_14.hpp>
#include <app/UseCase_14_Common.hpp>
#include <app/UseCase_14_Fields.hpp>

#include <app/UseCase_blas_algs.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 14: internal force computation.
// The function 'use_case_14_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk {
namespace app {

static const size_t spatial_dimension = 3;

bool use_case_blas_driver(MPI_Comm comm,
                        int num_threads,
                        int num_trials,
                        const std::string &working_directory,
                        const std::string &mesh_filename,
                        const std::string &mesh_type,
                        const std::string &thread_runner,
                        int bucket_size,
                        bool performance_test)
{
  bool output = !performance_test; // If running for performance measurements, turn off output

  if (stk::parallel_machine_rank(comm) == 0) {
    std::cout << " stk_mesh Use Case Blas - fill, axpby, dot, norm , begin" << std::endl ;
    std::cout << "Running '" << mesh_filename << "' case, num_trials = "
              << num_trials << std::endl;
  }


  const AlgorithmRunnerInterface* alg_runner = NULL ;
  if ( thread_runner.empty() ||
       thread_runner == std::string("NonThreaded") ) {
    alg_runner = stk::algorithm_runner_non_thread();
  }
  else if ( thread_runner == std::string("TPI") ) {
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

  //----------------------------------

  // Timing:
  //   [0] = stk::mesh::MetaData creation
  //   [1] = stk::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = fill and axpby
  //   [4] = dot and norm2

  double time_min[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double time_max[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double wtime = 0 ;

  //--------------------------------------------------------------------

  reset_malloc_stats();

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh performance use case BLAS" << std::endl
              << "  Number Processes = " << stk::parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  {
    wtime = stk::wall_time();

    //------------------------------------------------------------------
    // Declare the mesh meta data: element blocks and associated fields

    stk::mesh::fem::FEMMetaData meta_data(  spatial_dimension );
    stk::io::MeshData mesh_data;
    std::string filename = working_directory + mesh_filename;
    stk::io::create_input_mesh(mesh_type, filename, comm,
			       meta_data, mesh_data);
    stk::io::define_input_fields(mesh_data, meta_data);

    Fields fields;
    use_case_14_declare_fields(fields, meta_data.get_meta_data(meta_data));

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    meta_data.commit();

    //------------------------------------------------------------------

    time_max[0] = stk::wall_dtime( wtime );

    //------------------------------------------------------------------
    // stk::mesh::BulkData bulk data conforming to the meta data.
    stk::mesh::BulkData bulk_data(meta_data.get_meta_data(meta_data) , comm, bucket_size);
    stk::io::populate_bulk_data(bulk_data, mesh_data);

    //------------------------------------------------------------------
    // Create output mesh...  (input filename + ".out14")
    if (output) {
      filename = working_directory + mesh_filename + ".blas";
      stk::io::create_output_mesh(filename, comm, bulk_data, mesh_data);
      stk::io::define_output_fields(mesh_data, meta_data, true);
    }

    stk::app::use_case_14_initialize_nodal_data(bulk_data ,
                                                *fields.model_coordinates ,
                                                *fields.coordinates_field ,
                                                *fields.velocity_field,
                                                1.0 /*dt*/);

    time_max[1] = stk::wall_dtime( wtime );

    //------------------------------------------------------------------
    // Ready to run the algorithms:
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    time_max[2] = stk::wall_dtime( wtime );
    //------------------------------------------------------------------

    wtime = stk::wall_time();

    double dot1 = 0;

    for(int n=0; n<num_trials; ++n) {
      //
      // Call BLAS algs.
      //

      wtime = stk::wall_time();

      fill( *alg_runner, bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK , *fields.velocity_field, 0.2 );

      fill( *alg_runner, bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK , *fields.fint_field, 1.0 );

      axpby( *alg_runner, bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK ,
             0.01, *fields.model_coordinates , 1.0 , *fields.coordinates_field );

      axpby( *alg_runner, bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK ,
             0.1, *fields.coordinates_field, 1.0 , *fields.velocity_field );

      time_max[3] += stk::wall_dtime( wtime );

      dot1 = dot( *alg_runner, bulk_data, stk::mesh::fem::FEMMetaData::NODE_RANK ,
                  *fields.velocity_field, *fields.coordinates_field );

      double dot2 = dot( *alg_runner, bulk_data, stk::mesh::fem::FEMMetaData::NODE_RANK,
                         *fields.velocity_field, *fields.fint_field );

      double norm_1 = norm2(*alg_runner, bulk_data, stk::mesh::fem::FEMMetaData::NODE_RANK, *fields.velocity_field );

      double norm_2 = norm2(*alg_runner, bulk_data, stk::mesh::fem::FEMMetaData::NODE_RANK, *fields.coordinates_field );

      if ( stk::parallel_machine_rank( comm ) == 0 ) {
        std::cout << "    " << dot1 << "  " << dot2 << "  " << norm_1 << "  " << norm_2 << std::endl;
      }

      time_max[4] += stk::wall_dtime( wtime );

      if (output) {
        stk::io::process_output_request(mesh_data, bulk_data, n);
      }

    }//end for(..num_trials...

    if ( stk::parallel_machine_rank( comm ) == 0 ) {
      //Try to make sure the number gets printed out just the way we want it,
      //so we can use it as a pass/fail check for a regression test...
      std::cout.precision(6);
      std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
      std::cout << "Final dot1: " << dot1 << std::endl;
    }
    //------------------------------------------------------------------

#ifdef USE_GNU_MALLOC_HOOKS
    if (parallel_machine_rank(comm) == 0) {
      double net_alloc = alloc_MB() - freed_MB();
      std::cout << "Mesh creation:" << "\n   Total allocated: "
                << alloc_MB()<<"MB in "<<alloc_blks() << " blocks."
                << "\n   Total freed: " << freed_MB() << "MB in "
                << freed_blks() << " blocks."
                << "\n   Net allocated: "<<net_alloc << "MB."<<std::endl;
    }
#endif

    //------------------------------------------------------------------
  }

  time_max[8] = stk::wall_dtime( wtime );

  time_min[0] = time_max[0] ;
  time_min[1] = time_max[1] ;
  time_min[2] = time_max[2] ;
  time_min[3] = time_max[3] ;
  time_min[4] = time_max[4] ;
  time_min[5] = time_max[5] ;
  time_min[6] = time_max[6] ;
  time_min[7] = time_max[7] ;
  time_min[8] = time_max[8] ;

  stk::all_reduce( comm , stk::ReduceMax<9>( time_max ) & stk::ReduceMin<9>( time_min ) );

  time_max[3] /= num_trials ;
  time_max[4] /= num_trials ;
  time_max[5] /= num_trials ;
  time_max[6] /= num_trials ;

  time_min[3] /= num_trials ;
  time_min[4] /= num_trials ;
  time_min[5] /= num_trials ;
  time_min[6] /= num_trials ;

  //   [0] = stk::mesh::MetaData creation
  //   [1] = stk::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force

  if ( ! stk::parallel_machine_rank( comm ) ) {
    std::cout
      << "stk_mesh performance use case results:" << std::endl
      << "  Number of trials         = " << num_trials << std::endl
      << "  Meta-data setup          = " << time_min[0] << " : "
      << time_max[0] << " sec, min : max"
      << std::endl
      << "  Bulk-data generation     = " << time_min[1] << " : "
      << time_max[1] << " sec, min : max"
      << std::endl
      << "  Initialization           = " << time_min[2] << " : "
      << time_max[2] << " sec, min : max"
      << std::endl
      << "  fill & axpby (per-trial) = " << time_min[3] << " : "
      << time_max[3] << " sec, min : max"
      << std::endl
      << "  dot & norm2 (per-trial)  = " << time_min[4] << " : "
      << time_max[4] << " sec, min : max"
      << std::endl
      << "  Mesh destruction         = " << time_min[8] << " : "
      << time_max[8] << " sec, min : max"
      << std::endl
      << std::endl ;
  }

  return true;
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

} // namespace app
} // namespace stk

