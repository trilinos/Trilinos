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
#include <stk_mesh/base/Property.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

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

bool use_case_14_driver(MPI_Comm comm,
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
    std::cout << " stk_mesh Use Case #14 - element internal force, begin" << std::endl ;
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

  double dt  = 1.0e-03;
  double YM  = 1e7; // Young's Modulus
  double PR  = 0.33; // Poisson Ratio
  std::string material_model_name("ELASTIC");

  //----------------------------------
  // Compute other material constants.
  double BM = YM/(3.0-6.0*PR); // Bulk Modulus
  double SM = YM/(2.0+2.0*PR); // Shear Modulus
  double TM = 2.0 * SM;
  double LAMBDA = YM*PR/((1.0+PR)*(1.0-2.0*PR));
  double DM = TM + LAMBDA; // Dilatational Modulus

  lame::MatProps materialProperties;

  std::vector<double> youngs_modulus; youngs_modulus.push_back(YM);
  std::vector<double> poissons_ratio; poissons_ratio.push_back(PR);
  std::vector<double> bulk_modulus;     bulk_modulus.push_back(BM);
  std::vector<double> shear_modulus;   shear_modulus.push_back(SM);
  std::vector<double> two_mu;                 two_mu.push_back(TM);
  std::vector<double> lambda;                 lambda.push_back(LAMBDA);
  std::vector<double> dilatational_modulus; dilatational_modulus.push_back(DM);

  materialProperties["YOUNGS_MODULUS"] = youngs_modulus;
  materialProperties["POISSONS_RATIO"] = poissons_ratio;
  materialProperties["BULK_MODULUS"]   = bulk_modulus;
  materialProperties["SHEAR_MODULUS"]  = shear_modulus;
  materialProperties["TWO_MU"]         = two_mu;
  materialProperties["LAMBDA"]         = lambda;
  materialProperties["DILATATIONAL_MODULUS"] = dilatational_modulus;

  lame::Material * matmodel = lame::Elastic::createMaterial( materialProperties );

  //----------------------------------

  APSHex8ug hex_element;

  lame::matParams materialParameters;
  materialParameters.dt = dt;

  // Timing:
  //   [0] = stk::mesh::MetaData creation
  //   [1] = stk::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force
  //   [4] = Parallel swap-add of internal force

  double time_min[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double time_max[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double wtime = 0 ;

  //--------------------------------------------------------------------

  reset_malloc_stats();

  if ( 0 == stk::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh performance use case #14" << std::endl
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

    stk::mesh::fem::FEMMetaData meta_data( spatial_dimension );
    stk::io::MeshData mesh_data;
    std::string filename = working_directory + mesh_filename;
    stk::io::create_input_mesh(mesh_type, filename, comm,
				     meta_data, mesh_data);
    stk::io::define_input_fields(mesh_data, meta_data);

    stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_1", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_2", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_10", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_20", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_30", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }
    {
      stk::mesh::Part & block_hex        = meta_data.declare_part("block_40", meta_data.element_rank());
      stk::mesh::fem::set_cell_topology( block_hex, hex_top );
    }

    // Add property to each hex element block...
    {
      stk::mesh::Property<double>& delta_t = meta_data.get_meta_data(meta_data).declare_property<double>("dt");
      stk::mesh::Property<lame::MatProps>& mprops =
	      meta_data.get_meta_data(meta_data).declare_property<lame::MatProps>("materialProperties");

      // All parts of the meta data:
      const mesh::PartVector & all_parts = meta_data.get_parts();
      for (mesh::PartVector::const_iterator i = all_parts.begin();
           i != all_parts.end(); ++i) {

        mesh::Part * const part = *i ;
	if ( part->primary_entity_rank() == meta_data.element_rank()) {
          const CellTopologyData *const topology = meta_data.get_cell_topology(*part).getCellTopologyData();
	  if (topology->key == shards::Hexahedron<8>::key) {

	    meta_data.get_meta_data(meta_data).put_property( delta_t, *part);
	    double* delta_t_ptr = stk::mesh::property_data(delta_t, *part);
	    *delta_t_ptr = dt;

	    meta_data.get_meta_data(meta_data).put_property(mprops, *part);
	    lame::MatProps* mprops_ptr = stk::mesh::property_data(mprops, *part);
	    *mprops_ptr = materialProperties;
	  }
	}
      }
    }

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
      filename = working_directory + mesh_filename + ".out14";
      stk::io::create_output_mesh(filename, comm, bulk_data, mesh_data);
      stk::io::define_output_fields(mesh_data, meta_data, true);
    }

    stk::app::use_case_14_initialize_nodal_data(bulk_data ,
                                                *fields.model_coordinates ,
                                                *fields.coordinates_field ,
                                                *fields.velocity_field,
                                                dt);

    time_max[1] = stk::wall_dtime( wtime );

    //------------------------------------------------------------------
    // Ready to run the algorithms:
    //------------------------------------------------------------------

    stk::mesh::PartVector hex_parts;

    bool skip_topology_root_parts = true;
    stk::mesh::fem::get_parts_with_topology<shards::Hexahedron<8> >(bulk_data, hex_parts, skip_topology_root_parts);

    // Reference to all node and element buckets,
    // will select which ones we need later

    const std::vector< stk::mesh::Bucket * > & element_buckets = bulk_data.buckets( meta_data.element_rank());
    const std::vector< stk::mesh::Bucket * > & node_buckets = bulk_data.buckets( meta_data.node_rank() );

    // Selectors for buckets:

    stk::mesh::Selector select_owned_hexahedrons =
      meta_data.locally_owned_part() & selectUnion(hex_parts);


    stk::mesh::Selector select_used =
      meta_data.locally_owned_part() |
      meta_data.globally_shared_part() ;

    // Initialize bucket fields:

    use_case_14_initialize_element_fields(fields, select_owned_hexahedrons, element_buckets, YM, PR);

    //------------------------------------------------------------------
    time_max[2] = stk::wall_dtime( wtime );
    //------------------------------------------------------------------

    wtime = stk::wall_time();

    MyHexInternalForceAlg elem_alg(fields, materialParameters,
				   matmodel, meta_data.get_meta_data(meta_data) );
    MyNodalForceScatterAlg node_alg(fields, meta_data.element_rank());

    for(int n=0; n<num_trials; ++n) {
      //
      // Call Internal Force!!!
      //

      wtime = stk::wall_time();

      // Need to zero out the old accumulated internal force so that it
      // does not pollute the new accumulation.
      fill( *alg_runner, bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK , *fields.fint_field, 0.0 );

      alg_runner->run_parts( select_owned_hexahedrons , hex_parts, element_buckets , elem_alg );

      stk::mesh::PartVector empty_union_vector;
      alg_runner->run( select_used , empty_union_vector, node_buckets , node_alg );

      time_max[3] += stk::wall_dtime( wtime );

      stk::mesh::parallel_reduce( bulk_data , stk::mesh::sum(*fields.fint_field) );

      time_max[4] += stk::wall_dtime( wtime );

      if (output) {
        stk::io::process_output_request(mesh_data, bulk_data, n);
      }

    }//end for(..num_trials...
    delete matmodel;

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
      << "  Number trials       = " << num_trials << std::endl
      << "  Meta-data setup     = " << time_min[0] << " : "
      << time_max[0] << " sec, min : max"
      << std::endl
      << "  Bulk-data generation= " << time_min[1] << " : "
      << time_max[1] << " sec, min : max"
      << std::endl
      << "  Initialization      = " << time_min[2] << " : "
      << time_max[2] << " sec, min : max"
      << std::endl
      << "  Internal force      = " << time_min[3] << " : "
      << time_max[3] << " sec, min : max"
      << std::endl
      << "  Internal force (total) = " << time_min[3]*num_trials
      << std::endl
      << "  Swap-add            = " << time_min[4] << " : "
      << time_max[4] << " sec, min : max"
      << std::endl
      << "  Mesh destruction    = " << time_min[8] << " : "
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

MyHexInternalForceAlg::MyHexInternalForceAlg(
  const stk::app::Fields &    fields,
  lame::matParams &           matParameters,
  lame::Material *            matModel,
  stk::mesh::MetaData &       md)
  : hex_element(),
    materialParameters(matParameters),
    matmodel(matModel),
    m_fields(fields)
{
  // element-block property
  delta_t = md.get_property<double>("dt");
  m_materialProperties = md.get_property<lame::MatProps>("materialProperties");
}

void
MyHexInternalForceAlg::apply(
  stk::mesh::Bucket::iterator   ibegin ,
  stk::mesh::Bucket::iterator   iend ,
  const stk::mesh::PartVector & parts ) const
{
  enum { num_nodes = 8,
         num_nodes_X_3 = 24 };

  if (parts.size() != 1) {
    std::cout << "MyHexInternalForceAlg::apply ERROR, expect parts.size() == 1."<<std::endl;
  }

  const stk::mesh::Part& part = *(parts[0]);

  const int num_elements = iend - ibegin;

  double *mid_hg_op_ptr                 = stk::mesh::field_data( *m_fields.MidHourglassOp,          ibegin);
  double *material_eff_twomu_ptr        = stk::mesh::field_data( *m_fields.Material_eff_twomu,      ibegin);
  double *material_eff_bulk_modulus_ptr = stk::mesh::field_data( *m_fields.Material_eff_bulk_mod,   ibegin);
  double *mid_step_volume_ptr           = stk::mesh::field_data( *m_fields.Midstep_volume,          ibegin);
  double *element_time_step_ptr         = stk::mesh::field_data( *m_fields.Element_time_step,       ibegin);
  double *element_mass_ptr              = stk::mesh::field_data( *m_fields.Element_mass,            ibegin);
  double *hg_energy_ptr                 = stk::mesh::field_data( *m_fields.Hourglass_energy,        ibegin);
  double *internal_energy_ptr           = stk::mesh::field_data (*m_fields.Internal_energy,         ibegin);
  double *shear_modulus_ptr             = stk::mesh::field_data( *m_fields.Shear_Modulus,           ibegin);
  double *dilatational_modulus_ptr      = stk::mesh::field_data( *m_fields.Dilatational_Modulus,    ibegin );

  double *rotation_old_ptr              = stk::mesh::field_data( *m_fields.RotationOld,             ibegin);
  double *rotation_new_ptr              = stk::mesh::field_data( *m_fields.RotationNew,             ibegin);

  double *stretch_ptr                   = stk::mesh::field_data( *m_fields.Stretch,                 ibegin);
  double *strain_rate_ptr               = stk::mesh::field_data( *m_fields.StrainRate,              ibegin);
  double *stress_old_ptr                = stk::mesh::field_data( *m_fields.StressOld,               ibegin);
  double *stress_new_ptr                = stk::mesh::field_data( *m_fields.StressNew,               ibegin);
  double *rotated_stress_ptr            = stk::mesh::field_data( *m_fields.RotatedStress,           ibegin);

  double *hg_resistance_old_ptr         = stk::mesh::field_data( *m_fields.HourglassResistanceOld,  ibegin);
  double *hg_resistance_new_ptr         = stk::mesh::field_data( *m_fields.HourglassResistanceNew,  ibegin);

  double *vorticity_ptr                 = stk::mesh::field_data( *m_fields.Vorticity,               ibegin);

  const double **coord = (const double **)stk::mesh::field_data( *m_fields.coord_gather,            ibegin);
  double *force_new                     = stk::mesh::field_data( *m_fields.force_new_field,         ibegin);
  const double **velocity               = (const double **) stk::mesh::field_data( *m_fields.velocity_gather,  ibegin);

  // scratch space for gather:
  double elem_coord[num_nodes_X_3* maximum_entity_count ];
  double elem_vel[num_nodes_X_3* maximum_entity_count ];

  //--------------------------------
  { // Gather nodal data into contiguous arrays.

    // element-node pointer fields
    const double ** field_coord  = coord ;
    const double ** vnodes       = velocity ;

    const int num_elem_nodes = num_elements*num_nodes ;

    double* elem_coord_ptr = elem_coord;
    double* elem_vel_ptr   = elem_vel;

    for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
      const double * const f_coord = *field_coord ;
      const double * const f_vel   = *vnodes ;

      elem_coord_ptr[0] = f_coord[0];
      elem_coord_ptr[1] = f_coord[1];
      elem_coord_ptr[2] = f_coord[2];

      elem_vel_ptr[0] = f_vel[0];
      elem_vel_ptr[1] = f_vel[1];
      elem_vel_ptr[2] = f_vel[2];

      ++field_coord ;
      ++vnodes ;
      elem_coord_ptr += 3 ;
      elem_vel_ptr += 3 ;
    }

    //--------------------------------

    lame::matParams localMaterialParameters;

    double* dt = stk::mesh::property_data(*delta_t, part);
    lame::MatProps *materialProperties = stk::mesh::property_data(*m_materialProperties, part);

    double current_stable_time_step = *dt;
    localMaterialParameters.dt   = *dt ;
    localMaterialParameters.nelements   = num_elements ;
    localMaterialParameters.strain_rate = strain_rate_ptr;
    localMaterialParameters.stress_old  = stress_old_ptr;
    localMaterialParameters.stress_new  = stress_new_ptr;

    hex_element.internalForce( num_elements ,
                               *dt,
                               current_stable_time_step,
                               element_time_step_ptr,
                               *matmodel,
                               localMaterialParameters,
                               *materialProperties,
                               & elem_coord[0],
                               & elem_vel[0],
                               rotation_old_ptr,
                               rotation_new_ptr,
                               mid_step_volume_ptr,
                               vorticity_ptr,
                               stretch_ptr,
                               strain_rate_ptr,
                               mid_hg_op_ptr,
                               stress_old_ptr,
                               stress_new_ptr,
                               rotated_stress_ptr,
                               material_eff_bulk_modulus_ptr,
                               material_eff_twomu_ptr,
                               shear_modulus_ptr,
                               dilatational_modulus_ptr,
                               element_mass_ptr,
                               & force_new[0],
                               hg_energy_ptr,
                               internal_energy_ptr,
                               hg_resistance_old_ptr,
                               hg_resistance_new_ptr
                               );
  }
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------

MyNodalForceScatterAlg::MyNodalForceScatterAlg(const stk::app::Fields &fields,
                                               stk::mesh::EntityRank element_rank)
  : m_elementRank(element_rank)
{
  fint_field      = fields.fint_field;
  force_new_field = fields.force_new_field;
}


void
MyNodalForceScatterAlg::apply(
  stk::mesh::Bucket::iterator ibegin ,
  stk::mesh::Bucket::iterator iend ) const
{
  // Scatter internal force values.

  double *fint = stk::mesh::field_data( *fint_field, ibegin);
  size_t num_dof = 3;

  for ( stk::mesh::Bucket::iterator i = ibegin ; i != iend ; ++i ) {
    stk::mesh::Entity & node = *i;
    stk::mesh::PairIterRelation elem_relations = node.relations(m_elementRank);

    for( std::vector<stk::mesh::Relation>::const_iterator
         iter =  elem_relations.first;
         iter != elem_relations.second; ++iter) {

      stk::mesh::Entity* elem = iter->entity();

      int local_node = iter->identifier();

      int elem_offset = local_node*num_dof;

      const double* f_new = field_data( *force_new_field, *elem ) + elem_offset;

      fint[0] += f_new[0];
      fint[1] += f_new[1];
      fint[2] += f_new[2];
    }
    fint += num_dof;
  }
}

