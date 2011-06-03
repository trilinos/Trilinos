/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author K. H. Pierson  <hcedwar@sandia.gov>
 * @date   April 2009
 */

#include <stdexcept>
#include <sstream>
#include <vector>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/SimpleArrayOps.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/Comm.hpp>

// #include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <common/gnu_malloc_hooks.hpp>
#include <mesh/UseCase_AD.hpp>

//
// How do I modify this code to find this header?
//
//#include <sacado/Sacado.hpp>

// 2. Template the outer methods.
//
// 3. Write another driver to compute the derivative (vs. internal force
// call).
//
// 4.
// DFad -> dynamic FAD ( dynamic memory allocation)
//
// typedef Sacado::Fad::DFad<Real> FAD_Type;

//----------------------------------------------------------------------
// This file contains the implementation of use-case AD: internal force computation.
// The function 'use_case_AD_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace stk {
namespace app {

enum { SpatialDim = 3 };

inline stk::mesh::EntityRank get_element_rank(const stk::mesh::fem::FEMMetaData& meta_data)
{
  return meta_data.element_rank();
}

inline stk::mesh::EntityRank get_element_rank(const stk::mesh::Part& part)
{
  return get_element_rank(stk::mesh::fem::FEMMetaData::get(part));
}

VectorField &
declare_vector_field_on_all_nodes(
  stk::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk::mesh::put_field( meta_data.declare_field<VectorField>(s), stk::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() , n1 );
}


VectorField &
declare_vector_field_on_all_elements(
  stk::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk::mesh::put_field( meta_data.declare_field<VectorField>(s), get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


ScalarField &
declare_scalar_field_on_all_elements(
  stk::mesh::fem::FEMMetaData & meta_data , const std::string & s )
{
  return stk::mesh::put_field( meta_data.declare_field<ScalarField>(s), get_element_rank(meta_data) , meta_data.universal_part() );
}


SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  stk::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return put_field( meta_data.declare_field<SymmetricTensorField>(s), get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


template< typename Type , class T1 >
stk::mesh::Field<Type,T1> &
put_field_on_elements( stk::mesh::Field<Type,T1> & f , stk::mesh::Part & p , unsigned n1 )
{
  stk::mesh::put_field( f , get_element_rank(p) , p , n1 );
  return f ;
}


template< typename Type , class T1 >
stk::mesh::Field<Type,T1> & put_field_on_all_elements( stk::mesh::Field<Type,T1> & f , unsigned n1 )
{
  put_field_on_elements( f , stk::mesh::MetaData::get(f).universal_part() , n1 );
  return f ;
}

void use_case_AD_driver(
  MPI_Comm comm ,
  const std::string & file_name ,
  const unsigned box_size[] ,
  const unsigned box_sides[][2] ,
  const unsigned num_trials );

//--------------------------------------------------------------------
//
// main driver for use-case AD: element internal force.
//

void use_case_AD_driver( MPI_Comm comm , bool performance_test )
{
  int num_procs = stk::parallel_machine_size( comm );

  if ( ! stk::parallel_machine_rank( comm ) ) {
    std::cout << " stk_mesh Use Case #AD - element internal force, begin" << std::endl ;
  }

  if ( ! performance_test ) {
    // Quick & small test for correctness:
    const unsigned a_box_size[3] = { 10 , 10 , 10*num_procs };
    const unsigned a_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    int num_trials = 1;
    use_case_AD_driver( comm , "use_case_14.exo" , a_box_size , a_box_sides , num_trials );
  }
  else {
    int num_trials = 20 ; // 582 ;

    const unsigned a_box_size[3] = { 100 , 100 , 100 };
    const unsigned a_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

    std::cout << "Running 100x100x100 case, num_trials = " << num_trials << std::endl;

    use_case_AD_driver( comm , "use_case_14a.exo" , a_box_size , a_box_sides , num_trials );

//     const unsigned b_box_size[3] = { 100 , 125 , 10*num_procs };
//     const unsigned b_box_sides[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

//     std::cout << "Running 100x125x10*nprocs case, num_trials = " << num_trials << std::endl;

//     use_case_AD_driver( comm , "use_case_14b.exo" , b_box_size , b_box_sides , num_trials );
  }
}

//--------------------------------------------------------------------

inline
void zero_data( double * beg , double * const end )
{ while ( beg < end ) { *beg++ = 0 ; } }

template< class FieldType >
void zero_field_data( stk::mesh::BulkData & mesh , stk::mesh::EntityRank type , const FieldType & field )
{
  typedef stk::mesh::BucketArray< FieldType > array_type ;

  const std::vector<stk::mesh::Bucket*> & ks = mesh.buckets( type );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    stk::mesh::Bucket & bucket = **ik ;
    array_type data( field , bucket );
    zero_data( data.contiguous_data() , data.contiguous_data() + data.size() );
  }
}

//--------------------------------------------------------------------


void use_case_AD_driver(
  MPI_Comm comm ,
  const std::string & /*file_name */,
  const unsigned box_size[] ,
  const unsigned box_sides[][2] ,
  const unsigned num_trials )
{
  double dt  = 1.0e-02;
  double YM  = 1e7; // Young's Modulus
  double PR  = 0.33; // Poisson Ratio
  double element_mass_dummy_value = 1.0;
  std::string material_model_name("ELASTIC");

  double current_stable_time_step = dt;

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

  lame::Material<double> * matmodel = lame::Elastic<double>::createMaterial( materialProperties );

  //----------------------------------

  Hex8ug<double> hex_element;

  lame::matParams<FAD_Type> materialParameters;

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
    std::cout << "stk_mesh performance use case #AD" << std::endl
              << "  Number Processes = " << stk::parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  {
    wtime = stk::wall_time();

    //------------------------------------------------------------------
    // Declare the mesh meta data and bulk data

    stk::mesh::fem::FEMMetaData mesh_meta_data( SpatialDim );
    const stk::mesh::EntityRank element_rank = mesh_meta_data.element_rank();
    stk::mesh::BulkData mesh_bulk_data( mesh_meta_data.get_meta_data(mesh_meta_data), MPI_COMM_WORLD, 1000 );

    //--------------------------------
    // Element-block declarations typically occur when reading the
    // mesh-file meta-data, and thus won't usually appear in application code.
    // Declaring the element blocks and associating an element traits
    // with each element block.

    stk::mesh::Part & block_hex = mesh_meta_data.declare_part("block_1", element_rank);
    stk::mesh::fem::CellTopology hex_top(shards::getCellTopologyData<shards::Hexahedron<> >());
    stk::mesh::fem::set_cell_topology( block_hex, hex_top );

    //--------------------------------

    // Nodal vector fields
    VectorField &model_coordinates     = declare_vector_field_on_all_nodes( mesh_meta_data , "model_coordinates" , SpatialDim );
    VectorField &coordinates_field     = declare_vector_field_on_all_nodes( mesh_meta_data , "coordinates" , SpatialDim );
    VectorField &velocity_field        = declare_vector_field_on_all_nodes( mesh_meta_data , "velocity" ,    SpatialDim );
    VectorField &fint_field            = declare_vector_field_on_all_nodes( mesh_meta_data , "force_internal" , SpatialDim );

    // Element vector fields:
    VectorField &Vorticity             = declare_vector_field_on_all_elements( mesh_meta_data , "Vorticity" , SpatialDim );

    // Element scalar fields:
    ScalarField &Shear_Modulus         = declare_scalar_field_on_all_elements( mesh_meta_data , "shear_modulus" );
    ScalarField &Dilatational_Modulus  = declare_scalar_field_on_all_elements( mesh_meta_data , "dilatational_modulus");
    ScalarField &Material_eff_twomu    = declare_scalar_field_on_all_elements( mesh_meta_data , "material_effictive_two_mu");
    ScalarField &Material_eff_bulk_mod = declare_scalar_field_on_all_elements( mesh_meta_data , "material_effective_bulk_moduli");
    ScalarField &Midstep_volume        = declare_scalar_field_on_all_elements( mesh_meta_data , "mid_step_volume");
    ScalarField &Element_time_step     = declare_scalar_field_on_all_elements( mesh_meta_data , "element_time_step");
    ScalarField &Element_mass          = declare_scalar_field_on_all_elements( mesh_meta_data , "element_mass");
    ScalarField &Hourglass_energy      = declare_scalar_field_on_all_elements( mesh_meta_data , "hourglass_energy");
    ScalarField &Internal_energy       = declare_scalar_field_on_all_elements( mesh_meta_data , "internal_energy");

    // Element symmetric tensor fields:
    SymmetricTensorField &Stretch       = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "stretch" ,    6 );
    SymmetricTensorField &StrainRate    = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "StrainRate" , 6 );
    SymmetricTensorField &RotatedStress = declare_symmetric_tensor_field_on_all_elements( mesh_meta_data , "RotatedStress" , 6 );

    //--------------------------------
    // We don't like specifying the template type TWICE! very error prone.
    //--------------------------------
    // The multi-state fields don't have the 'declare and put' convenience functions (yet)
    //
    //  For clarity declare a integer to used for the number of states arguments.
    //
    const unsigned two_states = 2 ;

    // Element two state symmetric tensor field, on all elements (as two function calls):
    SymmetricTensorField &StressNew = mesh_meta_data.declare_field< SymmetricTensorField >( "Stress" , two_states );
    put_field_on_all_elements( StressNew , 6 );

    // Element two state full tensor field on all elements (as nested function calls):
    FullTensorField &RotationNew =
      put_field_on_all_elements( mesh_meta_data.declare_field< FullTensorField >("Rotation", two_states ) , 9 );

    //--------------------------------
    // Hourglass fields, these don't have the convenience functions for 'declare and put'

    HourglassArrayField &HourglassResistanceNew =
      put_field_on_all_elements( mesh_meta_data.declare_field< HourglassArrayField >("HourglassResistance" , two_states ) , 12 );

    HourglassOpField & MidHourglassOp =
      put_field_on_all_elements( mesh_meta_data.declare_field< HourglassOpField >("mid_hourglass_operator") , 32 );

    //--------------------------------
    // Declare aggressive "gather" fields which are an array of
    // pointers to the element's nodes' coordinate, velocity, and
    // internal force field data.
    //
    // The declarations specify element fields of the following form:
    //
    //     double * coord_gather   [ nodes_per_element ]
    //     double * velocity_gather[ nodes_per_element ]
    //     double * fint_gather    [ nodes_per_element ]
    //
    // where
    //
    //     coord_gather[i]    == field_data( coordinates_field , element_node[i] )
    //     velocity_gather[i] == field_data( velocity ,          element_node[i] )
    //     fint_gather[i]     == field_data( fint ,              element_node[i] )
    //
    // The number of nodes per element could vary, so the field is put on each element block
    // with a size of the number of nodes per element in that element block.

    ElementNodePointerField & coord_gather    = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "coord_gather" , coordinates_field );
    ElementNodePointerField & velocity_gather = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "velocity_gather" , velocity_field );
    ElementNodePointerField & fint_gather     = declare_element_node_pointer_field( mesh_meta_data.get_meta_data(mesh_meta_data), "fint_gather" , fint_field );

    put_field_on_elements( coord_gather ,    block_hex , shards::Hexahedron<> ::node_count );
    put_field_on_elements( velocity_gather , block_hex , shards::Hexahedron<> ::node_count );
    put_field_on_elements( fint_gather ,     block_hex , shards::Hexahedron<> ::node_count );

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    mesh_meta_data.commit();

    //------------------------------------------------------------------

    time_max[0] = stk::wall_dtime( wtime );

    //------------------------------------------------------------------

    // In a typical app, the mesh would be read from file at this point.
    // But in this use-case, we generate the mesh and initialize
    // field data to use-case defined values.

    use_case_AD_generate_mesh(
      mesh_bulk_data ,
      box_size ,
      model_coordinates ,
      coord_gather ,
      block_hex ,
      box_sides );

    use_case_AD_initialize_data(
      mesh_bulk_data ,
      model_coordinates ,
      coordinates_field ,
      velocity_field );

    time_max[1] = stk::wall_dtime( wtime );

  //------------------------------------------------------------------
  // Ready to run the algorithms:
  //------------------------------------------------------------------

  stk::mesh::Selector select_owned( stk::mesh::MetaData::get(mesh_bulk_data).locally_owned_part() );

  const std::vector< stk::mesh::Bucket * > & element_buckets =
    mesh_bulk_data.buckets(element_rank);

  unsigned maximum_bucket_size = 0 ;
  for ( std::vector< stk::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned(**k) ) {
    const stk::mesh::Bucket & bucket = **k ;
    if ( maximum_bucket_size < bucket.size() ) { maximum_bucket_size = bucket.size(); }
  }

  // Need both the the old and new states of these two-state fields:

  HourglassArrayField  & HourglassResistanceOld = HourglassResistanceNew.field_of_state( stk::mesh::StateOld );
  SymmetricTensorField & StressOld              = StressNew.field_of_state( stk::mesh::StateOld );
  FullTensorField      & RotationOld            = RotationNew.field_of_state( stk::mesh::StateOld );

  for ( std::vector< stk::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
    const stk::mesh::Bucket & bucket = **k ;

    const int num_elements = bucket.size();

    double * stretch        = stk::mesh::field_data( Stretch,    bucket.begin() );
    double * strain_rate    = stk::mesh::field_data( StrainRate, bucket.begin() );
    double * stress_new     = stk::mesh::field_data( StressNew,  bucket.begin() );
    double * stress_old     = stk::mesh::field_data( StressOld , bucket.begin() );
    double * rotated_stress = stk::mesh::field_data( RotatedStress , bucket.begin() );
    double * rotation_old   = stk::mesh::field_data( RotationOld , bucket.begin() );
    double * rotation_new   = stk::mesh::field_data( RotationNew , bucket.begin() );
    double * mass           = stk::mesh::field_data( Element_mass , bucket.begin() );
    double * hg_old         = stk::mesh::field_data( HourglassResistanceOld , bucket.begin() );
    double * hg_new         = stk::mesh::field_data( HourglassResistanceNew , bucket.begin() );
    double *vorticity_ptr   = stk::mesh::field_data( Vorticity, bucket.begin());
    double *mid_hg_op_ptr   = stk::mesh::field_data( MidHourglassOp, bucket.begin());

    for ( int i = 0 ; i < num_elements ; ++i ) {

      stk::Copy<32>( mid_hg_op_ptr , 0.0 ); mid_hg_op_ptr += 32 ;
      stk::Copy< 3>( vorticity_ptr , 0.0 ); vorticity_ptr += 3 ;
      stk::Copy< 6>( rotated_stress, 0.0 ); rotated_stress += 6 ;
      stk::Copy< 6>( strain_rate,    0.0 ); strain_rate += 6 ;
      stk::Copy<12>( hg_old,         0.0 ); hg_old += 12 ;
      stk::Copy<12>( hg_new,         0.0 ); hg_new += 12 ;
      stk::Copy< 6>( stress_new,     0.0 ); stress_new += 6 ;
      stk::Copy< 6>( stress_old,     0.0 ); stress_old += 6 ;

      mass[i] = element_mass_dummy_value;

      // initialize stretch to identity.
      stretch[0] = 1.0;
      stretch[1] = 1.0;
      stretch[2] = 1.0;
      stretch[3] = 0.0;
      stretch[4] = 0.0;
      stretch[5] = 0.0;
      stretch+=6;

      // initialize rotation to identity.
      rotation_old[0] = 1.0;
      rotation_old[1] = 1.0;
      rotation_old[2] = 1.0;
      rotation_old[3] = 0.0;
      rotation_old[4] = 0.0;
      rotation_old[5] = 0.0;
      rotation_old[6] = 0.0;
      rotation_old[7] = 0.0;
      rotation_old[8] = 0.0;
      rotation_old +=9;

      // initialize rotation to identity.
      rotation_new[0] = 1.0;
      rotation_new[1] = 1.0;
      rotation_new[2] = 1.0;
      rotation_new[3] = 0.0;
      rotation_new[4] = 0.0;
      rotation_new[5] = 0.0;
      rotation_new[6] = 0.0;
      rotation_new[7] = 0.0;
      rotation_new[8] = 0.0;
      rotation_new += 9;
    }
  }

  for ( std::vector< stk::mesh::Bucket * >::const_iterator
        k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
    const stk::mesh::Bucket & bucket = **k ;

    const unsigned num_elements = bucket.size();

    double * const sm = stk::mesh::field_data( Shear_Modulus , bucket.begin() );
    double * const dm = stk::mesh::field_data( Dilatational_Modulus , bucket.begin() );
    double * const twomu = stk::mesh::field_data( Material_eff_twomu , bucket.begin() );
    double * const bulk = stk::mesh::field_data( Material_eff_bulk_mod , bucket.begin() );
    double * const mv = stk::mesh::field_data( Midstep_volume , bucket.begin() );
    double * const hg_energy = stk::mesh::field_data( Hourglass_energy , bucket.begin() );
    double * const internal_energy = stk::mesh::field_data( Internal_energy , bucket.begin() );

    for ( unsigned i = 0 ; i < num_elements ; ++i ) {
      sm[i]              = SM;
      dm[i]              = DM;
      twomu[i]           = TM;
      bulk[i]            = BM;
      mv[i]              = 0.0;
      hg_energy[i]       = 0.0;
      internal_energy[i] = 0.0;
    }
  }

  //------------------------------------------------------------------
    time_max[2] = stk::wall_dtime( wtime );
  //------------------------------------------------------------------

  // Scratch space

  enum { num_dof = 24 };
  const int num_dof_max_bucket = num_dof * maximum_bucket_size ;
  std::vector< double > vel(                 num_dof_max_bucket );
  std::vector< double > element_coordinates( num_dof_max_bucket );
  std::vector< double > force_new(           num_dof_max_bucket );

  wtime = stk::wall_time();

  for(unsigned n=0; n<num_trials; ++n) {
    //
    // Call Internal Force!!!
    //

    wtime = stk::wall_time();

    // Need to zero out the old accumulated internal force so that its does not pollute the new accumulation.
    zero_field_data( mesh_bulk_data , stk::mesh::fem::FEMMetaData::NODE_RANK , fint_field );

    for ( std::vector< stk::mesh::Bucket * >::const_iterator
          k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned( **k ) ) {
      const stk::mesh::Bucket & bucket = **k ;

      const int num_elements = bucket.size();
      double *mid_hg_op_ptr = stk::mesh::field_data( MidHourglassOp, bucket.begin());
      double *material_eff_twomu_ptr = stk::mesh::field_data( Material_eff_twomu, bucket.begin());
      double *material_eff_bulk_modulus_ptr = stk::mesh::field_data( Material_eff_bulk_mod, bucket.begin());
      double *mid_step_volume_ptr = stk::mesh::field_data( Midstep_volume, bucket.begin());
      double *element_time_step_ptr = stk::mesh::field_data( Element_time_step, bucket.begin());
      double *element_mass_ptr = stk::mesh::field_data( Element_mass, bucket.begin());
      double *hg_energy_ptr = stk::mesh::field_data( Hourglass_energy, bucket.begin());
      double *internal_energy_ptr = stk::mesh::field_data (Internal_energy, bucket.begin());
      double *shear_modulus_ptr =  stk::mesh::field_data( Shear_Modulus, bucket.begin());
      double *dilatational_modulus_ptr =  stk::mesh::field_data( Dilatational_Modulus, bucket.begin() );

      double *rotation_old_ptr = stk::mesh::field_data( RotationOld, bucket.begin());
      double *rotation_new_ptr = stk::mesh::field_data( RotationNew, bucket.begin());

      double *stretch_ptr      = stk::mesh::field_data( Stretch, bucket.begin());
      double *strain_rate_ptr  = stk::mesh::field_data( StrainRate, bucket.begin());
      double *stress_old_ptr   = stk::mesh::field_data( StressOld, bucket.begin());
      double *stress_new_ptr   = stk::mesh::field_data( StressNew, bucket.begin());
      double *rotated_stress_ptr =  stk::mesh::field_data( RotatedStress, bucket.begin());

      double *hg_resistance_old_ptr =  stk::mesh::field_data( HourglassResistanceOld, bucket.begin());
      double *hg_resistance_new_ptr =  stk::mesh::field_data( HourglassResistanceNew, bucket.begin());

      double *vorticity_ptr    =  stk::mesh::field_data( Vorticity, bucket.begin());

      const double **coord    =
        (const double **) stk::mesh::field_data( coord_gather,    bucket.begin());
      const double **velocity =
        (const double **) stk::mesh::field_data( velocity_gather, bucket.begin());
      double **fint        = stk::mesh::field_data( fint_gather, bucket.begin());

      const int num_elem_nodes = num_elements * 8 ;

      //--------------------------------
      { // Gather nodal data into contiguous arrays.

        // element-node pointer fields
        const double ** field_coord  = coord ;
        const double ** vnodes       = velocity ;

        // scratch space for gather:
        double * elem_coord = & element_coordinates[0] ;
        double * elem_vel   = & vel[0] ;

        for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
          const double * const f_coord = *field_coord ;
          const double * const f_vel   = *vnodes ;

          elem_coord[0] = f_coord[0];
          elem_coord[1] = f_coord[1];
          elem_coord[2] = f_coord[2];

          elem_vel[0] = f_vel[0];
          elem_vel[1] = f_vel[1];
          elem_vel[2] = f_vel[2];

          ++field_coord ;
          ++vnodes ;
          elem_coord += 3 ;
          elem_vel += 3 ;
        }
      }

      //--------------------------------

      materialParameters.dt          = dt;
      materialParameters.nelements   = num_elements ;
      materialParameters.strain_rate = strain_rate_ptr;
      materialParameters.stress_old  = stress_old_ptr;
      materialParameters.stress_new  = stress_new_ptr;

      hex_element.internalForce( num_elements ,
        dt,
        current_stable_time_step,
        element_time_step_ptr,
        *matmodel,
        materialParameters,
        materialProperties,
        & element_coordinates[0],
        & vel[0],
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

      { // Scatter internal force values.
        double ** f = fint ;
        const double * f_new = & force_new[0] ;
        for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
          double * const node_f = *f ;
          node_f[0] += f_new[0];
          node_f[1] += f_new[1];
          node_f[2] += f_new[2];
          ++f ;
          f_new += 3 ;
        }
      }
    }

    time_max[3] += stk::wall_dtime( wtime );

    stk::mesh::parallel_reduce( mesh_bulk_data , stk::mesh::sum(fint_field) );

    time_max[4] += stk::wall_dtime( wtime );

  }//end for(..num_trials...

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
      << "  Box size            = { " << box_size[0] << " , "
                                      << box_size[1] << " , "
                                      << box_size[2] << " }" << std::endl
      << "  Box sides           = { { "
                                    << box_sides[0][0] << " , "
                                    << box_sides[0][1] << " } , { "
                                    << box_sides[1][0] << " , "
                                    << box_sides[1][1] << " } , { "
                                    << box_sides[2][0] << " , "
                                    << box_sides[2][1] << " } }"
                                    << std::endl
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
}

//--------------------------------------------------------------------
//----------------------------------------------------------------------

void use_case_AD_initialize_data(
  stk::mesh::BulkData & mesh ,
  const VectorField & node_model_corrdinates ,
  const VectorField & node_coordinates ,
  const VectorField & node_velocity )
{
  const double dt = 1.0e-2;

  const std::vector<stk::mesh::Bucket*> & buckets = mesh.buckets( stk::mesh::fem::FEMMetaData::NODE_RANK );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {
    stk::mesh::Bucket & bucket = **k ;
    const unsigned length   = bucket.size();

    double * X     = stk::mesh::field_data( node_model_corrdinates , bucket.begin() );
    double * coord = stk::mesh::field_data( node_coordinates , bucket.begin() );
    double * v     = stk::mesh::field_data( node_velocity , bucket.begin() );

    for ( unsigned i = 0 ; i < length ; ++i ) {
      v[0]     = X[0];
      v[1]     = X[1];
      v[2]     = X[2];
      coord[0] = X[0] + dt * v[0];
      coord[1] = X[1] + dt * v[1];
      coord[2] = X[2] + dt * v[2];
      //std::cout << "KHP: Xn= " << coord[0] << ", " << coord[1] << ", " << coord[2] << "\n";
      X += 3;
      v += 3;
      coord += 3;
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

void use_case_AD_generate_mesh(
  stk::mesh::BulkData & mesh ,
  const unsigned N[] ,
  const VectorField & node_coord ,
  const ElementNodePointerField & /*coord_gather */,
  stk::mesh::Part & hex_block ,
  const unsigned shell_side[][2] )
{
  stk::mesh::fem::FEMMetaData &fem = stk::mesh::fem::FEMMetaData::get(mesh);
  const stk::mesh::EntityRank element_rank = fem.element_rank();
  mesh.modification_begin();

  const unsigned parallel_size = mesh.parallel_size();
  const unsigned parallel_rank = mesh.parallel_rank();

  double t = 0 ;
  size_t num_hex = 0 ;
  size_t num_nodes = 0 ;
  size_t num_block = 0 ;
  int error_flag = 0 ;

  try {

    Iogn::GeneratedMesh gmesh( N[0], N[1], N[2], parallel_size, parallel_rank );

    if ( shell_side[0][0] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::MX ); }
    if ( shell_side[0][1] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::PX ); }
    if ( shell_side[1][0] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::MY ); }
    if ( shell_side[1][1] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::PY ); }
    if ( shell_side[2][0] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::MZ ); }
    if ( shell_side[2][1] ) { gmesh.add_shell_block( Iogn::GeneratedMesh::PZ ); }

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
	gmesh.element_map( i, elem_map );

        std::vector<int> elem_conn( num_elem * top_info.second );

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
            << " ERROR, stk::mesh::fem::FEMMetaData::NODE_RANK not found: "
            << node_map[i] << " = node_map[" << i << "]" ;
        throw std::runtime_error( msg.str() );
      }

      double * const data = stk::mesh::field_data( node_coord , *node );
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

  std::vector< unsigned > entity_count ;

  stk::mesh::Selector selector(stk::mesh::MetaData::get(mesh).locally_owned_part());
  stk::mesh::count_entities( selector, mesh , entity_count );

  stk::all_reduce( mesh.parallel() ,
                   stk::ReduceSum<SpatialDim>( & entity_count[0] ) &
                   stk::ReduceMax<1>( & dt ) );

  if ( 0 == parallel_rank ) {
    std::cout << "On " << parallel_size << " processors, meshed "
              << entity_count[ element_rank ] << " elements and "
              << entity_count[ stk::mesh::fem::FEMMetaData::NODE_RANK ] << " nodes "
              << " in " << dt << " sec"
              << std::endl ;
    std::cout.flush();
  }
}

} // namespace app
} // namespace stk

//----------------------------------------------------------------------

