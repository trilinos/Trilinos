/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <app/UseCase_14_Fields.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

namespace {

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::MetaData& meta_data)
{
  const stk_classic::mesh::fem::FEMMetaData &fem =  stk_classic::mesh::fem::FEMMetaData::get(meta_data);
  return fem.element_rank();
}

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::Part& part)
{
  return get_element_rank(stk_classic::mesh::MetaData::get(part));
}

CartesianField &
declare_vector_field_on_all_nodes(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), stk_classic::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() , n1 );
}


CartesianField &
declare_vector_field_on_all_elements(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s),
                               get_element_rank(meta_data),
                               meta_data.universal_part(),
                               n1);
}


ScalarField &
declare_scalar_field_on_all_elements(
  stk_classic::mesh::MetaData & meta_data , const std::string & s )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<ScalarField>(s),
                               get_element_rank(meta_data),
                               meta_data.universal_part() );
}


SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field(meta_data.declare_field<SymmetricTensorField>(s),
                              get_element_rank(meta_data),
                              meta_data.universal_part(),
                              n1);
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> &
put_field_on_elements( stk_classic::mesh::Field<Type,T1> & f , stk_classic::mesh::Part & p , unsigned n1 )
{
  stk_classic::mesh::put_field( f , get_element_rank(p) , p , n1 );
  return f ;
}


template< typename Type , class T1 , class T2 >
stk_classic::mesh::Field<Type,T1,T2> &
put_field_on_elements( stk_classic::mesh::Field<Type,T1,T2> & f , stk_classic::mesh::Part & p , unsigned n1 , unsigned n2 )
{
  stk_classic::mesh::put_field( f , get_element_rank(p) , p , n1 , n2 );
  return f ;
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> & put_field_on_all_elements( stk_classic::mesh::Field<Type,T1> & f , unsigned n1 )
{
  put_field_on_elements( f , stk_classic::mesh::MetaData::get(f).universal_part() , n1 );
  return f ;
}

} // namespace <unnamed>

void stk_classic::app::use_case_14_declare_fields(Fields &fields, stk_classic::mesh::MetaData &meta_data)
{
  // Nodal vector fields
  fields.model_coordinates     = &declare_vector_field_on_all_nodes(meta_data, "coordinates", SpatialDim);
  fields.coordinates_field     = &declare_vector_field_on_all_nodes(meta_data, "current_coordinates", SpatialDim);
  fields.velocity_field        = &declare_vector_field_on_all_nodes(meta_data, "velocity",    SpatialDim);
  fields.fint_field            = &declare_vector_field_on_all_nodes(meta_data, "force_internal", SpatialDim);

  // Element vector fields:
  fields.Vorticity             = &declare_vector_field_on_all_elements(meta_data, "Vorticity" , SpatialDim);

  // Element scalar fields:
  fields.Shear_Modulus         = &declare_scalar_field_on_all_elements(meta_data, "shear_modulus" );
  fields.Dilatational_Modulus  = &declare_scalar_field_on_all_elements(meta_data, "dilatational_modulus");
  fields.Material_eff_twomu    = &declare_scalar_field_on_all_elements(meta_data, "material_effictive_two_mu");
  fields.Material_eff_bulk_mod = &declare_scalar_field_on_all_elements(meta_data, "material_effective_bulk_moduli");
  fields.Midstep_volume        = &declare_scalar_field_on_all_elements(meta_data, "mid_step_volume");
  fields.Element_time_step     = &declare_scalar_field_on_all_elements(meta_data, "element_time_step");
  fields.Element_mass          = &declare_scalar_field_on_all_elements(meta_data, "element_mass");
  fields.Hourglass_energy      = &declare_scalar_field_on_all_elements(meta_data, "hourglass_energy");
  fields.Internal_energy       = &declare_scalar_field_on_all_elements(meta_data, "internal_energy");

  // Element symmetric tensor fields:
  fields.Stretch               = &declare_symmetric_tensor_field_on_all_elements(meta_data, "stretch",    6 );
  fields.StrainRate            = &declare_symmetric_tensor_field_on_all_elements(meta_data, "StrainRate", 6 );
  fields.RotatedStress         = &declare_symmetric_tensor_field_on_all_elements(meta_data, "RotatedStress", 6 );

  //--------------------------------
  // The multi-state fields don't have the 'declare and put' convenience functions (yet)
  //
  //  For clarity declare a integer to used for the number of states arguments.
  //
  const unsigned two_states = 2 ;

  // Element two state symmetric tensor field, on all elements (as two function calls):
  fields.StressNew = &meta_data.declare_field< SymmetricTensorField >( "Stress" , two_states );
  put_field_on_all_elements( *fields.StressNew , 6 );

  // Element two state full tensor field on all elements (as nested function calls):
  fields.RotationNew = &put_field_on_all_elements( meta_data.declare_field< FullTensorField >("Rotation", two_states ) , 9 );

  //--------------------------------
  // Hourglass fields, these don't have the convenience functions for 'declare and put'

  fields.HourglassResistanceNew =
    &put_field_on_all_elements( meta_data.declare_field< HourglassArrayField >("HourglassResistance" , two_states ) , 12 );

  fields.MidHourglassOp =
    &put_field_on_all_elements( meta_data.declare_field< HourglassOpField >("mid_hourglass_operator") , 32 );

  //--------------------------------
  // Declare aggressive "gather" fields which are an array of
  // pointers to the element's nodes' coordinate, velocity, and
  // internal force field data.
  //
  // The declarations specify element fields of the following form:
  //
  //     double * coord_gather   [ nodes_per_element ]
  //     double * velocity_gather[ nodes_per_element ]
  //
  // where
  //
  //     coord_gather[i]    == field_data( coordinates_field , element_node[i] )
  //     velocity_gather[i] == field_data( velocity ,          element_node[i] )
  //
  // The number of nodes per element could vary, so the field is put
  // on each element block with a size of the number of nodes per
  // element in that element block.

  fields.coord_gather    = &declare_element_node_pointer_field( meta_data , "coord_gather" ,    *fields.coordinates_field );
  fields.velocity_gather = &declare_element_node_pointer_field( meta_data , "velocity_gather" , *fields.velocity_field );

  //----------------------------------
  // Declare an element field with one value per connected node in
  // which to temporarily store nodal force that is calculated by
  // the internal force algorithm.

  fields.force_new_field = &meta_data.declare_field<ElementNodeVectorField>( "force_new_field" , 1 /* 1 state */ );

  // All parts of the meta data:
  {
    const mesh::PartVector & all_parts = meta_data.get_parts();
    for (mesh::PartVector::const_iterator i = all_parts.begin();
	 i != all_parts.end(); ++i) {
      mesh::Part * const part = *i ;
      if ( part->primary_entity_rank() == get_element_rank(meta_data) ) {
	put_field_on_elements(*fields.coord_gather,    *part, shards::Hexahedron<> ::node_count );
	put_field_on_elements(*fields.velocity_gather, *part, shards::Hexahedron<> ::node_count );
	put_field_on_elements(*fields.force_new_field, *part, SpatialDim, shards::Hexahedron<>::node_count );
      }
    }
  }
}

void stk_classic::app::use_case_14_initialize_nodal_data(stk_classic::mesh::BulkData & mesh ,
						 const CartesianField & model_coordinates ,
						 const CartesianField & coordinates ,
						 const CartesianField & velocity,
						 double dt)
{
  const std::vector<stk_classic::mesh::Bucket*> & buckets = mesh.buckets( stk_classic::mesh::fem::FEMMetaData::NODE_RANK);

  for ( std::vector<stk_classic::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {
    stk_classic::mesh::Bucket & bucket = **k ;
    const unsigned length   = bucket.size();

    double * X     = stk_classic::mesh::field_data( model_coordinates , bucket.begin() );
    double * coord = stk_classic::mesh::field_data( coordinates , bucket.begin() );
    double * v     = stk_classic::mesh::field_data( velocity , bucket.begin() );

    for ( unsigned i = 0 ; i < length ; ++i ) {
      v[0]     = X[0];
      v[1]     = X[1];
      v[2]     = X[2];
      coord[0] = X[0] + dt * v[0];
      coord[1] = X[1] + dt * v[1];
      coord[2] = X[2] + dt * v[2];
      X += 3;
      v += 3;
      coord += 3;
    }
  }
}

void stk_classic::app::use_case_14_initialize_element_fields(
  Fields &fields,
  const stk_classic::mesh::Selector & selector ,
  const std::vector< stk_classic::mesh::Bucket * > & element_buckets,
  double YM, double PR)
{
  double element_mass_dummy_value = 1.0;

  // Need both the the old and new states of these two-state fields:
  fields.HourglassResistanceOld = &fields.HourglassResistanceNew->field_of_state( stk_classic::mesh::StateOld );
  fields.StressOld              = &fields.StressNew->field_of_state( stk_classic::mesh::StateOld );
  fields.RotationOld            = &fields.RotationNew->field_of_state( stk_classic::mesh::StateOld );

  for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
	  k = element_buckets.begin();
          k != element_buckets.end() ; ++k ) if ( selector( **k ) ) {
    const stk_classic::mesh::Bucket & bucket = **k ;

    const int num_elements = bucket.size();

    double * stretch        = stk_classic::mesh::field_data( *fields.Stretch,    bucket.begin() );
    double * strain_rate    = stk_classic::mesh::field_data( *fields.StrainRate, bucket.begin() );
    double * stress_new     = stk_classic::mesh::field_data( *fields.StressNew,  bucket.begin() );
    double * stress_old     = stk_classic::mesh::field_data( *fields.StressOld , bucket.begin() );
    double * rotated_stress = stk_classic::mesh::field_data( *fields.RotatedStress , bucket.begin() );
    double * rotation_old   = stk_classic::mesh::field_data( *fields.RotationOld , bucket.begin() );
    double * rotation_new   = stk_classic::mesh::field_data( *fields.RotationNew , bucket.begin() );
    double * mass           = stk_classic::mesh::field_data( *fields.Element_mass , bucket.begin() );
    double * hg_old         = stk_classic::mesh::field_data( *fields.HourglassResistanceOld , bucket.begin() );
    double * hg_new         = stk_classic::mesh::field_data( *fields.HourglassResistanceNew , bucket.begin() );
    double * vorticity_ptr  = stk_classic::mesh::field_data( *fields.Vorticity, bucket.begin());
    double * mid_hg_op_ptr  = stk_classic::mesh::field_data( *fields.MidHourglassOp, bucket.begin());

    std::fill(strain_rate,    strain_rate    +  6*num_elements, 0.0);
    std::fill(stress_new,     stress_new     +  6*num_elements, 0.0);
    std::fill(stress_old,     stress_old     +  6*num_elements, 0.0);
    std::fill(rotated_stress, rotated_stress +  6*num_elements, 0.0);
    std::fill(mass,           mass           +    num_elements, element_mass_dummy_value);
    std::fill(hg_old,         hg_old         + 12*num_elements, 0.0);
    std::fill(hg_new,         hg_new         + 12*num_elements, 0.0);
    std::fill(vorticity_ptr,  vorticity_ptr  +  3*num_elements, 0.0);
    std::fill(mid_hg_op_ptr,  mid_hg_op_ptr  + 32*num_elements, 0.0);

    for ( int i = 0 ; i < num_elements ; ++i ) {
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

  //----------------------------------
  // Compute other material constants.
  double BM = YM/(3.0-6.0*PR); // Bulk Modulus
  double SM = YM/(2.0+2.0*PR); // Shear Modulus
  double TM = 2.0 * SM;
  double LAMBDA = YM*PR/((1.0+PR)*(1.0-2.0*PR));
  double DM = TM + LAMBDA; // Dilatational Modulus

  for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
	  k = element_buckets.begin();
          k != element_buckets.end() ; ++k ) if ( selector( **k ) ) {
    const stk_classic::mesh::Bucket & bucket = **k ;

    const unsigned num_elements = bucket.size();

    double * const sm              = stk_classic::mesh::field_data( *fields.Shear_Modulus , bucket.begin() );
    double * const dm              = stk_classic::mesh::field_data( *fields.Dilatational_Modulus , bucket.begin() );
    double * const twomu           = stk_classic::mesh::field_data( *fields.Material_eff_twomu , bucket.begin() );
    double * const bulk            = stk_classic::mesh::field_data( *fields.Material_eff_bulk_mod , bucket.begin() );
    double * const mv              = stk_classic::mesh::field_data( *fields.Midstep_volume , bucket.begin() );
    double * const hg_energy       = stk_classic::mesh::field_data( *fields.Hourglass_energy , bucket.begin() );
    double * const internal_energy = stk_classic::mesh::field_data( *fields.Internal_energy , bucket.begin() );

    std::fill(sm,              sm+num_elements,           SM);
    std::fill(dm,              dm+num_elements,           DM);
    std::fill(twomu,           twomu+num_elements,        TM);
    std::fill(bulk,            bulk+num_elements,         BM);
    std::fill(mv,              mv+num_elements,              0.0);
    std::fill(hg_energy,       hg_energy+num_elements,       0.0);
    std::fill(internal_energy, internal_energy+num_elements, 0.0);
  }
}

