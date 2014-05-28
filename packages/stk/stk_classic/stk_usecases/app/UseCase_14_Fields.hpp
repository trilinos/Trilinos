/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef USECASE_14_FIELDS_HPP
#define USECASE_14_FIELDS_HPP

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <app/UseCase_14_Common.hpp>

namespace stk {
  namespace mesh {
    class Bucket;
  }
  namespace app {

    struct Fields {

      // Nodal vector fields
      CartesianField *model_coordinates;
      CartesianField *coordinates_field;
      CartesianField *velocity_field;
      CartesianField *fint_field;

      // Element vector fields:
      CartesianField *Vorticity;

      // Element scalar fields:
      ScalarField *Shear_Modulus;
      ScalarField *Dilatational_Modulus;
      ScalarField *Material_eff_twomu;
      ScalarField *Material_eff_bulk_mod;
      ScalarField *Midstep_volume;
      ScalarField *Element_time_step;
      ScalarField *Element_mass;
      ScalarField *Hourglass_energy;
      ScalarField *Internal_energy;

      // Element symmetric tensor fields:
      SymmetricTensorField *Stretch;
      SymmetricTensorField *StrainRate;
      SymmetricTensorField *RotatedStress;

      // Element two state symmetric tensor field, on all elements (as two function calls):
      SymmetricTensorField *StressNew;

      // Element two state full tensor field on all elements (as nested function calls):
      FullTensorField *RotationNew;

      //--------------------------------
      // Hourglass fields, these don't have the convenience functions for 'declare and put'

      HourglassArrayField *HourglassResistanceNew;

      HourglassOpField * MidHourglassOp;

      ElementNodePointerField * coord_gather;
      ElementNodePointerField * velocity_gather;

      ElementNodeVectorField * force_new_field;

      HourglassArrayField  * HourglassResistanceOld;
      SymmetricTensorField * StressOld;
      FullTensorField      * RotationOld;
    };

    //--------------------------------------------------------------------
    void use_case_14_declare_fields(Fields &fields, stk::mesh::MetaData &meta_data);
    void use_case_14_initialize_element_fields(
      Fields &fields,
      const stk::mesh::Selector & selector ,
      const std::vector< stk::mesh::Bucket * > & element_buckets,
      double YM, double PR);

    void use_case_14_initialize_nodal_data(stk::mesh::BulkData & mesh ,
                                           const CartesianField & model_coordinates ,
                                           const CartesianField & coordinates ,
                                           const CartesianField & velocity,
                                           double dt);
  }
}
#endif
