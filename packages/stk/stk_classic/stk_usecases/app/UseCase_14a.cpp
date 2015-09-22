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

#include <assert.h>
#include <limits.h>

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

#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_algsup/AlgorithmRunner.hpp>

#include <stk_io/MeshReadWriteUtils.hpp>
#include <init/Ionit_Initializer.h>

#include <common/gnu_malloc_hooks.hpp>
#include <app/UseCase_14_Fields.hpp>
#include <app/UseCase_14_Common.hpp>

//----------------------------------------------------------------------
// This file contains the implementation of use-case 14: internal force computation.
// The function 'use_case_14a_driver' below is the equivalent of 'main'.
//----------------------------------------------------------------------

namespace {

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::fem::FEMMetaData& meta_data)
{
  return meta_data.element_rank();
}

inline stk_classic::mesh::EntityRank get_element_rank(const stk_classic::mesh::Part& part)
{
  return get_element_rank(stk_classic::mesh::fem::FEMMetaData::get(part));
}

}

namespace stk_classic {
namespace app {

typedef stk_classic::mesh::Field<double>                                ScalarField ;
typedef stk_classic::mesh::Field<int>                                   ScalarIntField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::Cartesian>          CartesianField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::FullTensor>         FullTensorField ;
typedef stk_classic::mesh::Field<double, stk_classic::mesh::SymmetricTensor>    SymmetricTensorField ;

int g_lockCollision = 0;

enum { SpatialDim = 3 };

CartesianField &
declare_vector_field_on_all_nodes(
  stk_classic::mesh::MetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), stk_classic::mesh::fem::FEMMetaData::NODE_RANK , meta_data.universal_part() , n1 );
}


CartesianField &
declare_vector_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<CartesianField>(s), get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


ScalarField &
declare_scalar_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<ScalarField>(s), get_element_rank(meta_data), meta_data.universal_part() );
}


ScalarIntField &
declare_scalar_int_field_on_all_nodes( stk_classic::mesh::MetaData & meta_data , const std::string & n )
{
  return put_field( meta_data.declare_field<ScalarIntField>(n) , stk_classic::mesh::fem::FEMMetaData::NODE_RANK, meta_data.universal_part() );
}


SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  stk_classic::mesh::fem::FEMMetaData & meta_data , const std::string & s , unsigned n1 )
{
  return stk_classic::mesh::put_field( meta_data.declare_field<SymmetricTensorField>(s) , get_element_rank(meta_data), meta_data.universal_part() , n1 );
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> &
put_field_on_elements( stk_classic::mesh::Field<Type,T1> & f , stk_classic::mesh::Part & p , unsigned n1 )
{
  stk_classic::mesh::put_field( f , get_element_rank(p) , p , n1 );
  return f ;
}


template< typename Type , class T1 , class T2>
stk_classic::mesh::Field<Type,T1,T2> &
put_field_on_elements( stk_classic::mesh::Field<Type,T1,T2> & f , stk_classic::mesh::Part & p , unsigned n1 , unsigned n2 )
{
  stk_classic::mesh::put_field( f , get_element_rank(p) , p , n1 , n2);
  return f ;
}


template< typename Type , class T1 >
stk_classic::mesh::Field<Type,T1> & put_field_on_all_elements( stk_classic::mesh::Field<Type,T1> & f , unsigned n1 )
{
  put_field_on_elements( f , stk_classic::mesh::MetaData::get(f).universal_part() , n1 );
  return f ;
}

// ------------------------------------------------------------------------

typedef stk_classic::mesh::Field<int *,stk_classic::mesh::ElementNode> ElementNodeLockField ;

extern int g_lockCollision;

inline
int
thread_lock(volatile int *addr)
{
//This function is not portable, as evidenced by the following '#if'.
//The test which exercises this is stk_usecases_rtest/app/UseCase_14a. That test has
//been disabled since it was aborting on some supported platforms.

#if (defined __i386__ || defined __x86_64__) && defined __GNUC__
  register int content = 1;

  asm volatile ("xchgl %0, %1\n\t"
                : "=r" (content),
                  "=m" (*addr)
                : "0" (content),
                "m" (*addr));
  return content;
#else
  std::abort();
  return 0;
#endif
}

inline
void
lock(
  int *         addr)
{
  volatile int *x = addr;

  while (thread_lock(addr) != 0) {
    ++g_lockCollision;
    while (*x)
      ;
  }
}

inline
void
unlock(
  int *         addr)
{
  *addr = 0;
}

class Lock
{
public:
  explicit Lock(int *addr)
    : m_addr(addr)
  {
    lock(m_addr);
  }

  ~Lock()
  {
    try {
      unlock(m_addr);
    }
    catch(...) {}
  }

private:
  Lock& operator=(const Lock&);
  Lock(const Lock&);
  int *         m_addr;
};


class MyHexInternalForceAlg2 {
public:
  MyHexInternalForceAlg2(lame::matParams& matParameters, lame::MatProps& matProperties,
                         lame::Material* matModel,
                         stk_classic::mesh::MetaData& md, int max_kernel_size)
    : hex_element(),
      materialParameters(matParameters),
      materialProperties(matProperties),
      matmodel(matModel)
  {
    std::cout << "Assembling MyHexInternalForceAlg2" << std::endl;
    // Nodal vector fields
    model_coordinates     = md.get_field<CartesianField>("model_coordinates");
    coordinates_field     = md.get_field<CartesianField>("coordinates");
    velocity_field        = md.get_field<CartesianField>("velocity");
    fint_field            = md.get_field<CartesianField>("force_internal");

    // Element vector fields:
    Vorticity             = md.get_field<CartesianField>("Vorticity");

    // Element scalar fields:
    Shear_Modulus         = md.get_field<ScalarField>("shear_modulus");
    Dilatational_Modulus  = md.get_field<ScalarField>("dilatational_modulus");
    Material_eff_twomu    = md.get_field<ScalarField>("material_effictive_two_mu");
    Material_eff_bulk_mod = md.get_field<ScalarField>("material_effective_bulk_moduli");
    Midstep_volume        = md.get_field<ScalarField>("mid_step_volume");
    Element_time_step     = md.get_field<ScalarField>("element_time_step");
    Element_mass          = md.get_field<ScalarField>("element_mass");
    Hourglass_energy      = md.get_field<ScalarField>("hourglass_energy");
    Internal_energy       = md.get_field<ScalarField>("internal_energy");

    // Element symmetric tensor fields:
    Stretch               = md.get_field<SymmetricTensorField>("stretch");
    StrainRate            = md.get_field<SymmetricTensorField>("StrainRate");
    RotatedStress         = md.get_field<SymmetricTensorField>("RotatedStress");

    // Element two state symmetric tensor field, on all elements (as two function calls):
    StressNew             = md.get_field<SymmetricTensorField>("Stress");

    // Element two state full tensor field on all elements (as nested function calls):
    RotationNew           = md.get_field<FullTensorField>("Rotation");

    //--------------------------------
    // Hourglass fields, these don't have the convenience functions for 'declare and put'

    HourglassResistanceNew = md.get_field<HourglassArrayField>("HourglassResistance");

    MidHourglassOp = md.get_field<HourglassOpField>("mid_hourglass_operator");

    coord_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("coord_gather");
    velocity_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("velocity_gather");
    fint_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("fint_gather");
    fint_lock_gather = md.get_field<stk_classic::mesh::ElementNodeLockField>("fint_lock_gather");

    force_new_field = md.get_field<ElementNodeVectorField>("force_new_field");

    HourglassResistanceOld = &HourglassResistanceNew->field_of_state( stk_classic::mesh::StateOld );
    StressOld              = &StressNew->field_of_state( stk_classic::mesh::StateOld );
    RotationOld            = &RotationNew->field_of_state( stk_classic::mesh::StateOld );
  }

  enum { maximum_entity_count = 1000 };

  void apply( stk_classic::mesh::Bucket::iterator ibegin ,
	      stk_classic::mesh::Bucket::iterator iend ) const
  {
    enum { nodes_per_elem = 8, num_nodes_X_3 = 24 };

    const int num_elements = iend - ibegin;
    double *mid_hg_op_ptr                 = stk_classic::mesh::field_data( *MidHourglassOp,          ibegin);
    double *material_eff_twomu_ptr        = stk_classic::mesh::field_data( *Material_eff_twomu,      ibegin);
    double *material_eff_bulk_modulus_ptr = stk_classic::mesh::field_data( *Material_eff_bulk_mod,   ibegin);
    double *mid_step_volume_ptr           = stk_classic::mesh::field_data( *Midstep_volume,          ibegin);
    double *element_time_step_ptr         = stk_classic::mesh::field_data( *Element_time_step,       ibegin);
    double *element_mass_ptr              = stk_classic::mesh::field_data( *Element_mass,            ibegin);
    double *hg_energy_ptr                 = stk_classic::mesh::field_data( *Hourglass_energy,        ibegin);
    double *internal_energy_ptr           = stk_classic::mesh::field_data (*Internal_energy,         ibegin);
    double *shear_modulus_ptr             = stk_classic::mesh::field_data( *Shear_Modulus,           ibegin);
    double *dilatational_modulus_ptr      = stk_classic::mesh::field_data( *Dilatational_Modulus,    ibegin );

    double *rotation_old_ptr              = stk_classic::mesh::field_data( *RotationOld,             ibegin);
    double *rotation_new_ptr              = stk_classic::mesh::field_data( *RotationNew,             ibegin);

    double *stretch_ptr      = stk_classic::mesh::field_data( *Stretch,  ibegin);
    double *strain_rate_ptr  = stk_classic::mesh::field_data( *StrainRate,  ibegin);
    double *stress_old_ptr   = stk_classic::mesh::field_data( *StressOld,  ibegin);
    double *stress_new_ptr   = stk_classic::mesh::field_data( *StressNew,  ibegin);
    double *rotated_stress_ptr =  stk_classic::mesh::field_data( *RotatedStress,  ibegin);

    double *hg_resistance_old_ptr =  stk_classic::mesh::field_data( *HourglassResistanceOld,  ibegin);
    double *hg_resistance_new_ptr =  stk_classic::mesh::field_data( *HourglassResistanceNew,  ibegin);

    double *vorticity_ptr    =  stk_classic::mesh::field_data( *Vorticity,  ibegin);

    const double **coord    =
      (const double **) stk_classic::mesh::field_data( *coord_gather,     ibegin);
    const double **velocity =
      (const double **) stk_classic::mesh::field_data( *velocity_gather,  ibegin);
    double* force_new = stk_classic::mesh::field_data( *force_new_field,  ibegin);
   double **fint        = stk_classic::mesh::field_data( *fint_gather,  ibegin);
   int **fint_lock        = stk_classic::mesh::field_data( *fint_lock_gather,  ibegin);

    const int num_elem_nodes = nodes_per_elem * num_elements ;

    //--------------------------------

    double elem_coord[ maximum_entity_count * num_nodes_X_3 ];
    double elem_vel[   maximum_entity_count * num_nodes_X_3 ];

    { // Gather nodal data into contiguous arrays.

      // element-node pointer fields
      const double ** field_coord  = coord ;
      const double ** vnodes       = velocity ;

      double* elem_coord_ptr = &elem_coord[0];
      double* elem_vel_ptr   = &elem_vel[0];

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
    }

    //--------------------------------

    double dt = materialParameters.dt;
    double current_stable_time_step = dt;
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
                               &elem_coord[0],
                               &elem_vel[0],
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
      int ** f_lock = fint_lock ;
      const double * f_new = & force_new[0] ;
      for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
        double * const node_f = *f ;
        {
          int * const node_f_lock = *f_lock ;
          Lock __lock(node_f_lock);

          node_f[0] += f_new[0];
          node_f[1] += f_new[1];
          node_f[2] += f_new[2];
        }
        ++f ;
        ++f_lock ;
        f_new += 3 ;
      }
    }
  }

private:
  MyHexInternalForceAlg2& operator=(const MyHexInternalForceAlg2&);
  MyHexInternalForceAlg2(const MyHexInternalForceAlg2&);

  APSHex8ug hex_element;
  lame::matParams& materialParameters;
  lame::MatProps& materialProperties;
  lame::Material* matmodel;

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
  ElementNodePointerField * fint_gather;
  ElementNodeLockField * fint_lock_gather;

  ElementNodeVectorField * force_new_field;

  HourglassArrayField  * HourglassResistanceOld;
  SymmetricTensorField * StressOld;
  FullTensorField      * RotationOld;
};

class MyHexInternalForceAlg3 {
public:
  MyHexInternalForceAlg3(lame::matParams& matParameters, lame::MatProps& matProperties,
                         lame::Material* matModel,
                         stk_classic::mesh::MetaData& md, int max_bucket_size)
    : m_lock(0),
      hex_element(),
      materialParameters(matParameters),
      materialProperties(matProperties),
      matmodel(matModel)
  {
    std::cout << "Assembling MyHexInternalForceAlg3" << std::endl;

    // Nodal vector fields
    model_coordinates     = md.get_field<CartesianField>("model_coordinates");
    coordinates_field     = md.get_field<CartesianField>("coordinates");
    velocity_field        = md.get_field<CartesianField>("velocity");
    fint_field            = md.get_field<CartesianField>("force_internal");

    // Element vector fields:
    Vorticity             = md.get_field<CartesianField>("Vorticity");

    // Element scalar fields:
    Shear_Modulus         = md.get_field<ScalarField>("shear_modulus");
    Dilatational_Modulus  = md.get_field<ScalarField>("dilatational_modulus");
    Material_eff_twomu    = md.get_field<ScalarField>("material_effictive_two_mu");
    Material_eff_bulk_mod = md.get_field<ScalarField>("material_effective_bulk_moduli");
    Midstep_volume        = md.get_field<ScalarField>("mid_step_volume");
    Element_time_step     = md.get_field<ScalarField>("element_time_step");
    Element_mass          = md.get_field<ScalarField>("element_mass");
    Hourglass_energy      = md.get_field<ScalarField>("hourglass_energy");
    Internal_energy       = md.get_field<ScalarField>("internal_energy");

    // Element symmetric tensor fields:
    Stretch               = md.get_field<SymmetricTensorField>("stretch");
    StrainRate            = md.get_field<SymmetricTensorField>("StrainRate");
    RotatedStress         = md.get_field<SymmetricTensorField>("RotatedStress");

    // Element two state symmetric tensor field, on all elements (as two function calls):
    StressNew             = md.get_field<SymmetricTensorField>("Stress");

    // Element two state full tensor field on all elements (as nested function calls):
    RotationNew           = md.get_field<FullTensorField>("Rotation");

    //--------------------------------
    // Hourglass fields, these don't have the convenience functions for 'declare and put'

    HourglassResistanceNew = md.get_field<HourglassArrayField>("HourglassResistance");

    MidHourglassOp = md.get_field<HourglassOpField>("mid_hourglass_operator");

    coord_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("coord_gather");
    velocity_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("velocity_gather");
    fint_gather = md.get_field<stk_classic::mesh::ElementNodePointerField>("fint_gather");

    force_new_field = md.get_field<ElementNodeVectorField>("force_new_field");

    HourglassResistanceOld = &HourglassResistanceNew->field_of_state( stk_classic::mesh::StateOld );
    StressOld              = &StressNew->field_of_state( stk_classic::mesh::StateOld );
    RotationOld            = &RotationNew->field_of_state( stk_classic::mesh::StateOld );
  }

//  void run(const stk_classic::mesh::Bucket& bucket, int ibegin, int iend )
  void operator()(const stk_classic::mesh::Bucket::iterator & ibegin, const stk_classic::mesh::Bucket::iterator iend )
  {
    enum { nodes_per_elem = 8, num_nodes_X_3 = 24 };

    const int num_elements = iend - ibegin;
    double *mid_hg_op_ptr                 = stk_classic::mesh::field_data( *MidHourglassOp,          ibegin);
    double *material_eff_twomu_ptr        = stk_classic::mesh::field_data( *Material_eff_twomu,      ibegin);
    double *material_eff_bulk_modulus_ptr = stk_classic::mesh::field_data( *Material_eff_bulk_mod,   ibegin);
    double *mid_step_volume_ptr           = stk_classic::mesh::field_data( *Midstep_volume,          ibegin);
    double *element_time_step_ptr         = stk_classic::mesh::field_data( *Element_time_step,       ibegin);
    double *element_mass_ptr              = stk_classic::mesh::field_data( *Element_mass,            ibegin);
    double *hg_energy_ptr                 = stk_classic::mesh::field_data( *Hourglass_energy,        ibegin);
    double *internal_energy_ptr           = stk_classic::mesh::field_data (*Internal_energy,         ibegin);
    double *shear_modulus_ptr             = stk_classic::mesh::field_data( *Shear_Modulus,           ibegin);
    double *dilatational_modulus_ptr      = stk_classic::mesh::field_data( *Dilatational_Modulus,    ibegin );

    double *rotation_old_ptr              = stk_classic::mesh::field_data( *RotationOld,             ibegin);
    double *rotation_new_ptr              = stk_classic::mesh::field_data( *RotationNew,             ibegin);

    double *stretch_ptr      = stk_classic::mesh::field_data( *Stretch,  ibegin);
    double *strain_rate_ptr  = stk_classic::mesh::field_data( *StrainRate,  ibegin);
    double *stress_old_ptr   = stk_classic::mesh::field_data( *StressOld,  ibegin);
    double *stress_new_ptr   = stk_classic::mesh::field_data( *StressNew,  ibegin);
    double *rotated_stress_ptr =  stk_classic::mesh::field_data( *RotatedStress,  ibegin);

    double *hg_resistance_old_ptr =  stk_classic::mesh::field_data( *HourglassResistanceOld,  ibegin);
    double *hg_resistance_new_ptr =  stk_classic::mesh::field_data( *HourglassResistanceNew,  ibegin);

    double *vorticity_ptr    =  stk_classic::mesh::field_data( *Vorticity,  ibegin);

    const double **coord    =
      (const double **) stk_classic::mesh::field_data( *coord_gather,     ibegin);
    const double **velocity =
      (const double **) stk_classic::mesh::field_data( *velocity_gather,  ibegin);
    double* force_new = stk_classic::mesh::field_data( *force_new_field,  ibegin);
    double **fint        = stk_classic::mesh::field_data( *fint_gather,  ibegin);

    const int num_elem_nodes = nodes_per_elem * num_elements ;

    double elem_coord[24000];
    double elem_vel[24000];

    if (iend - ibegin > 1000)
      throw std::runtime_error("Must increase MyHexInternalForceAlg2 TLS size");

    //--------------------------------
    { // Gather nodal data into contiguous arrays.

      // element-node pointer fields
      const double ** field_coord  = coord ;
      const double ** vnodes       = velocity ;

      double* elem_coord_ptr = &elem_coord[0];
      double* elem_vel_ptr   = &elem_vel[0];

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
    }

    //--------------------------------

    double dt = materialParameters.dt;
    double current_stable_time_step = dt;
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
                               &elem_coord[0],
                               &elem_vel[0],
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
      Lock __lock(&m_lock);

      double ** f = fint ;
      const double * f_new = & force_new[0] ;
      for ( int i = 0 ; i < num_elem_nodes ; ++i ) {
        double * const node_f = *f ;
        {
          node_f[0] += f_new[0];
          node_f[1] += f_new[1];
          node_f[2] += f_new[2];
        }
        ++f ;
        f_new += 3 ;
      }
    }
  }

private:
  MyHexInternalForceAlg3& operator=(const MyHexInternalForceAlg3&);
  MyHexInternalForceAlg3(const MyHexInternalForceAlg3&);

  int                   m_lock;

  APSHex8ug hex_element;
  lame::matParams& materialParameters;
  lame::MatProps& materialProperties;
  lame::Material* matmodel;

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
  ElementNodePointerField * fint_gather;

  ElementNodeVectorField * force_new_field;

  HourglassArrayField  * HourglassResistanceOld;
  SymmetricTensorField * StressOld;
  FullTensorField      * RotationOld;
};

bool use_case_14a_driver(MPI_Comm comm,
                        int num_threads,
                        int num_trials,
                        const std::string &working_directory,
                        const std::string &mesh_filename,
                        const std::string &mesh_type,
                        const std::string &thread_runner,
                        int max_entity_per_bucket,
                        bool performance_test)
{
  bool output = !performance_test; // If running for performance measurements, turn off output

  if (stk_classic::parallel_machine_rank(comm) == 0) {
    std::cout << " stk_mesh Use Case #14a - element internal force, begin" << std::endl ;
  }

  std::cout << "Running '" << mesh_filename << "' case, num_trials = "
            << num_trials << std::endl;

  const AlgorithmRunnerInterface* alg_runner = NULL ;
  if ( thread_runner.empty() ||
       thread_runner == std::string("NonThreaded") ) {
    alg_runner = stk_classic::algorithm_runner_non_thread();
  }
  else if ( thread_runner == std::string("TPI") ) {
    alg_runner = stk_classic::algorithm_runner_tpi(num_threads);
  }
  else if ( thread_runner == std::string("TBB") ) {
    alg_runner = stk_classic::algorithm_runner_tbb(num_threads);
  }

  if (alg_runner != NULL) {
    if (stk_classic::parallel_machine_rank(comm) == 0)
      std::cout << "Using " << thread_runner << " algorithm runner, num_threads = " << num_threads << std::endl;
  } else {
    std::cout << "ERROR, failed to obtain requested AlgorithmRunner '" << thread_runner << "'." << std::endl;
    return false;
  }

  double dt  = 1.0e-03;
  double YM  = 1e7; // Young's Modulus
  double PR  = 0.33; // Poisson Ratio
  double element_mass_dummy_value = 1.0;
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
  //   [0] = stk_classic::mesh::MetaData creation
  //   [1] = stk_classic::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force
  //   [4] = Parallel swap-add of internal force

  double time_min[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double time_max[9] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
  double wtime = 0 ;

  //--------------------------------------------------------------------

  reset_malloc_stats();

  if ( 0 == stk_classic::parallel_machine_rank( comm ) ) {
    std::cout << "stk_mesh performance use case #14a" << std::endl
              << "  Number Processes = " << stk_classic::parallel_machine_size( comm )
              << std::endl ;
    std::cout.flush();
  }

  //--------------------------------------------------------------------

  // Initialize IO system.  Registers all element types and storage
  // types and the exodusII default database type.
  Ioss::Init::Initializer init_db;

  {
    wtime = stk_classic::wall_time();

    //------------------------------------------------------------------
    // Declare the mesh meta data: element blocks and associated fields

    stk_classic::mesh::fem::FEMMetaData meta_data( SpatialDim );
    stk_classic::io::MeshData mesh_data;
    std::string filename = working_directory + mesh_filename;
    stk_classic::io::create_input_mesh(mesh_type, filename, comm,
			       meta_data, mesh_data);
    stk_classic::io::define_input_fields(mesh_data, meta_data);

    stk_classic::mesh::Part & universal              = meta_data.universal_part();
    const stk_classic::mesh::EntityRank node_rank    = meta_data.node_rank();
    const stk_classic::mesh::EntityRank element_rank = meta_data.element_rank();
    // Nodal vector fields
    CartesianField &model_coordinates  = stk_classic::mesh::put_field(meta_data.declare_field<CartesianField>("coordinates"), 
                                                                                       node_rank, universal, SpatialDim);
    CartesianField &coordinates_field  = stk_classic::mesh::put_field(meta_data.declare_field<CartesianField>("current_coordinates"), 
                                                                                       node_rank, universal, SpatialDim);
    CartesianField &velocity_field     = stk_classic::mesh::put_field(meta_data.declare_field<CartesianField>("velocity"),    
                                                                                       node_rank, universal, SpatialDim);
    CartesianField &fint_field         = stk_classic::mesh::put_field(meta_data.declare_field<CartesianField>("force_internal"), 
                                                                                       node_rank, universal, SpatialDim);

    ScalarIntField &fint_lock_field    = stk_classic::mesh::put_field(meta_data.declare_field<ScalarIntField>("force_internal_lock"),
                                                                                       node_rank, universal);

    // Element vector fields:
    CartesianField &Vorticity          = stk_classic::mesh::put_field(meta_data.declare_field<CartesianField>("Vorticity"), 
                                                                                       element_rank, universal, SpatialDim);

    // Element scalar fields:
    ScalarField &Shear_Modulus         = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("shear_modulus" ), element_rank, universal);
    ScalarField &Dilatational_Modulus  = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("dilatational_modulus"), element_rank, universal);
    ScalarField &Material_eff_twomu    = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("material_effictive_two_mu"), element_rank, universal);
    ScalarField &Material_eff_bulk_mod = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("material_effective_bulk_moduli"), element_rank, universal);
    ScalarField &Midstep_volume        = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("mid_step_volume"), element_rank, universal);
                                         stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("element_time_step"), element_rank, universal);
    ScalarField &Element_mass          = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("element_mass"), element_rank, universal);
    ScalarField &Hourglass_energy      = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("hourglass_energy"), element_rank, universal);
    ScalarField &Internal_energy       = stk_classic::mesh::put_field(meta_data.declare_field<ScalarField>("internal_energy"), element_rank, universal);

    // Element symmetric tensor fields:
    SymmetricTensorField &Stretch      = declare_symmetric_tensor_field_on_all_elements(meta_data, "stretch",    6 );
    SymmetricTensorField &StrainRate   = declare_symmetric_tensor_field_on_all_elements(meta_data, "StrainRate", 6 );
    SymmetricTensorField &RotatedStress= declare_symmetric_tensor_field_on_all_elements(meta_data, "RotatedStress", 6 );

    //--------------------------------
    // The multi-state fields don't have the 'declare and put' convenience functions (yet)
    //
    //  For clarity declare a integer to used for the number of states arguments.
    //
    const unsigned two_states = 2 ;

    // Element two state symmetric tensor field, on all elements (as two function calls):
    SymmetricTensorField &StressNew = meta_data.declare_field< SymmetricTensorField >( "Stress" , two_states );
    put_field_on_all_elements( StressNew , 6 );

    // Element two state full tensor field on all elements (as nested function calls):
    FullTensorField &RotationNew =
      put_field_on_all_elements( meta_data.declare_field< FullTensorField >("Rotation", two_states ) , 9 );

    //--------------------------------
    // Hourglass fields, these don't have the convenience functions for 'declare and put'

    HourglassArrayField &HourglassResistanceNew =
      put_field_on_all_elements( meta_data.declare_field< HourglassArrayField >("HourglassResistance" , two_states ) , 12 );

    HourglassOpField & MidHourglassOp =
      put_field_on_all_elements( meta_data.declare_field< HourglassOpField >("mid_hourglass_operator") , 32 );

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

    ElementNodePointerField & coord_gather    = declare_element_node_pointer_field(meta_data.get_meta_data(meta_data) , "coord_gather" , coordinates_field );
    ElementNodePointerField & velocity_gather = declare_element_node_pointer_field(meta_data.get_meta_data(meta_data) , "velocity_gather" , velocity_field );
    ElementNodePointerField & fint_gather     = declare_element_node_pointer_field(meta_data.get_meta_data(meta_data) , "fint_gather" , fint_field );
    ElementNodeLockField & fint_lock_gather   = declare_element_node_lock_field(meta_data.get_meta_data(meta_data), "fint_lock_gather" , fint_lock_field );

    //----------------------------------
    // Declare an element field with one value per connected node in
    // which to temporarily store nodal force that is calculated by
    // the internal force algorithm.

    ElementNodeVectorField & force_new_field = meta_data.declare_field<ElementNodeVectorField>( "force_new_field" , 1 /* 1 state */ );

    // All parts of the meta data:
    {
      const mesh::PartVector & all_parts = meta_data.get_parts();
      for (mesh::PartVector::const_iterator i = all_parts.begin();
           i != all_parts.end(); ++i) {
        mesh::Part * const part = *i ;
        if ( part->primary_entity_rank() == get_element_rank(meta_data) ) {
          put_field_on_elements(coord_gather,    *part, shards::Hexahedron<> ::node_count );
          put_field_on_elements(velocity_gather, *part, shards::Hexahedron<> ::node_count );
          put_field_on_elements( fint_gather ,   *part , shards::Hexahedron<> ::node_count );
          put_field_on_elements( fint_lock_gather , *part , shards::Hexahedron<> ::node_count );
          put_field_on_elements(force_new_field, *part, SpatialDim, shards::Hexahedron<>::node_count );
        }
      }
    }

    //--------------------------------
    // Commit (finalize) the meta data.  Is now ready to be used
    // in the creation and management of mesh bulk data.

    meta_data.commit();

    //------------------------------------------------------------------

    time_max[0] = stk_classic::wall_dtime( wtime );

    //------------------------------------------------------------------
    // stk_classic::mesh::BulkData bulk data conforming to the meta data.
    stk_classic::mesh::BulkData bulk_data(meta_data.get_meta_data(meta_data) , comm, max_entity_per_bucket);
    stk_classic::io::populate_bulk_data(bulk_data, mesh_data);

    //------------------------------------------------------------------
    // Create output mesh...  (input filename + ".out14")
    if (output) {
      filename = working_directory + mesh_filename + ".out14a";
      stk_classic::io::create_output_mesh(filename, comm, bulk_data, mesh_data);
      stk_classic::io::define_output_fields(mesh_data, meta_data, true);
    }

    stk_classic::app::use_case_14_initialize_nodal_data(bulk_data ,
						model_coordinates ,
						coordinates_field ,
						velocity_field, dt);

    time_max[1] = stk_classic::wall_dtime( wtime );

    //------------------------------------------------------------------
    // Ready to run the algorithms:
    //------------------------------------------------------------------

    const std::vector< stk_classic::mesh::Bucket * >
      & element_buckets = bulk_data.buckets( get_element_rank(meta_data) );

    stk_classic::mesh::Selector
      select_owned_buckets( stk_classic::mesh::MetaData::get(bulk_data).locally_owned_part() );

    size_t maximum_bucket_size = 0 ;
    size_t minimum_bucket_size = INT_MAX ;

    for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
          k = element_buckets.begin(); k != element_buckets.end() ; ++k )
      if ( select_owned_buckets( **k ) ) {
        const stk_classic::mesh::Bucket & bucket = **k ;
        if ( maximum_bucket_size < bucket.size() ) { maximum_bucket_size = bucket.size(); }
        if ( minimum_bucket_size > bucket.size() ) { minimum_bucket_size = bucket.size(); }
    }

    // Need both the the old and new states of these two-state fields:

    HourglassArrayField  & HourglassResistanceOld = HourglassResistanceNew.field_of_state( stk_classic::mesh::StateOld );
    SymmetricTensorField & StressOld              = StressNew.field_of_state( stk_classic::mesh::StateOld );
    FullTensorField      & RotationOld            = RotationNew.field_of_state( stk_classic::mesh::StateOld );

    for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
            k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned_buckets( **k ) ) {
      const stk_classic::mesh::Bucket & bucket = **k ;

      const int num_elements = bucket.size();

      double * stretch        = stk_classic::mesh::field_data( Stretch,    bucket.begin() );
      double * strain_rate    = stk_classic::mesh::field_data( StrainRate, bucket.begin() );
      double * stress_new     = stk_classic::mesh::field_data( StressNew,  bucket.begin() );
      double * stress_old     = stk_classic::mesh::field_data( StressOld , bucket.begin() );
      double * rotated_stress = stk_classic::mesh::field_data( RotatedStress , bucket.begin() );
      double * rotation_old   = stk_classic::mesh::field_data( RotationOld , bucket.begin() );
      double * rotation_new   = stk_classic::mesh::field_data( RotationNew , bucket.begin() );
      double * mass           = stk_classic::mesh::field_data( Element_mass , bucket.begin() );
      double * hg_old         = stk_classic::mesh::field_data( HourglassResistanceOld , bucket.begin() );
      double * hg_new         = stk_classic::mesh::field_data( HourglassResistanceNew , bucket.begin() );
      double * vorticity_ptr  = stk_classic::mesh::field_data( Vorticity, bucket.begin());
      double * mid_hg_op_ptr  = stk_classic::mesh::field_data( MidHourglassOp, bucket.begin());

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
    if (stk_classic::parallel_machine_rank(comm) == 0)
      std::cout << "KHP: Done element field init\n";

    for ( std::vector< stk_classic::mesh::Bucket * >::const_iterator
            k = element_buckets.begin(); k != element_buckets.end() ; ++k ) if ( select_owned_buckets( **k ) ) {
      const stk_classic::mesh::Bucket & bucket = **k ;

      const unsigned num_elements = bucket.size();

      double * const sm = stk_classic::mesh::field_data( Shear_Modulus , bucket.begin() );
      double * const dm = stk_classic::mesh::field_data( Dilatational_Modulus , bucket.begin() );
      double * const twomu = stk_classic::mesh::field_data( Material_eff_twomu , bucket.begin() );
      double * const bulk = stk_classic::mesh::field_data( Material_eff_bulk_mod , bucket.begin() );
      double * const mv = stk_classic::mesh::field_data( Midstep_volume , bucket.begin() );
      double * const hg_energy = stk_classic::mesh::field_data( Hourglass_energy , bucket.begin() );
      double * const internal_energy = stk_classic::mesh::field_data( Internal_energy , bucket.begin() );

      std::fill(sm,              sm+num_elements,           SM);
      std::fill(dm,              dm+num_elements,           DM);
      std::fill(twomu,           twomu+num_elements,        TM);
      std::fill(bulk,            bulk+num_elements,         BM);
      std::fill(mv,              mv+num_elements,              0.0);
      std::fill(hg_energy,       hg_energy+num_elements,       0.0);
      std::fill(internal_energy, internal_energy+num_elements, 0.0);
    }

    zero_field_data( bulk_data , stk_classic::mesh::fem::FEMMetaData::NODE_RANK , fint_lock_field );

    //------------------------------------------------------------------
    time_max[2] = stk_classic::wall_dtime( wtime );
    //------------------------------------------------------------------

    wtime = stk_classic::wall_time();

    MyHexInternalForceAlg2 elem_alg(materialParameters, materialProperties, matmodel,
                                    meta_data.get_meta_data(meta_data), maximum_bucket_size);
//     MyHexInternalForceAlg3 elem_alg(materialParameters, materialProperties, matmodel,
//                                     meta_data, maximum_bucket_size);

    for(int n=0; n<num_trials; ++n) {
      //
      // Call Internal Force!!!
      //

      wtime = stk_classic::wall_time();

      // Need to zero out the old accumulated internal force so that it
      // does not pollute the new accumulation.
      zero_field_data( bulk_data , stk_classic::mesh::fem::FEMMetaData::NODE_RANK , fint_field );

      mesh::PartVector empty_union_vector;
      alg_runner->run( select_owned_buckets , empty_union_vector, element_buckets , elem_alg );

      time_max[3] += stk_classic::wall_dtime( wtime );

      stk_classic::mesh::parallel_reduce( bulk_data , stk_classic::mesh::sum(fint_field) );

      time_max[4] += stk_classic::wall_dtime( wtime );

      if (output) {
        stk_classic::io::process_output_request(mesh_data, bulk_data, n);
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

  time_max[8] = stk_classic::wall_dtime( wtime );

  time_min[0] = time_max[0] ;
  time_min[1] = time_max[1] ;
  time_min[2] = time_max[2] ;
  time_min[3] = time_max[3] ;
  time_min[4] = time_max[4] ;
  time_min[5] = time_max[5] ;
  time_min[6] = time_max[6] ;
  time_min[7] = time_max[7] ;
  time_min[8] = time_max[8] ;

  stk_classic::all_reduce( comm , stk_classic::ReduceMax<9>( time_max ) & stk_classic::ReduceMin<9>( time_min ) );

  time_max[3] /= num_trials ;
  time_max[4] /= num_trials ;
  time_max[5] /= num_trials ;
  time_max[6] /= num_trials ;

  time_min[3] /= num_trials ;
  time_min[4] /= num_trials ;
  time_min[5] /= num_trials ;
  time_min[6] /= num_trials ;

  //   [0] = stk_classic::mesh::MetaData creation
  //   [1] = stk_classic::mesh::BulkData creation
  //   [2] = Initialization
  //   [3] = Internal force

  if ( ! stk_classic::parallel_machine_rank( comm ) ) {
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
      << "  Lock collisions     = " << g_lockCollision
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
} // namespace stk_classic
