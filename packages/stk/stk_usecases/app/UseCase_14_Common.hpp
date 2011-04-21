/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef USECASE_14_COMMON_HPP
#define USECASE_14_COMMON_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

typedef stk::mesh::Field<double>                                  ScalarField;
typedef stk::mesh::Field<double, stk::mesh::Cartesian>            CartesianField;
typedef stk::mesh::Field<double, stk::mesh::SymmetricTensor>      SymmetricTensorField;
typedef stk::mesh::Field<double, stk::mesh::FullTensor>           FullTensorField;

//--------------------------------------------------------------------

template< class FieldType >
void zero_field_data( stk::mesh::BulkData & mesh ,
		      stk::mesh::EntityRank type , const FieldType & field )
{
  typedef stk::mesh::BucketArray< FieldType > array_type ;

  const std::vector<stk::mesh::Bucket*> & ks = mesh.buckets( type );

  for ( std::vector<stk::mesh::Bucket*>::const_iterator ik = ks.begin() ; ik != ks.end() ; ++ik ) {
    stk::mesh::Bucket & bucket = **ik ;
    array_type data( field , bucket );
    std::fill(data.contiguous_data(), data.contiguous_data() + data.size(), 0 );
  }
}

//--------------------------------------------------------------------

using std::sqrt;

#define INLINE /* inline */

#define TEMPLATE_SCALAR /* template<typename Scalar> */

typedef double Scalar;

//----------------------------------------------------------------------
// Tags for my special kinds of field types,
// which belong in a separate application-shared file...

struct HourglassOp : public shards::ArrayDimTag {
  const char * name() const  ;
  static const HourglassOp & tag();
private:
  HourglassOp() {}
  HourglassOp( const HourglassOp & );
  HourglassOp & operator = ( const HourglassOp & );
};

inline
const HourglassOp & HourglassOp::tag() {
  static const HourglassOp self ;
  return self ;
}

inline
const char * HourglassOp::name() const {
  static const char n[] = "HourglassOp" ;
  return n ;
}

struct HourglassModes : public shards::ArrayDimTag {
  const char * name() const  ;
  static const HourglassModes & tag();
private:
  HourglassModes() {}
  HourglassModes( const HourglassModes & );
  HourglassModes & operator = ( const HourglassModes & );
};

inline
const HourglassModes & HourglassModes::tag()
{ static const HourglassModes self ; return self ; }

inline
const char * HourglassModes::name() const
{ static const char n[] = "HourglassModes" ; return n ; }

//----------------------------------------------------------------------

typedef stk::mesh::Field<double, HourglassOp>     HourglassOpField;
typedef stk::mesh::Field<double, HourglassModes>  HourglassArrayField;
typedef stk::mesh::Field<double, stk::mesh::Cartesian, stk::mesh::ElementNode> ElementNodeVectorField;

enum { SpatialDim = 3 };


// Specification for the aggressive gather pointer-field for elements.

typedef stk::mesh::Field<double*,stk::mesh::ElementNode> ElementNodePointerField ;

namespace lame {
typedef std::map<std::string, std::vector<double> > MatProps;

/**
 *  The struture, matParams, has all of the information that
 *  might be passed back and forth from the application code
 *  and the constitutive model.  This provides a very simple,
 *  easily modified interface between the application code
 *  and the constitutive model.
 *
 *  The following are the more frequently used parameters.
 *
 *  @param nelements    is the number of material points sent
 *                      to LAME for evaluation
 *  @param dt           is the size of the current time step
 *  @param time         is the current solution time
 *  @param dtrnew       is the suggested ratio for the new time step
 *  @param energy_dep   is the total ammount of energy deposited
 *  @param stress_old   are the stress components at \f$t_{n}\f$
 *  @param stress_new   are the stress components at \f$t_{n+1}\f$
 *  @param state_old    are the state variables at \f$t_{n}\f$
 *  @param state_new    are the state variables at \f$t_{n+1}\f$
 *  @param strain_rate  are the components of the current rate
 *                      of deformation
 *  @param left_stretch are the components of the left stretch
 *                      tensor at \f$t_{n+1}\f$
 *  @param rotation     are the components of the rotation
 *                      tensor at \f$t_{n+1}\f$
 *  @param temp_old     is the temperature at \f$t_{n}\f$
 *  @param temp_new     is the temperature at \f$t_{n+1}\f$
 */
struct matParams {
public:
  matParams()
    : nelements           (0),
      nintg               (0),
      dt                  (0.0),
      time                (0.0),
      dtrnew              (1.0),
      energy_dep          (0.0),
      strain_rate         (NULL),
      stress_old          (NULL),
      stress_new          (NULL),
      state_old           (NULL),
      state_new           (NULL),
      temp_old            (NULL),
      temp_new            (NULL),
      left_stretch        (NULL),
      rotation            (NULL),
      entropy             (NULL),
      energy_balance_term (NULL),
      material_properties (NULL),
      tangent_moduli      (NULL),
      ym_old              (NULL),
      ym_new              (NULL),
      nU_new              (NULL),
      nU_old              (NULL),
      bulk_scaling        (NULL),
      shear_scaling       (NULL),
      tangent_flag        (true),
      tangent_comp        (true)
  {}

  int nelements;
  int nintg;
  double dt;
  double time;
  double dtrnew;
  double energy_dep;
  double * strain_rate;
  double * stress_old;
  double * stress_new;
  double * state_old;
  double * state_new;
  double * temp_old;
  double * temp_new;
  double * left_stretch;
  double * rotation;
  double * entropy;
  double * energy_balance_term;
  double * material_properties;
  double * tangent_moduli;
  double * ym_old;
  double * ym_new;
  double * nU_new;
  double * nU_old;
  double * bulk_scaling;
  double * shear_scaling;
  bool tangent_flag;
  bool tangent_comp;
};


class Material
{
public:
  static void reportError( int & int_val, std::string & message ){}

  explicit Material(const MatProps & props)
    : properties() {
    num_material_properties = 0;
    num_state_vars   = 0;
    num_scratch_vars = 0;

    double lambda = getMaterialProperty("LAMBDA",props);
    shrmod = getMaterialProperty("SHEAR_MODULUS",props);
    datmod = lambda + 2.0*shrmod;
    shrmod0 = shrmod;
    datmod0 = datmod;

    double lPlusS = lambda + shrmod;

    if ( lPlusS > 0 ) {
      youngs   = shrmod*(3.0*lambda+2.0*shrmod)/lPlusS;
      poissons = 0.5*lambda/lPlusS;
    } else {
      youngs   = std::numeric_limits<double>::quiet_NaN();
      poissons = std::numeric_limits<double>::quiet_NaN();
    }
  }

  virtual ~Material()
  {}

  int getNumStateVars(){
    return num_state_vars;
  }

  virtual int initialize( matParams * p );

  virtual int getStress( matParams * p ) = 0;

  virtual int loadStepInit( matParams * p ) {
    return 0;
  }

protected:
  std::vector<double> properties;
  double shrmod;
  double datmod;
  double youngs;
  double poissons;
  double shrmod0;
  double datmod0;
  int num_material_properties;
  int num_scratch_vars;
  int num_state_vars;

  double getMaterialProperty( const std::string & name, const MatProps & props);

private:
  Material( const Material & );
  Material & operator= ( const Material & );
};


/**
 * This is the class for the elastic constitutive model.  The
 * elastic model is a hypoelastic model with the following constitutive description
 *
 * \f$
 *     \sigma_{ij}^{n+1} = \sigma_{ij}^{n}
 *     + \Delta t \left( \lambda \delta_{ij} D_{kk}
 *                       + 2\mu D_{ij} \right)
 * \f$
 */
class Elastic: public Material
{
public:
  static Material *createMaterial(const MatProps & props ) {
    return new Elastic(props);
  }

  explicit Elastic( const MatProps & props )
    : Material(props)
  {
    num_material_properties = 2;
    properties.resize(num_material_properties);
    properties[0] = getMaterialProperty("YOUNGS_MODULUS",props);
    properties[1] = getMaterialProperty("POISSONS_RATIO",props);
  }

  ~Elastic() {
  }

  virtual int initialize( matParams * p );

  virtual int getStress( matParams * p );

private:
  Elastic( const Elastic & );
  Elastic & operator= ( const Elastic & );
};

} // namespace lame

namespace APS {
namespace Hex {

typedef int		Int;
typedef double		Real;

TEMPLATE_SCALAR
Int elem_ug3dh8_mi_compute_stretch(const Int nelem,
                                          const Scalar dt,
                                          const Scalar *const cordel,
                                          const Scalar *const vel,
                                          const Scalar *const rotation_old,
                                          Scalar *const mid_vol,
                                          Scalar *const vorticity,
                                          Scalar *const rotation_new,
                                          Scalar *const stretch,
                                          Scalar *const rotated_stretching,
                                          Scalar *const mid_hgop);

} // namespace Hex
} // namespace APS


class APSHex8ug
{
public:
  APSHex8ug()
  {}
  ~APSHex8ug()
  {}

  int num_nodes() const { return 8; }

  int compute_stretch( const int num_elements,
                       const double dt,
                       const double *coordinates,
                       const double *velocity,
                       const double *rotation_old,
                       double *volume,
                       double *vorticity_tensor,
                       double *rotation_new,
                       double *stretch,
                       double *strain_rate,
                       double *mid_hgop,
                       bool debug = false
                       ) const {
    int err = APS::Hex::elem_ug3dh8_mi_compute_stretch(num_elements,
                                                       dt,
                                                       coordinates,
                                                       velocity,
                                                       rotation_old,
                                                       volume,
                                                       vorticity_tensor,
                                                       rotation_new,
                                                       stretch,
                                                       strain_rate,
                                                       mid_hgop);

    return err;
  }

  int  internalForce(  const int num_elements,
                       const double dt,
                       double current_stable_time_step,
                       double element_time_step,
                       lame::Material  &material_model,
                       lame::matParams &materialParameters,
                       lame::MatProps  &materialProperties,
                       std::vector<double> &coordinates,
                       std::vector<double> &velocity,
                       std::vector<double> &rotation_old, std::vector<double> &rotation_new,
                       std::vector<double> &midstep_volume,
                       std::vector<double> &vorticity_tensor,
                       std::vector<double> &stretch,
                       std::vector<double> &strain_rate,
                       std::vector<double> &mid_hgop,
                       std::vector<double> &stress_old, std::vector<double> &stress_new,
                       std::vector<double> &rotated_stress,
                       std::vector<double> &material_eff_bulk_mod,
                       std::vector<double> &material_eff_twomu,
                       std::vector<double> &shrmod,
                       std::vector<double> &dilmod,
                       std::vector<double> &element_mass,
                       std::vector<double> &force_new,
                       std::vector<double> &hourglass_energy,
                       std::vector<double> &internal_energy,
                       std::vector<double> &hg_resistance_old, std::vector<double> &hg_resistance_new
                       ) const;

    int  internalForce(  const int num_elements,
                       const double dt,
                       double current_stable_time_step,
                       double* const element_time_step,
                       lame::Material  &material_model,
                       lame::matParams &materialParameters,
                       lame::MatProps  &materialProperties,
                       double *const coordinates,
                       double *const velocity,
                       double *const rotation_old, double* const rotation_new,
                       double *const midstep_volume,
                       double *const vorticity_tensor,
                       double *const stretch,
                       double *const strain_rate,
                       double *const mid_hgop,
                       double *const stress_old, double *const stress_new,
                       double *const rotated_stress,
                       double *const material_eff_bulk_mod,
                       double *const material_eff_twomu,
                       double *const shrmod,
                       double *const dilmod,
                       double *const element_mass,
                       double *const force_new,
                       double *const hourglass_energy,
                       double *const internal_energy,
                       double *const hg_resistance_old, double* const hg_resistance_new
                       ) const;
 private:
  APSHex8ug& operator=(const APSHex8ug&);
  APSHex8ug(const APSHex8ug&);
};

#endif
