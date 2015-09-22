/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_io_util_Gears_hpp
#define stk_io_util_Gears_hpp

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

namespace stk_classic {
namespace mesh {
class MetaData;
class BulkData;
}

namespace io {
namespace util {

struct GearFields {

  enum { SpatialDimension = 3 };

  typedef stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>            CartesianField ;
  typedef stk_classic::mesh::Field<double,stk_classic::mesh::Cylindrical>          CylindricalField ;

  CylindricalField & gear_coord ;
  CartesianField   & model_coord ;

  GearFields( stk_classic::mesh::MetaData & S );
  GearFields( stk_classic::mesh::fem::FEMMetaData & S );

private:
  GearFields();
  GearFields( const GearFields & );
  GearFields & operator = ( const GearFields & );
};

class Gear {
public:
  Gear( stk_classic::mesh::fem::FEMMetaData & S ,
        const std::string & name ,
        const GearFields & gear_fields ,
        const double   center[] ,
        const double   rad_min ,
        const double   rad_max ,
        const size_t   rad_num ,
        const double   z_min ,
        const double   z_max ,
        const size_t   z_num ,
        const size_t   angle_num ,
        const int      turn_direction );

  void mesh( stk_classic::mesh::BulkData &M );
  void turn( double turn_angle ) const ;

  stk_classic::mesh::fem::FEMMetaData *m_mesh_fem_meta_data ;
  stk_classic::mesh::MetaData & m_mesh_meta_data ;
  stk_classic::mesh::BulkData * m_mesh ;
  stk_classic::mesh::Part & m_gear ;
  stk_classic::mesh::Part & m_surf ;
  const GearFields::CylindricalField  & m_gear_coord ;
  const GearFields::CartesianField    & m_model_coord ;

private:

  Gear( const Gear & );
  Gear & operator = ( const Gear & );

  double m_center[3] ;
  double m_z_min ;
  double m_z_max ;
  double m_z_inc ;
  double m_rad_min ;
  double m_rad_max ;
  double m_rad_inc ;
  double m_ang_inc ;
  size_t m_rad_num ;
  size_t m_z_num ;
  size_t m_angle_num ;
  int    m_turn_dir ;

  stk_classic::mesh::Entity &create_node( const std::vector<stk_classic::mesh::Part*> &parts ,
                                  stk_classic::mesh::EntityId node_id_base ,
                                  size_t iz ,
                                  size_t ir ,
                                  size_t ia ) const ;
};

}
}
}

#endif

