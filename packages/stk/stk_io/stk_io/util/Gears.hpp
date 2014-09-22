// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_io_util_Gears_hpp
#define stk_io_util_Gears_hpp

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

namespace stk {
namespace mesh {
class MetaData;
class BulkData;
}

namespace io {
namespace util {

struct GearFields {

  enum { SpatialDimension = 3 };

  typedef stk::mesh::Field<double,stk::mesh::Cartesian>            CartesianField ;
  typedef stk::mesh::Field<double,stk::mesh::Cylindrical>          CylindricalField ;

  CylindricalField & gear_coord ;
  CartesianField   & model_coord ;

  GearFields( stk::mesh::MetaData & S );

private:
  GearFields();
  GearFields( const GearFields & );
  GearFields & operator = ( const GearFields & );
};

class Gear {
public:
  Gear( stk::mesh::MetaData & S ,
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

  void mesh( stk::mesh::BulkData &M );
  void turn( double turn_angle ) const ;

  stk::mesh::MetaData & m_mesh_meta_data ;
  stk::mesh::BulkData * m_mesh ;
  stk::mesh::Part & m_gear ;
  stk::mesh::Part & m_surf ;
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

  stk::mesh::Entity create_node( const std::vector<stk::mesh::Part*> &parts ,
                                  stk::mesh::EntityId node_id_base ,
                                  size_t iz ,
                                  size_t ir ,
                                  size_t ia ) const ;
};

}
}
}

#endif

