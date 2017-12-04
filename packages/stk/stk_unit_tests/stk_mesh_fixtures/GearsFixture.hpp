// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#ifndef STK_MESH_FIXTURES_GEARS_FIXTURE_HPP
#define STK_MESH_FIXTURES_GEARS_FIXTURE_HPP

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/CoordinateSystems.hpp>  // for Cartesian
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <vector>                       // for vector
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { namespace fixtures { class Gear; } } }




namespace stk {
namespace mesh {
namespace fixtures {


struct GearParams {
  GearParams()
    : element_size(0.1),
      radius_min(0.6),
      radius_max(1.05),
      height_min(-0.4),
      height_max(0.4)
  {}

  GearParams(double input_element_size,
             double input_radius_min,
             double input_radius_max,
             double input_height_min,
             double input_height_max)
    : element_size(input_element_size),
      radius_min(input_radius_min),
      radius_max(input_radius_max),
      height_min(input_height_min),
      height_max(input_height_max)
  {}

  double element_size;
  double radius_min;
  double radius_max;
  double height_min;
  double height_max;
};

class GearsFixture{
 public:

  typedef Field< double , Cylindrical>  CylindricalField ;
  typedef Field< double , Cartesian>    CartesianField ;

  enum { SpatialDimension = 3 };

  GearsFixture( stk::ParallelMachine pm, size_t num_gears,
                GearParams gear_params=GearParams());
  ~GearsFixture();

  void generate_mesh();

  const size_t NUM_GEARS;

  MetaData meta_data;
  BulkData bulk_data;

  Part & cylindrical_coord_part;
  Part & hex_part;
  Part & wedge_part;

  CartesianField    & cartesian_coord_field ;
  CartesianField    & displacement_field ;
  CartesianField    & translation_field ;
  CylindricalField  & cylindrical_coord_field;

  Gear & get_gear(size_t i) {
    return * m_gears[i];
  }

  const Gear & get_gear(size_t i) const{
    return * m_gears[i];
  }

  void communicate_model_fields();

 private:

  std::vector<Gear *> m_gears;

  GearsFixture( const GearsFixture & );
  GearsFixture & operator = ( const GearsFixture & );
};

/// \brief Distribute gears across processors
void distribute_gear_across_processors(Gear & gear, GearsFixture::CylindricalField & cylindrical_coord_field);
unsigned destination_processor(const Gear & gear, double rad, double angle, double height, unsigned p_rank, unsigned p_size);

} // fixtures
} // mesh
} // stk

#endif //STK_MESH_FIXTURES_GEARS_FIXTURE_HPP

