/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_GEARS_FIXTURE_HPP
#define STK_MESH_FIXTURES_GEARS_FIXTURE_HPP

#include <vector>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/DefaultFEM.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>


namespace stk {
namespace mesh {
namespace fixtures {

class Gear;

class GearsFixture{
 public:

  typedef Field< double , Cylindrical>  CylindricalField ;
  typedef Field< double , Cartesian>    CartesianField ;

  enum { SpatialDimension = 3 };

  GearsFixture( stk::ParallelMachine pm, size_t num_gears);
  ~GearsFixture();

  void generate_mesh();

  const size_t NUM_GEARS;

  MetaData  meta_data;
  DefaultFEM fem_data;
  BulkData  bulk_data;

  const EntityRank element_rank;

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

  void update_displacement_field();

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

