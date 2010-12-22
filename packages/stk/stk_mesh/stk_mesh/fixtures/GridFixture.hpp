/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/DefaultFEM.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

class GridFixture
{
public:
  GridFixture(stk::ParallelMachine pm);

  ~GridFixture();

  MetaData& meta_data() { return m_meta_data; }
  BulkData& bulk_data() { return m_bulk_data; }
  DefaultFEM& fem() { return m_fem; }

  Part* quad_part() const { return & m_quad_part; }
  Part* dead_part() const { return & m_dead_part; }

  void generate_grid();

  const unsigned m_spatial_dimension;

  MetaData      m_meta_data;
  BulkData      m_bulk_data;
  DefaultFEM    m_fem;
  Part &        m_quad_part;
  Part &        m_dead_part;
};

} // fixtures
} // mesh
} // stk

#endif

