/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_GRID_MESH_FIXTURE_HPP

#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_mesh/fixtures/FixtureNodeSharing.hpp>

namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {
namespace fixtures {

class GridFixture
{
public:
  GridFixture(stk::ParallelMachine pm);

  ~GridFixture();

  MetaData& fem_meta() { return m_fem_meta; }
  BulkData& bulk_data() { return m_bulk_data; }

  Part* quad_part() const { return & m_quad_part; }
  Part* dead_part() const { return & m_dead_part; }

  void generate_grid();

  const unsigned m_spatial_dimension;

  MetaData  m_fem_meta;
  BulkData  m_bulk_data;
  Part &    m_quad_part;
  Part &    m_dead_part;

private:
  NodeToProcsMMap m_nodes_to_procs;

  void fill_node_map(unsigned num_nodes, unsigned num_quad_faces, int p_rank);
};

} // fixtures
} // mesh
} // stk

#endif

