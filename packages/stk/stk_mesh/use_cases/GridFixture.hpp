#ifndef STK_MESH_UNITTEST_GRID_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_GRID_MESH_FIXTURE_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

class GridFixture
{
public:
  GridFixture(stk::ParallelMachine pm);

  ~GridFixture();

  stk::mesh::MetaData& meta_data() { return *m_meta_data; }
  stk::mesh::BulkData& bulk_data() { return *m_bulk_data; }

  stk::mesh::Part* quad_part() const { return m_quad_part; }
  stk::mesh::Part* boundary_part() const { return m_boundary_part; }
  stk::mesh::Part* dead_part() const { return m_dead_part; }
  
private:
  void generate_grid();

  stk::mesh::MetaData* m_meta_data;
  stk::mesh::BulkData* m_bulk_data;
  stk::mesh::Part* m_quad_part;
  stk::mesh::Part* m_boundary_part;
  stk::mesh::Part* m_dead_part;
};

#endif

