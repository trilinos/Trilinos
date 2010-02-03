#ifndef STK_MESH_UNITTEST_GRID_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_GRID_MESH_FIXTURE_HPP

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

class GridMeshFixture
{
public:
  GridMeshFixture(stk::ParallelMachine pm);

  ~GridMeshFixture();

  stk::mesh::MetaData& meta_data() { return *m_meta_data; }
  stk::mesh::BulkData& bulk_data() { return *m_bulk_data; }

  // intentionally exposed to the public
  std::vector<stk::mesh::EntityId> m_node_ids;
  std::vector<stk::mesh::EntityId> m_quad_face_ids;
  std::vector<stk::mesh::EntityId> m_shell_face_ids;

  const std::vector<stk::mesh::Entity*>& get_predefined_closure() const
  {
    return m_closure;
  }

private:
  void generate_grid();

  stk::mesh::MetaData* m_meta_data;
  stk::mesh::BulkData* m_bulk_data;
  stk::mesh::Part* m_quad_part;
  stk::mesh::Part* m_shell_part;
  std::vector<stk::mesh::Entity*> m_closure;
};

#endif

