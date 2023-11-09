#ifndef FIELD_BASE_H
#define FIELD_BASE_H

#include "mesh_entity.hpp"
#include <array>
#include <memory>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {

// defines number of nodes on each *dimension* entity
struct FieldShape
{
    explicit FieldShape(const int vert = 0, const int edge = 0, const int el = 0)
      : count{vert, edge, el}
    {}

    std::array<int, 3> count;

    int get_num_nodes(int dim) const { return count. at(dim); }
};

inline bool operator==(const FieldShape& lhs, const FieldShape& rhs)
{
  return lhs.count == rhs.count;
}

inline bool operator!=(const FieldShape& lhs, const FieldShape& rhs)
{
  return !(lhs == rhs);
}

namespace impl {


using FieldShape = FieldShape;
using EntityCount = FieldShape;

class FieldManager; // forward declaration

class FieldBase
{
  public:
    virtual ~FieldBase(){};

  protected:
    // function to be called when a new entity is created
    virtual void add_entity(int dim) = 0;

    // function to be called when mesh is condensing its array
    virtual void condense_arrays(const std::vector<::stk::middle_mesh::mesh::MeshEntityPtr>& verts,
                                 const std::vector<::stk::middle_mesh::mesh::MeshEntityPtr>& edges,
                                 const std::vector<::stk::middle_mesh::mesh::MeshEntityPtr>& elements) = 0;

    friend FieldManager;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
