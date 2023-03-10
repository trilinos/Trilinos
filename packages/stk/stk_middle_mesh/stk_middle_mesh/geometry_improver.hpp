#ifndef GEOMETRY_IMPROVER_H
#define GEOMETRY_IMPROVER_H

#include "mesh.hpp"
#include "mesh_relational_data.hpp"

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

// base class for classes that improve the geometry somehow (generally
// by using data from both meshes to produce a mesh_in that has lower
// geometric approximation error)
class GeometryImprover
{
  public:
    GeometryImprover(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                     std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<MeshRelationalData> relationalData)
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_meshIn(meshIn)
      , m_relationalData(relationalData)
    {}

    virtual ~GeometryImprover() {}

    virtual void run() = 0;

  protected:
    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_relationalData;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif