#ifndef MIDDLE_GRID_TRIANGULATOR_H
#define MIDDLE_GRID_TRIANGULATOR_H

#include "mesh_relational_data.hpp"
#include "predicates/point_classifier_normal_wrapper.hpp"
#include "bounding_box_search_opts.hpp"

#include <map>
#include <set>

namespace stk {
namespace middle_mesh {
namespace nonconformal4 {
namespace impl {

class MiddleGridTriangulator
{
  public:
    MiddleGridTriangulator(std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
                           std::shared_ptr<mesh::Mesh> meshIn, std::shared_ptr<MeshRelationalData> relationalData,
                           std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> pointClassifier,
                           const search::BoundingBoxSearchOpts& searchOpts = search::BoundingBoxSearchOpts())
      : m_mesh1(mesh1)
      , m_mesh2(mesh2)
      , m_meshIn(meshIn)
      , m_relationalData(relationalData)
      , m_classifier(pointClassifier)
      , m_searchOpts(searchOpts)
    {
      assert(stk::middle_mesh::mesh::count_valid(m_meshIn->get_edges()) > 0);
    }

    void triangulate();

  private:
    /*
        //TODO: DEBUGGING
        int getFakeVertId(mesh::MeshEntityPtr vert_in)
        {
          for (size_t i=0; i < m_relational_data->fake_verts_to_verts_in.size(); ++i)
            if (vert_in == m_relational_data->fake_verts_to_verts_in[i])
              return i;

          throw std::runtime_error("unable to find vert_in fakevert id");
        }
    */

    mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_to_mesh2_element_maps();

    std::shared_ptr<mesh::Mesh> m_mesh1;
    std::shared_ptr<mesh::Mesh> m_mesh2;
    std::shared_ptr<mesh::Mesh> m_meshIn;
    std::shared_ptr<MeshRelationalData> m_relationalData;
    std::shared_ptr<predicates::impl::PointClassifierNormalWrapper> m_classifier;
    search::BoundingBoxSearchOpts m_searchOpts;
    const bool m_output = false;
};

} // namespace impl

} // namespace nonconformal4
} // namespace middle_mesh
} // namespace stk
#endif
