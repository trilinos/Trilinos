#ifndef STK_MIDDLE_MESH_SEARCH_MESH_BOUNDING_BOX_NORMAL_H
#define STK_MIDDLE_MESH_SEARCH_MESH_BOUNDING_BOX_NORMAL_H

#include "search_mesh_element_bounding_box_base.hpp"
#include "mesh.hpp"
#include "field.hpp"
#include <stk_search/CoarseSearch.hpp>
#include "stk_search/Box.hpp"

#include <limits>
#include <map>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class SearchMeshElementBoundingBoxNormal : public SearchMeshElementBoundingBoxBase {
  public:

    SearchMeshElementBoundingBoxNormal(std::shared_ptr<mesh::Mesh> inputMesh, MPI_Comm unionComm, mesh::FieldPtr<utils::Point> averagedNormalField, double normalFac = 1)
      : m_mesh(inputMesh)
      , m_unionComm(unionComm)
      , m_averagedNormalField(averagedNormalField)
      , m_normalFac(normalFac)
      {}

    void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const override;

    std::shared_ptr<mesh::Mesh> get_mesh() const override { return m_mesh; }

  private:
    using Point         = stk::search::Point<double>;

    Box create_bounding_box(mesh::MeshEntityPtr element) const;

    std::shared_ptr<mesh::Mesh> m_mesh;
    MPI_Comm m_unionComm;
    mesh::FieldPtr<utils::Point> m_averagedNormalField;
    double m_normalFac;
};

}
}
}
}

#endif