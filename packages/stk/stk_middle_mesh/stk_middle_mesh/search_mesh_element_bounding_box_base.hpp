#ifndef STK_MIDDLE_MESH_SEARCH_MESH_ELEMENT_BOUNDING_BOX_BASE_H
#define STK_MIDDLE_MESH_SEARCH_MESH_ELEMENT_BOUNDING_BOX_BASE_H

#include "mesh.hpp"
#include <stk_search/CoarseSearch.hpp>
#include "stk_search/Box.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {


class SearchMeshElementBoundingBoxBase {
  public:
    using Entity        = mesh::MeshEntity;
    using EntityKey     = int;

    using EntityProc    = stk::search::IdentProc<EntityKey, unsigned>;
    using EntityProcVec = std::vector<EntityProc>;
    using Box           = stk::search::Box<double>;
    using BoundingBox   = std::pair<Box, EntityProc>;

    virtual ~SearchMeshElementBoundingBoxBase() = default;

    virtual void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const = 0;

    virtual std::shared_ptr<mesh::Mesh> get_mesh() const = 0;
};

}
}
}
}

#endif