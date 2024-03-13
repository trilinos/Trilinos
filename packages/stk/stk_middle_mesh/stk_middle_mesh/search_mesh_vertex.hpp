#ifndef STK_MIDDLE_MESH_VERTEX_SEARCH_MESH_H
#define STK_MIDDLE_MESH_VERTEX_SEARCH_MESH_H

#include "mesh.hpp"
#include "field.hpp"
#include "predicates/average_normal_field.hpp"
#include "stk_search/Box.hpp"
#include "stk_search/IdentProc.hpp"
#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class SearchMeshVertex
{
  public:
    using Entity = MeshEntity;
    using EntityKey = int;

    using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
    using EntityProcVec = std::vector<EntityProc>;
    using Box = stk::search::Box<double>;
    using BoundingBox = std::pair<Box, EntityProc>;

    //TODO: take the AveragedNormalField as an argument
    explicit SearchMeshVertex(std::shared_ptr<Mesh> mesh, MPI_Comm unionComm, double boxNormalFac=1) :
      m_mesh(mesh),
      m_normalField(predicates::impl::AveragedNormalField(mesh).get_field()),
      m_boxNormalFac(boxNormalFac),
      m_unionComm(unionComm)
    {}


    void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const;


  private:
    using Point = stk::search::Point<double>;

    //TODO: what we really want is a line segment class
    Box create_bounding_box(MeshEntityPtr vert) const;

    std::shared_ptr<Mesh> m_mesh;
    FieldPtr<utils::Point> m_normalField;
    double m_boxNormalFac;
    MPI_Comm m_unionComm;

};


}
}
}
}

#endif