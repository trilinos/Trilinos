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
    using EntityKey = int64_t;

    using EntityProc = stk::search::IdentProc<EntityKey, unsigned>;
    using EntityProcVec = std::vector<EntityProc>;
    using Box = stk::search::Box<double>;
    using BoundingBox = std::pair<Box, EntityProc>;

    explicit SearchMeshVertex(std::shared_ptr<Mesh> mesh, MPI_Comm unionComm, double boxNormalFac=1) :
      m_mesh(mesh),
      m_normalField(predicates::impl::AveragedNormalField(mesh).get_field()),
      m_boxNormalFac(boxNormalFac),
      m_unionComm(unionComm)
    {}


    void fill_bounding_boxes(std::vector<BoundingBox>& boundingBoxes) const
    {
      int proc = utils::impl::comm_rank(m_unionComm); // Local search only
      boundingBoxes.reserve(m_mesh->get_vertices().size());
      for (auto& vert : m_mesh->get_vertices())
        if (vert)
        {
          Box box = create_bounding_box(vert);
          EntityProc entityProc(vert->get_id(), proc);
          boundingBoxes.push_back(BoundingBox(box, entityProc));
        }
    }


  private:
    using Point = stk::search::Point<double>;

    //TODO: what we really want is a line segment class
    Box create_bounding_box(MeshEntityPtr vert) const 
    {
      utils::Point normal = (*m_normalField)(vert, 0, 0);
      utils::Point pt1 =  normal * m_boxNormalFac + vert->get_point_orig(0);
      utils::Point pt2 = -normal * m_boxNormalFac + vert->get_point_orig(0);

      Point minCorner, maxCorner;
      for (int d=0; d < 3; ++d)
      {
        minCorner[d] = std::min(pt1[d], pt2[d]);
        maxCorner[d] = std::max(pt1[d], pt2[d]);
      }

      return Box(minCorner, maxCorner);
    }


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