#ifndef PREDICATES_TRIANGLE_COORD_UTILS
#define PREDICATES_TRIANGLE_COORD_UTILS

#include "intersection_common.hpp"

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class TriangleCoordUtils
{
  public:
    utils::Point compute_xyz_coords(const PointRecordForTriangle& record, bool allowExterior=false)
    {
      if (!allowExterior)
      {
        assert(record.type != PointClassification::Exterior);
      }
      assert(record.el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);

      if (record.type == PointClassification::Vert)
      {
        std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> verts;
        get_downward(record.el, 0, verts.data());
        return verts[record.id]->get_point_orig(0);
      } else
      {
        // TODO: for edges, it would be more consistent with how mesh1 views the world
        //       to compute the edge xi and compute the xyz coordinates from that.  This
        //       would ensure the point is *exactly* on the line instead of approximately
        //       on the line
        return compute_tri_coords_from_xi_3d(record.el, record.m_ptXi);
      }
    }

    utils::Point compute_xi_coords(const PointRecordForTriangle& record, bool allowExterior=false)
    {
      if (!allowExterior)
      {
        assert(record.type != PointClassification::Exterior);
      }
      assert(record.el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);
      return record.m_ptXi;
    }

    double get_edge_xi(const PointRecordForTriangle& record)
    {
      assert(record.type == PointClassification::Edge);
      auto el     = record.el;
      auto type   = el->get_type();
      auto orient = el->get_down_orientation(record.id);
      return get_edge_xi(type, record.id, record.m_ptXi, orient);
    }

    double get_edge_xi(mesh::MeshEntityType type, const int id, const utils::Point& ptXi,
                       mesh::EntityOrientation orient)
    {
      double xi = 0;
      if (type == stk::middle_mesh::mesh::MeshEntityType::Quad)
      {
        switch (id)
        {
          case 0: {
            xi = ptXi.x;
            break;
          }
          case 1: {
            xi = ptXi.y;
            break;
          }
          case 2: {
            xi = 1 - ptXi.x;
            break;
          }
          case 3: {
            xi = 1 - ptXi.y;
            break;
          }
          default:
            throw std::invalid_argument("invalid id value");
        }
      } else if (type == stk::middle_mesh::mesh::MeshEntityType::Triangle)
      {
        double xi3 = 1 - ptXi.x - ptXi.y;
        switch (id)
        {
          case 0: {
            xi = ptXi.x;
            break;
          }
          case 1: {
            xi = ptXi.y;
            break;
          }
          case 2: {
            xi = xi3;
            break;
          }
          default:
            throw std::invalid_argument("invalid id value");
        }
      } else
      {
        throw std::invalid_argument("el must be either a triangle or quad");
      }

      if (orient == mesh::EntityOrientation::Reversed)
        xi = 1 - xi;

      return xi;
    }

    double compute_orthogonal_dist(const PointRecordForTriangle& record1, const int id)
    {
      return compute_orthogonal_dist(record1.el->get_type(), record1.m_ptXi, id);
    }

    double compute_orthogonal_dist(mesh::MeshEntityType type, const utils::Point& ptXi, const int id)
    {
      if (type == stk::middle_mesh::mesh::MeshEntityType::Triangle)
        return compute_orthogonal_dist_triangle(ptXi, id);
      else if (type == stk::middle_mesh::mesh::MeshEntityType::Quad)
        return compute_orthogonal_dist_quad(ptXi, id);
      else
        throw std::invalid_argument("type must be either Quad or Triangle");
    }

    double compute_orthogonal_dist_quad(const utils::Point& ptXi, const int id)
    {
      switch (id)
      {
        case 0: {
          return std::abs(ptXi.y);
        }
        case 1: {
          return std::abs(1 - ptXi.x);
        }
        case 2: {
          return std::abs(1 - ptXi.y);
        }
        case 3: {
          return std::abs(ptXi.x);
        }
        default:
          throw std::invalid_argument("invalid id");
      }
    }

    double compute_orthogonal_dist_triangle(const utils::Point& ptXi, const int id)
    {
      switch (id)
      {
        case 0: {
          return std::abs(ptXi.y);
        }
        case 1: {
          utils::Point r(-1, 1);
          auto x     = ptXi - utils::Point(1, 0);
          double fac = dot(r, x) / dot(r, r);
          utils::Point xParallel(fac * r.x, fac * r.y);
          auto xPerp = x - xParallel;
          return std::sqrt(dot(xPerp, xPerp));
        }
        case 2: {
          return std::abs(ptXi.x);
        }
        default:
          throw std::invalid_argument("invalid id");
      }
    }

    utils::Point orthogonal_proj(mesh::MeshEntityType type, const utils::Point& ptXi, const int id)
    {
      if (type == stk::middle_mesh::mesh::MeshEntityType::Triangle)
        return orthogonal_proj_triangle(ptXi, id);
      else if (type == stk::middle_mesh::mesh::MeshEntityType::Quad)
        return orthogonal_proj_quad(ptXi, id);
      else
        throw std::invalid_argument("type must be either Quad or Triangle");
    }

    utils::Point orthogonal_proj_quad(const utils::Point& ptXi, const int id)
    {
      switch (id)
      {
        case 0: {
          return utils::Point(ptXi.x, 0);
        }
        case 1: {
          return utils::Point(1, ptXi.y);
        }
        case 2: {
          return utils::Point(ptXi.x, 1);
        }
        case 3: {
          return utils::Point(0, ptXi.y);
        }
        default:
          throw std::invalid_argument("invalid id");
      }
    }

    utils::Point orthogonal_proj_triangle(const utils::Point& ptXi, const int id)
    {
      // for the unit right triangle referene element, xi1 = x, xi2 = y, so
      // for those cases the orthogonal projection is easy
      switch (id)
      {
        case 0: {
          return utils::Point(ptXi.x, 0);
        }
        case 1: {
          utils::Point r(-1, 1);
          auto x     = ptXi - utils::Point(1, 0);
          double fac = dot(r, x) / dot(r, r);
          return utils::Point(fac * r.x, fac * r.y);
        }
        case 2: {
          return utils::Point(0, ptXi.y);
        }
        default:
          throw std::invalid_argument("invalid id");
      }
    }

    PointRecordForTriangle create_record(mesh::MeshEntityPtr tri, int edgeId, double edgeXi)
    {
      assert(tri->get_type() == mesh::MeshEntityType::Triangle);
      assert(edgeId >= 0 && edgeId < 3);
      
      utils::Point ptXi;
      if (edgeId == 0)
      {
        ptXi = utils::Point(edgeXi, 0);
      } else if (edgeId == 1)
      {
        ptXi = utils::Point(1 - edgeXi, edgeXi);
      } else // edgeId == 2
      {
        ptXi = utils::Point(0, 1 - edgeXi);
      }

      return PointRecordForTriangle(PointClassification::Edge, edgeId, tri, ptXi);
    }

    PointRecordForTriangle create_record(mesh::MeshEntityPtr tri, int vertId)
    {
      assert(tri->get_type() == mesh::MeshEntityType::Triangle);
      assert(vertId >= 0 && vertId < 3);

      utils::Point ptXi;
      switch (vertId)
      {
        case 0: { ptXi = utils::Point(0, 0); break; }
        case 1: { ptXi = utils::Point(1, 0); break; }
        case 2: { ptXi = utils::Point(0, 1); break; }
        default:
          throw std::runtime_error("incorrect vertexId");
      }

      return PointRecordForTriangle(PointClassification::Vert, vertId, tri, ptXi);
    }

    // given ptXi in the coordinate system of triFrom, converts it to the coordate system of triTo
    PointRecordForTriangle classify_onto(const PointRecordForTriangle& record, mesh::MeshEntityPtr el)
    {
      assert(record.type == PointClassification::Vert || record.type == PointClassification::Edge);
      assert(el->get_type() == mesh::MeshEntityType::Triangle);

      if (record.el == el)
        return record;

      int dim = record.type == PointClassification::Vert ? 0 : 1;

      mesh::MeshEntityPtr entityOld = get_entity(record);

      std::array<mesh::MeshEntityPtr, mesh::MAX_DOWN> entitiesNew;
      mesh::get_downward(el, dim, entitiesNew.data());
      int entityNewLocalId = -1;
      for (int i=0; i < 3; ++i)
        if (entitiesNew[i] == entityOld)
        {
          entityNewLocalId = i;
          break;
        }

      if (record.type == PointClassification::Vert)
      {
        if (entityNewLocalId == -1)
          throw std::runtime_error("the specified vertex is not shared between record.el and el");

        return create_record(el, entityNewLocalId);
      } else // record.type == Edge
      {
        if (entityNewLocalId == -1)
          throw std::runtime_error("the specified edge is not shared between record.el and el");

        double edgeXiOld;
        if (record.id == 0)
          edgeXiOld = record.m_ptXi.x;
        else if (record.id == 1)
          edgeXiOld = 1 - record.m_ptXi.x;
        else
          edgeXiOld = 1 - record.m_ptXi.y;

        double edgeXiNew = 1 - edgeXiOld;

        return create_record(el, entityNewLocalId, edgeXiNew);
      }
    }

    // returns a measure of how far outside the element the point is
    double compute_exterior_deviation(const PointRecordForTriangle& r)
    {
      assert(r.type == PointClassification::Exterior);

      std::array<double, 3> xi = {1 - r.m_ptXi.x - r.m_ptXi.y, r.m_ptXi.x, r.m_ptXi.y};
      double dev = 0.0;
      for (int i=0; i < 3; ++i)
      {
        if (xi[i] > 1)
          dev += std::max(xi[i] - 1, 0.0);
        else if (xi[i] < 0)
          dev += -xi[i];
      }

      return dev;
    }
};

} // namespace impl
} // namespace predicates

} // namespace middle_mesh
} // namespace stk
#endif