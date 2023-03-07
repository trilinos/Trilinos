#ifndef ELEMENT_OPERATIONS_2D
#define ELEMENT_OPERATIONS_2D

#include "mesh.hpp"
#include "mesh_entity.hpp"
#include "newton2.hpp"
#include "projection.hpp"
#include <array>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class ElementOperations2D
{
  public:
    ElementOperations2D()
      : m_defaultProj(utils::Point(0, 0, 0), utils::Point(1, 0, 0), utils::Point(0, 1, 0))
    {
      set_projection(&m_defaultProj);
    }

    void set_projection(utils::impl::Projection* proj) { m_proj = proj; }

    utils::impl::Projection* get_projection() { return m_proj; }

    utils::Point compute_xi_coords(mesh::MeshEntityPtr el, const utils::Point& pt)
    {
      switch (el->get_type())
      {
        case stk::middle_mesh::mesh::MeshEntityType::Triangle: {
          return compute_tri_xi_coords(el, pt);
        }
        case stk::middle_mesh::mesh::MeshEntityType::Quad: {
          return compute_quad_xi_coords(el, pt);
        }
        default:
          throw std::invalid_argument("unsupported entity type");
      }
    }

    utils::Point compute_tri_xi_coords(mesh::MeshEntityPtr tri, const utils::Point& pt)
    {
      std::array<utils::Point, 3> pts;
      std::array<mesh::MeshEntityPtr, 3> verts;
      get_downward(tri, 0, verts.data());

      for (int i = 0; i < 3; ++i)
        pts[i] = m_proj->project_plane_coords(verts[i]->get_point_orig(0));

      auto ptProj = m_proj->project_plane_coords(pt);

      return compute_tri_xi_coords(pts, ptProj);
    }

    utils::Point compute_tri_xi_coords(const std::array<utils::Point, 3>& triVerts, const utils::Point& pt)
    {
      auto& n1 = triVerts[0];
      auto& n2 = triVerts[1];
      auto& n3 = triVerts[2];

      utils::impl::Mat2x2<double> a = {n2.get_x() - n1.get_x(), n3.get_x() - n1.get_x(), n2.get_y() - n1.get_y(),
                                       n3.get_y() - n1.get_y()};

      double dx[2] = {pt.get_x() - n1.get_x(), pt.get_y() - n1.get_y()};
      double xi[2];
      matsolve2x2(a, xi, dx);

      return utils::Point(xi[0], xi[1]);
    }

    utils::Point compute_quad_xi_coords(mesh::MeshEntityPtr quad, const utils::Point& pt)
    {
      assert(quad->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad);
      std::array<mesh::MeshEntityPtr, 4> verts;
      get_downward(quad, 0, verts.data());

      std::array<utils::Point, 4> pts;
      for (int i = 0; i < 4; ++i)
        pts[i] = m_proj->project_plane_coords(verts[i]->get_point_orig(0));

      utils::Point ptProj = m_proj->project_plane_coords(pt);

      return compute_quad_xi_coords(pts, ptProj);
    }

    utils::Point compute_quad_xi_coords(const std::array<utils::Point, 4>& verts, const utils::Point& pt)
    {
      auto& n1 = verts[0];
      auto& n2 = verts[1];
      auto& n3 = verts[2];
      auto& n4 = verts[3];

      double lagXi[2], lagEta[2];

      auto func = [&](const double x[2], double rhs[2]) {
        auto ptI = mesh::compute_quad_coords_from_xi(x, n1, n2, n3, n4);
        rhs[0]   = ptI.x - pt.x;
        rhs[1]   = ptI.y - pt.y;
      };

      auto jac = [&](const double x[2], utils::impl::Mat2x2<double>& a) {
        // polynomial derivatives are constant, don't compute them
        mesh::compute_lagrange_vals(x[0], lagXi);
        mesh::compute_lagrange_vals(x[1], lagEta);

        a(0, 0) = (n2.x - n1.x) * lagEta[0] + (n3.x - n4.x) * lagEta[1];
        a(0, 1) = (n4.x - n1.x) * lagXi[0] + (n3.x - n2.x) * lagXi[1];
        a(1, 0) = (n2.y - n1.y) * lagEta[0] + (n3.y - n4.y) * lagEta[1];
        a(1, 1) = (n4.y - n1.y) * lagXi[0] + (n3.y - n2.y) * lagXi[1];
      };

      double xi[2] = {0, 0};
      // double xi[2] = {0.297497, -1.00811};
      utils::impl::Newton2 newton(1e-14);
      newton.solve(func, jac, xi);

      return utils::Point(xi[0], xi[1]);
    }

    utils::Point compute_edge_coords(mesh::MeshEntityPtr edge, const double xi)
    {
      assert(edge->get_type() == stk::middle_mesh::mesh::MeshEntityType::Edge);
      auto p1 = m_proj->compute_plane_coords(edge->get_down(0)->get_point_orig(0));
      auto p2 = m_proj->compute_plane_coords(edge->get_down(1)->get_point_orig(0));
      return mesh::compute_edge_coords(p1, p2, xi);
    }

    double compute_area(mesh::MeshEntityPtr el)
    {
      if (el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad)
        return compute_quad_area(el);
      else if (el->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle)
        return compute_tri_area(el);
      else
        throw std::invalid_argument("unsupported entity type");
    }

    double compute_quad_area(mesh::MeshEntityPtr quad)
    {
      assert(quad->get_type() == stk::middle_mesh::mesh::MeshEntityType::Quad);

      mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
      get_downward(quad, 0, verts);
      return compute_tri_area(m_proj->compute_plane_coords(verts[0]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[1]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[2]->get_point_orig(0))) +
             compute_tri_area(m_proj->compute_plane_coords(verts[0]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[2]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[3]->get_point_orig(0)));
    }

    double compute_tri_area(mesh::MeshEntityPtr tri)
    {
      assert(tri->get_type() == stk::middle_mesh::mesh::MeshEntityType::Triangle);
      mesh::MeshEntityPtr verts[mesh::MAX_DOWN];
      get_downward(tri, 0, verts);
      return compute_tri_area(m_proj->compute_plane_coords(verts[0]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[1]->get_point_orig(0)),
                              m_proj->compute_plane_coords(verts[2]->get_point_orig(0)));
    }

    double compute_tri_area(const utils::Point& pt1, const utils::Point& pt2, const utils::Point& pt3)
    {
      return (pt1.x * (pt2.y - pt3.y) + pt2.x * (pt3.y - pt1.y) + pt3.x * (pt1.y - pt2.y)) / 2;
    }

    double compute_area(std::shared_ptr<mesh::Mesh> mesh)
    {
      double area = 0;
      for (auto& el : mesh->get_elements())
        if (el)
          area += compute_area(el);

      return area;
    }

  private:
    utils::impl::Projection* m_proj = nullptr;
    utils::impl::Projection m_defaultProj;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif
