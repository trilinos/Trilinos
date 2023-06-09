#ifndef POINT_CLASSIFIER_NORMAL_INTERPOLATION_REVERSE_SOLVER_H
#define POINT_CLASSIFIER_NORMAL_INTERPOLATION_REVERSE_SOLVER_H

#include <stk_middle_mesh/change_of_basis.hpp>
#include <stk_middle_mesh/mat2x2.hpp>
#include <stk_middle_mesh/mesh_entity.hpp>
#include <stk_middle_mesh/newton_scalar.hpp>
#include <stk_middle_mesh/projection.hpp>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

// solves for gamma (the distance along the normal vector) by solving
// a rootfinding problem f(gamma, lambda(gamma)) = 0.  Note that
// the lambdas are treated as functions of gamma.
class ReverseClassificationGammaSolver
{
  public:
    ReverseClassificationGammaSolver(const std::array<utils::Point, 3>& verts, const utils::Point& pt,
                                     const std::array<utils::Point, 4>& normals, double convergenceTol = 1e-13,
                                     int itermax = 100000)
      : m_x0(verts[0])
      , m_dx1(verts[1] - verts[0])
      , m_dx2(verts[2] - verts[0])
      , m_d0(normals[0])
      , m_dn1(normals[1] - normals[0])
      , m_dn2(normals[2] - normals[0])
      , m_pt(pt)
      , m_newton(convergenceTol, itermax)
    {
      std::array<utils::Point, 3> basis;
      basis[0] = m_dx1 / std::sqrt(dot(m_dx1, m_dx1));
      basis[1] = m_dx2 - dot(m_dx2, basis[0]) * basis[0];
      basis[1] = basis[1] / std::sqrt(dot(basis[1], basis[1]));
      basis[2] = cross(basis[0], basis[1]);

      /*
      std::array<utils::Point, 3> basis;
      utils::Point dx1 = verts[1] - verts[0];
      utils::Point normal = (normals[0] + normals[1] + normals[1])/3;
      normal = normal / std::sqrt(dot(normal, normal));
      basis[0] = normal;
      basis[1] = dx1 - dot(dx1, basis[0]) * basis[0];
      basis[2] = cross(basis[0], basis[1]);
      */

      // the matrix for computing the lambdas can be singular
      // if all the points are in the xz plane.  Do a change
      // of basis so it is full rank.
      utils::impl::ChangeOfBasis q(basis);
      m_x0  = q.project_forward(m_x0);
      m_dx1 = q.project_forward(m_dx1);
      m_dx2 = q.project_forward(m_dx2);
      m_d0  = q.project_forward(m_d0);
      m_dn1 = q.project_forward(m_dn1);
      m_dn2 = q.project_forward(m_dn2);
      m_pt  = q.project_forward(m_pt);
    }

    // returns an array [lambda1, lambda2, gamma]
    std::array<double, 3> solve(bool& didConverge);

  private:
    std::array<double, 2> solve_for_lambda(double gamma);

    void solve_for_lambdad_gamma(double gamma, std::array<double, 2>& lambdas, std::array<double, 2>& lambdasDot);

    double compute_f_of_gamma(double gamma);

    double compute_df_dgamma(double gamma);

    utils::Point m_x0;
    utils::Point m_dx1;
    utils::Point m_dx2;

    utils::Point m_d0;
    utils::Point m_dn1;
    utils::Point m_dn2;
    utils::Point m_pt;
    utils::impl::NewtonScalar m_newton;
};

} // namespace impl
} // namespace predicates

} // namespace middle_mesh
} // namespace stk
#endif
